import pandas as pd
import numpy as np
import dask.dataframe as dd
import os
from functools import partial
from ..utils import assign_to_use_column
import logging
dask_logger = logging.getLogger("dask")
import time
from typing import Union


def _scaffold_breaker(group, slop):
    group = group.sort_values("breakpoint", ascending=True)
    breakpoints = group["breakpoint"].values
    previous_iteration = group["previous_iteration"].values
    previous_iteration = np.concatenate([np.repeat(previous_iteration, 2), [previous_iteration[-1]]])
    assert breakpoints.shape[0] == np.unique(breakpoints).shape[0], group[["scaffold_index", "breakpoint"]]
    length = group["length"].unique()[0]
    scaffold_starts = np.vstack([breakpoints - slop + 1, breakpoints + slop + 1]).T.flatten()
    scaffold_starts = np.concatenate([[1], scaffold_starts])
    scaffold_ends = np.vstack([breakpoints - slop, breakpoints + slop]).T.flatten()
    scaffold_ends = np.concatenate([scaffold_ends, [length + scaffold_starts[0] - 1]])
    ends = scaffold_ends - scaffold_starts + 1
    final = np.vstack([
        np.repeat(1, scaffold_ends.shape[0]).T,  # start
        ends.T,  # end = length
        scaffold_starts.T,
        scaffold_ends.T,
        previous_iteration.T
    ]).T
    final = np.unique(final, axis=0)
    return final


def _fai_checker(fai):
    """This function will check that the derived scaffolds are coherent, ie, that each broken scaffolds will have final
    segments that *completely cover* the original scaffold, without any internal overlapping."""
    dev = fai.query("to_use == True & derived_from_split==True").compute()
    orig = fai.loc[dev.orig_scaffold_index.unique()].compute()
    imputed = (dev.groupby("orig_scaffold_index")["orig_end"].sum() - dev.groupby("orig_scaffold_index")[
        "orig_start"].sum() + dev.groupby("orig_scaffold_index").size()).to_frame("imputed")
    merged = pd.merge(imputed, orig[["length"]], left_index=True, right_index=True)
    assert merged.shape[0] == imputed.shape[0] == orig.shape[0]
    return (merged.query("imputed != length").shape[0] == 0)


def calculate_broken_scaffolds(breaks: pd.DataFrame, fai: str, save_dir: Union[str, None], slop: float,
                               cores=1) -> (dd.DataFrame, Union[str, None], np.ndarray):

    """
    This function will take the position of the breaks found by `find_10x_breaks` and
    determine how to split up the scaffolds, *while keeping track of the origin of each single
    resulting scaffold.*
    The slop parameter determines how much to keep around the breaking point.
    """

    if isinstance(fai, str):
        fai = dd.read_parquet(fai, infer_divisions=True)
    else:
        assert isinstance(fai, dd.DataFrame)

    original_types = dict(fai.dtypes)
    if "derived_from_split" not in fai.columns:  # Not initialised
        fai["derived_from_split"] = False
        fai["previous_iteration"] = fai.index
    original_indices = fai.index.values.compute()

    broken = breaks.copy()
    assert fai.index.dtype == broken.index.dtype == int
    broken = dd.merge(fai, broken.drop("length", axis=1, errors="ignore"), on="scaffold_index", how="right").compute()

    # Remove spurious breaks
    if broken.index.name == "scaffold_index":
        broken = broken.reset_index(drop=False)
    broken = broken.assign(idx=np.arange(1, broken.shape[0] + 1, dtype=int)).set_index("idx")
    broken = broken.sort_values(["scaffold_index", "breakpoint"], ascending=[True, True])
    assert broken.shape[0] == broken[["scaffold_index", "breakpoint"]].drop_duplicates().shape[0]
    assert "scaffold_index" in broken.columns
    grouped = broken.groupby("scaffold_index")
    check = (grouped["breakpoint"].shift(0) - grouped["breakpoint"].shift(1) <= 3 * slop)
    while check.any():
        broken = broken.loc[broken.index.values[~check]]
        grouped = broken.groupby("scaffold_index")
        check = (grouped["breakpoint"].shift(0) - grouped["breakpoint"].shift(1) <= 2 * slop)
    assert broken["scaffold_index"].unique().shape[0] == breaks.index.unique().shape[0]

    # Now we are guaranteed that all the breaks are at the correct distance.
    # Next up: ensure that none of these breaks happened in non-original scaffolds.
    # If they do, we need to change the coordinates.
    broken.loc[:, "original_breakpoint"] = broken.loc[:, ["breakpoint", "orig_start"]].sum(axis=1) - 1
    assert broken.shape[0] == broken[["orig_scaffold_index", "original_breakpoint"]].drop_duplicates().shape[0]
    broken["previous_iteration"] = broken["scaffold_index"]

    # Now go back to the fai. We need to keep memory of the index we start with.
    # Also we need to *remove* from the FAI the scaffolds we are breaking further, I think?
    right = fai[["scaffold"]].rename(columns={"scaffold": "orig_scaffold"})
    broken = dd.merge(broken.reset_index(drop=False), right, left_on="orig_scaffold_index", right_on="scaffold_index")
    broken = broken.compute().drop("idx", axis=1).set_index("scaffold_index")
    broken = broken.drop("scaffold", axis=1).rename(columns={"orig_scaffold": "scaffold"})
    assert "scaffold" in broken.columns, broken.columns
    # broken["derived_from_split"] = False
    sloppy = partial(_scaffold_breaker, slop=slop)
    dask_logger.debug("%s calculate_broken_scaffolds -  Starting scaffold breaking", time.ctime())
    maxid = fai.index.values.max()
    broken = dd.from_pandas(broken, npartitions=cores)
    assert "orig_scaffold_index" in broken.columns
    grouped_breaks = broken.groupby("scaffold_index")
    new_scaffolds = np.vstack(grouped_breaks.apply(sloppy, meta=np.int).compute().values)
    # _broken is a numpy array list of the form
    # start, end, orig_start, orig_ends, orig_scaffold_index
    new_scaffolds = pd.DataFrame().assign(
        start=new_scaffolds[:, 0], end=new_scaffolds[:, 1], length=new_scaffolds[:, 1],
        scaffold_start=new_scaffolds[:, 2], scaffold_end=new_scaffolds[:, 3], previous_iteration=new_scaffolds[:, 4],
        derived_from_split=True, scaffold_index=np.arange(maxid + 1, maxid + new_scaffolds.shape[0] + 1, dtype=np.int)
    )
    assert new_scaffolds.query("length < @slop").shape[0] == 0
    new_scaffolds = new_scaffolds.astype(dict((key, int) for key in new_scaffolds.columns
                                              if key != "derived_from_split"))
    orig_shape = new_scaffolds.shape[0]
    assert "scaffold_index" in new_scaffolds.columns
    new_scaffolds = dd.merge(new_scaffolds,
                       broken[["orig_scaffold_index", "previous_iteration", "scaffold", "orig_start"]
                       ].reset_index(drop=True).drop_duplicates(),
                       how="left",
                       on="previous_iteration")
    new_scaffolds["orig_start"] = new_scaffolds["scaffold_start"] + new_scaffolds["orig_start"] - 1
    new_scaffolds["orig_end"] = new_scaffolds["end"] + new_scaffolds["orig_start"] - 1
    new_scaffolds["length_check"] = new_scaffolds.eval("orig_end - orig_start + 1")
    assert new_scaffolds.query("length_check != length").shape[0].compute() == 0
    assert new_scaffolds.query("length_check < @slop", local_dict={"slop": slop}).shape[0].compute() == 0
    new_scaffolds = new_scaffolds.drop("length_check", axis=1)
    new_shape = new_scaffolds.shape[0].compute()
    assert new_shape == orig_shape, (orig_shape, new_shape)
    assert "scaffold_index" in new_scaffolds.columns, (new_scaffolds.columns, broken.columns)
    new_scaffolds = new_scaffolds.reset_index(drop=True).set_index("scaffold_index")
    assert (np.sort(new_scaffolds["previous_iteration"].unique().values.compute()) ==
            np.sort(breaks.index.unique().values)).all() == True
    assert "previous_iteration" in new_scaffolds.columns, new_scaffolds.columns
    assert "orig_scaffold_index" in new_scaffolds.columns, new_scaffolds.columns
    new_scaffolds = new_scaffolds[fai.columns.intersection(new_scaffolds.columns).values]
    new_scaffolds["scaffold"] = (new_scaffolds["scaffold"].astype(str) + ":"
                                 + new_scaffolds["orig_start"].astype(str) +
                                 "-" + new_scaffolds["orig_end"].astype(str))

    dask_logger.debug("%s calculate_broken_scaffolds -  Finished scaffold breaking", time.ctime())
    dask_logger.debug("%s calculate_broken_scaffolds -  Merging with the original FAI", time.ctime())
    try:
        nparts = fai.npartitions
        fai = dd.concat([fai.reset_index(drop=False),
                         new_scaffolds.reset_index(drop=False)]).astype(
            {"scaffold_index": np.int}).set_index("scaffold_index")
        fai = fai.repartition(npartitions=nparts)
    except NotImplementedError:
        print(fai.head())
        print("###\n###")
        print(broken.head())
        raise
    assert fai.index.name == "scaffold_index", fai.head()
    assert "orig_scaffold_index" in fai.columns, fai.columns
    fai = fai.astype(original_types)
    fai = fai.astype({"previous_iteration": int})
    dask_logger.debug("%s calculate_broken_scaffolds -  Finished, returning the FAI", time.ctime())
    fai = assign_to_use_column(fai)
    if save_dir is not None:
        fai_name = os.path.join(save_dir, "fai")
        dd.to_parquet(fai, fai_name, compression="gzip", compute=True, engine="pyarrow", schema="infer")
    else:
        fai_name = None
    # new_indices = fai.index.compute().difference(original_indices).values
    # Final check
    valid = _fai_checker(fai)
    if valid is False:
        dask_logger.critical("%s Bungled FAI - we have created overlapping fragments. Check %s",
                             time.ctime(), fai_name)
        import sys
        sys.exit(1)

    return fai, fai_name
