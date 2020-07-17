import pandas as pd
import numpy as np
import dask.dataframe as dd
import os
from functools import partial
import dask.array as da
import logging
dask_logger = logging.getLogger("dask")
import time


def _scaffold_breaker(group, slop):
    breakpoints = group["breakpoint"].values
    previous_iteration = group["previous_iteration"].values
    previous_iteration = np.concatenate([np.repeat(previous_iteration, 2), [previous_iteration[-1]]])
    assert breakpoints.shape[0] == np.unique(breakpoints).shape[0], group[["scaffold_index", "breakpoint"]]
    # starts, ends = [1], [breakpoints[0] - slop]
    length = group["length"].unique()[0]
    orig_starts = np.vstack([breakpoints - slop + 1, breakpoints + slop + 1]).T.flatten()
    orig_starts = np.concatenate([[group["orig_start"].min()], orig_starts])
    orig_ends = np.vstack([breakpoints - slop, breakpoints + slop]).T.flatten()
    orig_ends = np.concatenate([orig_ends, [length + orig_starts[0] - 1]])
    ends = orig_ends - orig_starts + 1
    final = np.vstack([
        np.repeat(1, orig_ends.shape[0]).T,  # start
        ends.T,  # end = length
        orig_starts.T,
        orig_ends.T,
        previous_iteration.T
    ]).T
    final = np.unique(final, axis=0)
    return final


def calculate_broken_scaffolds(breaks: pd.DataFrame, fai: str, save_dir: str, slop: float, cores=1) -> dict:

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

    if "derived_from_split" not in fai.columns:  # Not initialised
        fai["derived_from_split"] = False
        fai["previous_iteration"] = fai.index

    broken = breaks.copy()
    assert fai.index.dtype == broken.index.dtype == int
    broken = dd.merge(fai, broken.drop("length", axis=1, errors="ignore"), on="scaffold_index", how="right")
    broken = broken.compute()
    if broken.index.name == "scaffold_index":
        broken = broken.reset_index(drop=False)
    broken["idx"] = np.arange(1, broken.shape[0] + 1)
    broken = broken.set_index("idx")
    assert broken.shape[0] == broken[["scaffold_index", "breakpoint"]].drop_duplicates().shape[0]
    assert "scaffold_index" in broken.columns

    # First step: find out if any breakpoint has to be removed.
    grouped = broken.groupby("scaffold_index")
    # Check whether any two breakpoints are at less than "slop" distance.
    check = (grouped["breakpoint"].shift(0) - grouped["breakpoint"].shift(1) <=2 * slop)
    while check.any():
        broken = broken.loc[broken.index.values[~check]]
        grouped = broken.groupby("scaffold_index")
        check = (grouped["breakpoint"].shift(0) - grouped["breakpoint"].shift(1) <= 2 * slop)

    # Now we are guaranteed that all the breaks are at the correct distance.
    # Next up: ensure that none of these breaks happened in non-original scaffolds.
    # If they do, we need to change the coordinates.
    assert broken.shape[0] == broken[["scaffold_index", "breakpoint"]].drop_duplicates().shape[0]
    broken["previous_iteration"] = broken["scaffold_index"]
    bait = (broken["derived_from_split"] == True)
    broken["breakpoint"] = broken["breakpoint"].mask(bait, broken["breakpoint"] + broken["orig_start"] - 1)

    # Now go back to the fai. We need to keep memory of the index we start with.
    # Also we need to *remove* from the FAI the scaffolds we are breaking further, I think?
    right = fai[["scaffold"]].rename(columns={"scaffold": "orig_scaffold"})
    broken = dd.merge(broken.reset_index(drop=False), right, left_on="orig_scaffold_index",
                      right_on="scaffold_index")
    broken = broken.drop("idx", axis=1).set_index("scaffold_index").compute()
    broken["scaffold"] = broken["scaffold"].mask(bait, broken["orig_scaffold"])

    broken = broken.drop("orig_scaffold", axis=1)
    broken["scaffold_index"] = broken.index
    broken.index = broken.index.rename("old_scaffold_index")
    broken["scaffold_index"] = broken["scaffold_index"].mask(bait, broken["orig_scaffold_index"])
    broken = broken.drop_duplicates(subset=["scaffold_index", "breakpoint"], keep="first")

    # broken["derived_from_split"] = False
    sloppy = partial(_scaffold_breaker, slop=slop)
    dask_logger.warning("%s calculate_broken_scaffolds -  Starting scaffold breaking", time.ctime())
    maxid = fai.index.values.max()
    broken = dd.from_pandas(broken, npartitions=cores)
    _broken = broken.reset_index(drop=True).set_index(
        "scaffold_index").groupby("scaffold_index").apply(sloppy, meta=np.int).compute().values
    _broken = np.vstack(_broken)
    # _broken is a numpy array list of the form
    # start, end, orig_start, orig_ends, orig_scaffold_index
    if len(_broken.shape) != 2 or _broken.shape[1] != 5:
        dask_logger.critical("Wrong number of columns in _broken")
        dask_logger.critical(_broken)
        raise AssertionError
    assert _broken.shape[1] == 5
    _broken = pd.DataFrame().assign(
        start=_broken[:, 0],
        end=_broken[:, 1],
        length=_broken[:, 1],
        orig_start=_broken[:, 2],
        orig_end=_broken[:, 3],
        previous_iteration=_broken[:, 4],
        # previous_iteration=_broken[:, 5],
        derived_from_split=True,
        scaffold_index=np.arange(maxid + 1, maxid + _broken.shape[0] + 1, dtype=np.int)
    )
    orig_shape = _broken.shape[0]
    assert "scaffold_index" in _broken.columns
    _broken = dd.merge(_broken,
                       broken[["orig_scaffold_index", "previous_iteration", "scaffold"]
                       ].reset_index(drop=True).drop_duplicates(),
                       how="left",
                       on="previous_iteration")
    new_shape = _broken.shape[0].compute()
    assert new_shape == orig_shape, (orig_shape, new_shape)
    assert "scaffold_index" in _broken.columns, (_broken.columns, broken.columns)
    _broken = _broken.reset_index(drop=True).set_index("scaffold_index")
    _broken = _broken[fai.columns]
    broken = broken.astype({"orig_start": int, "orig_end": int})
    _broken["scaffold"] = (_broken["scaffold"].astype(str) + ":" + _broken["orig_start"].astype(str) +
                           "-" + _broken["orig_end"].astype(str))

    dask_logger.warning("%s calculate_broken_scaffolds -  Finished scaffold breaking", time.ctime())
    dask_logger.warning("%s calculate_broken_scaffolds -  Merging with the original FAI", time.ctime())
    try:
        nparts = fai.npartitions
        fai = dd.concat([fai.reset_index(drop=False),
                         _broken.reset_index(drop=False)]).astype(
            {"scaffold_index": np.int}).set_index("scaffold_index")
        fai = fai.repartition(npartitions=nparts)
    except NotImplementedError:
        print(fai.head())
        print("###\n###")
        print(broken.head())
        raise
    assert fai.index.name == "scaffold_index", fai.head()
    dask_logger.warning("%s calculate_broken_scaffolds -  Finished, returning the FAI", time.ctime())
    fai_name = os.path.join(save_dir, "fai")
    dd.to_parquet(fai, fai_name, compression="gzip", compute=True, engine="pyarrow")
    return {"fai": fai_name}
