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
    assert breakpoints.shape[0] == np.unique(breakpoints).shape[0], group[["scaffold_index", "breakpoint"]]
    # starts, ends = [1], [breakpoints[0] - slop]
    length = group["length"].unique()[0]
    orig_starts = np.vstack([breakpoints - slop + 1, breakpoints + slop + 1]).T.flatten()
    orig_starts = np.concatenate([[1], orig_starts])
    orig_ends = np.vstack([breakpoints - slop, breakpoints + slop]).T.flatten()
    orig_ends = np.concatenate([orig_ends, [length]])
    lengths = ends = orig_ends - orig_starts + 1
    final = np.vstack([
        np.repeat(1, orig_ends.shape[0]).T,  # start
        ends.T,  # end = length
        orig_starts.T,
        orig_ends.T,
        np.repeat(int(np.unique(group.index.values)[0]), orig_ends.shape[0]).T
    ]).T
    # df["scaffold"] = df["scaffold"] + ":" + df["orig_start"].astype(str) + "-" + df["orig_end"].astype(str)

    return final


def calculate_broken_scaffolds(breaks: pd.DataFrame, fai: str, save_dir: str, slop: float) -> dict:

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
    fai["derived_from_split"] = False

    broken = breaks.copy()
    broken = dd.merge(fai, broken.drop("length", axis=1, errors="ignore"),
                      on="scaffold_index", how="right")
    if broken.index.name == "scaffold_index":
        broken = broken.reset_index(drop=False)
    broken["idx"] = da.from_array(
        np.arange(1, broken.shape[0].compute() + 1),
        chunks=tuple(broken.map_partitions(len).compute().values.tolist()))

    broken = broken.set_index("idx", sorted=True)
    assert broken.shape[0].compute() == broken[
        ["scaffold_index", "breakpoint"]].drop_duplicates().shape[0].compute()
    assert "scaffold_index" in broken.columns

    # Broken structure:
    # scaffold_index
    # Index(['scaffold', 'length', 'orig_scaffold_index', 'start', 'orig_start',
    #        'end', 'orig_end', 'derived_from_split', 'breakpoint', 'n', 'd', 'nbin',
    #        'mn', 'r', 'b'],
    #       dtype='object')

    # First step: find out if any breakpoint has to be removed.
    grouped = broken.groupby("scaffold_index")
    # Check whether any two breakpoints are at less than "slop" distance.
    check = (grouped["breakpoint"].apply(
        lambda group: group.shift(0) - group.shift(1), meta=np.float) <= 2 * slop).compute()
    while check.any():
        broken = broken.loc[broken.index.compute().values[~check]]
        grouped = broken.groupby("scaffold_index")
        check = (grouped["breakpoint"].apply(
            lambda group: group.shift(0) - group.shift(1), meta=np.float) <= 2 * slop).compute()
    # Now we are guaranteed that all the breaks are at the correct distance.
    # Next up: ensure that none of these breaks happened in non-original scaffolds.
    # If they do, we need to change the coordinates.
    assert broken.shape[0].compute() == broken[
        ["scaffold_index", "breakpoint"]].drop_duplicates().shape[0].compute()
    bait = (broken["derived_from_split"] == True)
    broken["breakpoint"] = broken["breakpoint"].mask(bait,
                                                     broken["breakpoint"] + broken["orig_start"] - 1)

    # Now go back to the fai. We need to keep memory of the index we start with.
    # Also we need to *remove* from the FAI the scaffolds we are breaking further, I think?
    right = fai[["scaffold"]].rename(columns={"scaffold": "orig_scaffold"})
    broken = dd.merge(broken.reset_index(drop=False), right,
                      left_on="orig_scaffold_index", right_index=True).persist().set_index("scaffold_index")
    broken["scaffold"] = broken["scaffold"].mask(bait, broken["orig_scaffold"])
    broken = broken.drop("orig_scaffold", axis=1)
    broken["scaffold_index"] = broken.index
    broken.index = broken.index.rename("old_scaffold_index")
    broken["scaffold_index"] = broken["scaffold_index"].mask(bait, broken["orig_scaffold_index"])
    broken = broken.persist()
    broken = broken.drop_duplicates(subset=["scaffold_index", "breakpoint"], keep="first").persist()
    assert broken.shape[0].compute() == broken[["scaffold_index", "breakpoint"]].drop_duplicates().shape[0].compute()
    # broken["derived_from_split"] = False
    sloppy = partial(_scaffold_breaker, slop=slop)
    dask_logger.warning("%s Starting scaffold breaking", time.ctime())
    maxid = fai.index.values.max()
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
        ends=_broken[:, 1],
        orig_start=_broken[:, 2],
        orig_end=_broken[:, 3],
        orig_scaffold_index=_broken[:, 4],
        derived_from_split=True,
        scaffold_index=np.arange(maxid + 1, maxid + _broken.shape[0] + 1, dtype=np.int)
    )
    assert "scaffold_index" in _broken.columns
    _broken = dd.merge(_broken, broken,
                       how="left",
                       suffixes=("", "_old"),
                       on="orig_scaffold_index").persist()
    assert "scaffold_index" in _broken.columns, (_broken.columns, broken.columns)
    _broken = _broken.reset_index(drop=True).set_index("scaffold_index")
    _broken = _broken[fai.columns]
    _broken["scaffold"] = (_broken["scaffold"] + ":" + _broken["orig_start"].astype(str) +
                           "-" + _broken["orig_end"].astype(str))

    dask_logger.warning("%s Finished scaffold breaking", time.ctime())
    # _broken = dd.from_pandas(broken, npartitions=100).groupby("scaffold_index").apply(sloppy).compute()
    # assert (_broken["start"].compute() == 1).all()
    # assert _broken["derived_from_split"].compute().all()

    dask_logger.warning("%s Merging with the original FAI", time.ctime())
    try:
        fai = dd.concat([fai.reset_index(drop=False),
                         _broken.reset_index(drop=False)]).astype(
            {"scaffold_index": np.int}).set_index("scaffold_index")
    except NotImplementedError:
        print(fai.head())
        print("###\n###")
        print(broken.head())
        raise
    # assert fai[fai["derived_from_split"] == True].shape[0].compute() >= _broken.shape[0].compute()
    assert fai.index.name == "scaffold_index", fai.head()
    fai = fai.persist()
    dask_logger.warning("%s Finished, returning the FAI", time.ctime())
    fai_name = os.path.join(save_dir, "fai")
    dd.to_parquet(fai, fai_name, compression="gzip", compute=True, engine="pyarrow")
    return {"fai": fai_name}
