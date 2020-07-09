import pandas as pd
import numpy as np
import dask.dataframe as dd
import os
from functools import partial
import dask.array as da


def _scaffold_breaker(group, slop):
    breakpoints = group["breakpoint"].values
    assert breakpoints.shape[0] == np.unique(breakpoints).shape[0], group
    # starts, ends = [1], [breakpoints[0] - slop]
    starts, ends = [], []
    cstart = 1
    length = group["length"].unique()[0]
    mlength = len(breakpoints) - 1
    for index, bp in enumerate(breakpoints):
        starts.extend([cstart, bp - slop])
        ends.extend([bp - slop - 1, bp + slop])
        cstart = bp + slop + 1
        if index == mlength:
            starts.append(cstart)
            ends.append(length)

    df = pd.DataFrame().assign(
        start=1,
        orig_start=starts,
        orig_end=ends,
        derived_from_split=True,
        scaffold=group["scaffold"].unique()[0],
        orig_scaffold_index=group["scaffold_index"].unique()[0]
    )
    df = df.eval("end=orig_end - orig_start + 1").eval("length=end")
    df["scaffold"] = df["scaffold"] + ":" + df["orig_start"].astype(str) + "-" + df["orig_end"].astype(str)
    df["start"] = 1
    assert (df["start"] == 1).all()
    assert df.index.name is None
    assert "scaffold_index" not in df.columns
    return df


def calculate_broken_scaffolds(breaks: pd.DataFrame, fai: str, save_dir: str,
                               slop: float) -> dict:

    """
    This function will take the position of the breaks found by `find_10x_breaks` and
    determine how to split up the scaffolds, *while keeping track of the origin of each single
    resulting scaffold.*
    The slop parameter determines how much to keep around the breaking point.
    """

    fai = dd.read_parquet(fai, infer_divisions=True)
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
        lambda group: group.shift(0) - group.shift(1), meta=np.float) < 2 * slop).compute()
    while check.any():
        broken = broken.loc[broken.index.compute().values[~check]]
        grouped = broken.groupby("scaffold_index")
        check = (grouped["breakpoint"].apply(
            lambda group: group.shift(0) - group.shift(1), meta=np.float) < 2 * slop).compute()
    # Now we are guaranteed that all the breaks are at the correct distance.
    # Next up: ensure that none of these breaks happened in non-original scaffolds.
    # If they do, we need to change the coordinates.
    assert broken.shape[0].compute() == broken[
        ["scaffold_index", "breakpoint"]].drop_duplicates().shape[0].compute()
    bait = (broken["derived_from_split"] == True)
    broken["breakpoint"] = broken["breakpoint"].mask(bait,
                                                     broken["breakpoint"] - broken["orig_start"])

    # Now go back to the fai. We need to keep memory of the index we start with.
    # Also we need to *remove* from the FAI the scaffolds we are breaking further, I think?
    right = fai[["scaffold"]].rename(columns={"scaffold": "orig_scaffold"})
    broken = dd.merge(broken.reset_index(drop=False), right,
                      left_on="orig_scaffold_index", right_index=True).persist()
    broken["scaffold"] = broken["scaffold"].mask(bait,
                                                 broken["orig_scaffold"])
    broken = broken.drop("orig_scaffold", axis=1)
    assert broken.shape[0] == broken[["scaffold_index", "breakpoint"]].drop_duplicates().shape[0]
    # to_remove = broken.loc[bait, "scaffold_index"].to_numpy()
    # fai = fai.drop(to_remove, axis=0)
    broken["scaffold_index"] = broken.mask(bait, broken["orig_scaffold_index"])
    assert broken.shape[0] == broken[["scaffold_index", "breakpoint"]].drop_duplicates().shape[0]
    broken["derived_from_split"] = False
    sloppy = partial(_scaffold_breaker, slop=slop)
    _broken = broken.groupby("scaffold_index").apply(sloppy)
    # _broken = dd.from_pandas(broken, npartitions=100).groupby("scaffold_index").apply(sloppy).compute()
    assert (_broken["start"] == 1).all()
    maxid = fai.index.values.max()
    _broken["scaffold_index"] = np.arange(maxid + 1, maxid + _broken.shape[0] + 1, dtype=np.int)
    assert _broken.index.name is None
    _broken = _broken.reset_index(drop=True).set_index("scaffold_index")[fai.columns]
    assert _broken["derived_from_split"].all()
    assert (_broken["start"] == 1).all()
    try:
        fai = dd.concat([fai.reset_index(drop=False),
                         _broken.reset_index(drop=False)]).set_index("scaffold_index")
    except NotImplementedError:
        print(fai.head())
        print("###\n###")
        print(broken.head())
        raise
    assert fai[fai["derived_from_split"] == True].shape[0].compute() >= _broken.shape[0]
    assert fai.index.name == "scaffold_index", fai.head()
    fai_name = os.path.join(save_dir, "fai")
    import joblib
    joblib.dump(fai.compute(), "fai_test.pkl", compress=("gzip", 6))
    dd.to_parquet(fai, fai_name, compute=True, compression="gzip", engine="pyarrow")
    return {"fai": fai_name}
