import pandas as pd
import numpy as np
import dask.dataframe as dd
import os
from functools import partial


def _calculate_coordinates(broken, slop, maxid):
    ids = (maxid + np.arange(1, 3 * (broken.shape[0] + 1) - 2, dtype=np.int)).reshape(broken.shape[0], 3)
    # add scaffold_index1, 2, 3
    for index in range(3):
        broken.loc[:, "scaffold_index" + str(index + 1)] = ids[:, index]
    broken.loc[:, "start1"] = 1
    #   br[, end1 := pmax(0, br - slop - 1)]
    broken.loc[:, "end1"] = np.maximum(0, (broken["breakpoint"] - slop - 1).astype(int), dtype=np.int)
    #   br[, start2 := pmax(1, br - slop)]
    broken.loc[:, "start2"] = np.maximum(1, (broken["breakpoint"] - slop).astype(int), dtype=np.int)
    #   br[, end2 := pmin(br + slop - 1, length)]
    broken.loc[:, "end2"] = np.minimum(broken["length"], (broken["breakpoint"] + slop - 1).astype(int), dtype=np.int)
    #   br[, start3 := pmin(length + 1, br + slop)]
    broken.loc[:, "start3"] = np.minimum((broken["length"] + 1).astype(int),
                                         (broken["breakpoint"] + slop).astype(int), dtype=np.int)
    broken.loc[:, "end3"] = broken["length"].astype(np.int)
    return broken


def _create_children_dataframes(broken):
    dfs = []
    for index in range(1, 4):
        index = str(index)
        lkey, skey, ekey = "length" + index, "start" + index, "end" + index
        broken.loc[:, lkey] = 1 + broken[ekey] - broken[skey]
        # Assign a new name WHICH KEEPS TRACK OF WHERE WE ARE COMING FROM
        broken.loc[:, "scaffold" + index] = (
                broken["scaffold"] + ":" + broken[skey].astype(str) + "-" + broken[ekey].astype(str))
        # FAI has the following keys:
        # scaffold_index scaffold  length  orig_scaffold_index  start  orig_start  end  orig_end derived_from_split
        df = pd.DataFrame().assign(
            scaffold_index=broken["scaffold_index" + index],
            scaffold=broken["scaffold" + index],
            length=broken[lkey],
            orig_scaffold_index=broken["orig_scaffold_index"],
            start=1,
            end=broken[lkey],
            orig_start=broken["orig_start"] + broken[skey] - 1,
            orig_end=broken["orig_start"] + broken[skey] + broken[lkey] - 2,
            derived_from_split=True
        )
        for col in df.columns:
            if col in ["scaffold", "derived_from_split"]:
                continue
            df.loc[:, col] = pd.to_numeric(df[col], downcast="signed")

        df = df.loc[df["length"] > 0, :].copy()
        dfs.append(df)
    return pd.concat(dfs).reset_index(drop=True).set_index("scaffold_index")


def _scaffold_breaker(group, slop):
    breakpoints = group["breakpoint"].values
    starts, ends = [], []
    cstart = 1
    length = group["length"].unique()[0]
    for bp in breakpoints:
        starts.extend([cstart, bp - slop])
        cend = np.min([bp + slop, length])
        ends.extend([bp - slop + 1, np.min([bp + slop, length])])
        cstart = cend + 1
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
                      on="scaffold_index", how="right").compute()
    if broken.index.name == "scaffold_index":
        broken = broken.reset_index(drop=False)
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
    check = (grouped["breakpoint"].shift(0) - grouped["breakpoint"].shift(1)) < 2 * slop
    while check.any():
        broken = broken.loc[~check]
        grouped = broken.groupby("scaffold_index")
        check = (grouped["breakpoint"].shift(0) - grouped["breakpoint"].shift(1)) < 2 * slop
    # Now we are guaranteed that all the breaks are at the correct distance.
    # Next up: ensure that none of these breaks happened in non-original scaffolds.
    # If they do, we need to change the coordinates.
    bait = (broken["derived_from_split"] == True)
    broken.loc[bait, "breakpoint"] = broken["orig_start"] + broken["breakpoint"]
    # Now go back to the fai. We need to keep memory of the index we start with.
    # Also we need to *remove* from the FAI the scaffolds we are breaking further, I think?
    broken.loc[bait, "scaffold"] = fai.loc[broken.loc[bait, "orig_scaffold_index"].values,
                                           "scaffold"].compute().to_numpy()
    # to_remove = broken.loc[bait, "scaffold_index"].to_numpy()
    # fai = fai.drop(to_remove, axis=0)
    broken.loc[bait, "scaffold_index"] = broken.loc[bait, "orig_scaffold_index"]
    broken["derived_from_split"] = False
    sloppy = partial(_scaffold_breaker, slop=slop)
    _broken = dd.from_pandas(broken, npartitions=100).groupby("scaffold_index").apply(sloppy).compute()
    assert (_broken["start"] == 1).all()
    maxid = fai.index.values.max()
    _broken.loc[:, "scaffold_index"] = np.arange(maxid + 1, maxid + _broken.shape[0] + 1,
                                                 dtype=np.int)
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
