import pandas as pd
import numpy as np
import dask.dataframe as dd
import os


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
        # br[, .(scaffold=scaffold1, length=length1, orig_scaffold, orig_start=orig_start+start1-1,
        #        orig_end=orig_start+end1-1, old_scaffold)],
        # br[, .(scaffold=scaffold2, length=length2, orig_scaffold, orig_start=orig_start+start2-1,
        #        orig_end=orig_start+end2-1, old_scaffold)],
        # br[, .(scaffold=scaffold3, length=length3, orig_scaffold, orig_start=orig_start+start3-1,
        #        orig_end=orig_start+end3-1, old_scaffold)],

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
    return pd.concat(dfs).reset_index(drop=True)


def calculate_broken_scaffolds(breaks: pd.DataFrame, fai: dd.DataFrame, save_dir: str,
                               slop: float) -> dict:

    """
    This function will take the position of the breaks found by `find_10x_breaks` and
    determine how to split up the scaffolds, *while keeping track of the origin of each single
    resulting scaffold.*
    The slop parameter determines how much to keep around the breaking point.
    """

    fai = dd.read_parquet(fai)
    fai["derived_from_split"] = False

    broken = breaks.copy()
    broken = dd.merge(fai, broken.drop("length", axis=1, errors="ignore"),
                      on="scaffold_index", how="right").compute()
    cycle = 0
    while broken.shape[0] > 0:
        # breaks: scaffold_index  length   break   n       d  nbin          mn         r         b
        # fai: scaffold_index scaffold  length  orig_scaffold_index  start  orig_start     end  orig_end
        initial_shape = fai.shape[0].compute()
        maxid = fai.index.compute().values.max()  # The new scaffolds must have an index after this one.
        # TODO this is wrong, we need to go back to the original scaffold here.
        broken.loc[:, "original_breakpoint"] = broken["orig_start"] + broken["breakpoint"] - 1
        broken = broken.sort_values(["orig_scaffold_index", "breakpoint"])
        if any(_ not in broken.columns for _ in breaks.columns):
            raise KeyError([(_, _ in broken.columns) for _ in breaks.columns])
        broken_next_cycle = broken.loc[broken["orig_scaffold_index"].duplicated()][breaks.columns]
        broken = broken.loc[~broken["orig_scaffold_index"].duplicated()]
        broken = _calculate_coordinates(broken, slop, maxid)
        broken = _create_children_dataframes(broken)
        # Now, for all those scaffolds that have multiple break points, we have to:
        # - remove the rightmost one from broken
        # - reset the name
        # - merge with the broken_next_cycle
        if broken_next_cycle.shape[0] > 0:
            key = "orig_scaffold_index"
            bait = broken[key].isin(broken_next_cycle.index.values)
            # For next cycle
            fnc = broken.loc[bait].groupby(key).tail(1).reset_index(drop=True)
            # Reset the scaffold index and scaffold name
            fnc.loc[:, "scaffold_index"] = fnc["orig_scaffold_index"]
            fnc.loc[:, "scaffold"] = fai[["scaffold"]].merge(fnc[["scaffold_index"]], on="scaffold_index")[
                ["scaffold"]]
            to_keep = broken.loc[bait].groupby(key).head(2).reset_index(drop=True)
            other_bait = to_keep["scaffold_index"].isin(fnc["scaffold_index"])
            broken = pd.concat([broken.loc[~bait],
                                to_keep.loc[~other_bait]]).reset_index(drop=True)
            broken_next_cycle = fnc.merge(broken_next_cycle.drop("length", axis=1), on="scaffold_index")
            broken_next_cycle.loc[:, "breakpoint"] = broken_next_cycle["breakpoint"] - broken_next_cycle["orig_start"]

        to_merge = dd.from_pandas(broken.sort_values("scaffold_index"), npartitions=1)
        fai = dd.concat([fai.reset_index(drop=False), to_merge]).set_index("scaffold_index")
        new_shape = fai.shape[0].compute()
        if initial_shape < new_shape:
            print("Finished cycle", cycle, ", increased scaffolds from", initial_shape, "to", new_shape)
        cycle += 1
        # TODO: now we have to move forward the coordinates for the remaining scaffolds.
        broken = broken_next_cycle.copy()

    fai = fai.drop_duplicates()
    _check = fai[["orig_scaffold_index", "orig_start", "length"]].drop_duplicates().shape[0].compute()
    assert fai.shape[0].compute() == _check

    fai_name = os.path.join(save_dir, "fai")
    dd.to_parquet(fai, fai_name, compute=True, compression="gzip", engine="pyarrow")
    return {"fai": fai_name}
