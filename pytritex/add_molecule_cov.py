import pandas as pd
import numpy as np
import multiprocessing as mp


def _group_analyser(group, scaffold, binsize):
    # df.apply(lambda row: list(range(int(row["bin1"] + binsize), int(row["bin2"]), binsize)),
    #          axis=1).explode().reset_index(drop=True).to_frame("bin").groupby("bin").size().to_frame("N").reset_index(
    #     drop=False).assign(scaffold="a")
    return group.apply(lambda row: list(range(int(row["bin1"] + binsize), int(row["bin2"]), binsize)),
                       axis=1).explode().reset_index(drop=True).to_frame("bin").groupby(
        "bin").size().to_frame("n").reset_index().assign(scaffold=scaffold)


def add_molecule_cov(assembly: dict, scaffolds=None, binsize=200, cores=1):
    info = assembly["info"]
    if "molecules" not in assembly or assembly["molecules"] is None or assembly["molecules"].shape[0] == 0:
        raise KeyError("The assembly object does not have a molecule table; aborting")

    if "mr_10x" in info.columns:
        raise KeyError("Assembly['info'] already has a mr_10x column; aborting")

    if scaffolds is None:
        # scaffolds = info.loc[:, "scaffold"]
        mol = assembly["molecules"].copy()
        null = True
    else:
        info = info.merge(scaffolds, on="scaffold", how="right")
        mol = assembly["molecules"].merge(scaffolds, on="scaffold", how="right")
        null = False

    f = mol
    f.loc[:, "bin1"] = f["start"] // binsize * binsize
    f.loc[:, "bin2"] = f["end"] // binsize * binsize
    fgrouped = f.groupby("scaffold")
    pool = mp.Pool(cores)
    ff = pd.concat(pool.starmap(_group_analyser,
                               [(fgrouped.get_group(group), group, binsize)
                                for group in fgrouped.groups.keys()]))

    if ff.shape[0] > 0:
        # info[,.(scaffold, length)][ff, on = "scaffold"]->ff
        ff = info.loc[:, ["scaffold", "length"]].merge(ff, on="scaffold", how="right")
        ff.loc[:, "d"] = np.min([
            ff["bin"], ((ff["length"] - ff["bin"]) // binsize) * binsize
        ], axis=0)
        ff = ff.set_index("scaffold").merge(
            ff.groupby("scaffold").size().to_frame("nbin"),
            left_index=True, right_index=True, how="left").reset_index(drop=False)
        ff = ff.set_index("d").merge(
            ff.groupby("d").agg(mn=("n", "mean")), left_index=True, right_index=True).reset_index(drop=False)
        ff.loc[:, "r"] = np.log2(ff["n"] / ff["mn"])
        info_mr = info.merge(ff.groupby("scaffold").agg(mr_10x=("r", "min")), left_on="scaffold", right_index=True,
                        how="left")
    else:
        info_mr = info.copy()
        info_mr.loc[:, "mr_10x"] = np.nan

    if null is True:
        assembly["info"] = info_mr
        assembly["molecule_cov"] = ff
        assembly["mol_binsize"] = binsize
        return assembly
    else:
        return {"info": info_mr, "molecule_cov": ff}
