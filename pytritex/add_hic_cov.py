import pandas as pd
import numpy as np
import multiprocessing as mp


def _group_analyser(group, binsize):
    assert group["scaffold"].unique().shape[0] == 1, group["scaffold"]
    scaffold = group["scaffold"].head(1)

    return group.apply(lambda row: list(range(int(row["bin1"] + binsize), int(row["bin2"]), int(binsize))),
                       axis=1).explode().reset_index(drop=True).to_frame("bin").groupby(
        "bin").size().to_frame("n").reset_index().assign(scaffold=scaffold)


# Calculate physical coverage with Hi-C links in sliding windows along the scaffolds
def add_hic_cov(assembly, scaffolds=None, binsize=1e3, binsize2=1e5, minNbin=50, innerDist=1e5, cores=1):
    """
    :param assembly:
    :param scaffolds:
    :param binsize:
    :param binsize2:
    :param minNbin:
    :param innerDist:
    :param cores:
    :return:
    """

    info = assembly["info"]
    if "mr" in info.columns or "mri" in info.columns:
        raise KeyError("Assembly[info] already has mr and/or mri columns; aborting.")
    fpairs = assembly["fpairs"]
    if scaffolds is None:
        # scaffolds = info["scaffold"]
        null = True
    else:
        info = info[info.scaffold.isin(scaffolds)]
        fpairs = fpairs[fpairs.scaffold1.isin(scaffolds)]
        null = False

    bait = (fpairs.scaffold1 == fpairs.scaffold2) & (fpairs.pos1 < fpairs.pos2)
    f = pd.DataFrame(
        {"scaffold": fpairs.loc[bait]["scaffold1"],
         "bin1": (fpairs.loc[bait]["pos1"] // binsize) * binsize,
         "bin2": (fpairs.loc[bait]["pos2"] // binsize) * binsize,
         }
    ).loc[lambda df: df.bin2 - df.bin1 > 2 * binsize, :]
    f.loc[:, "bin1"] = f.loc[:, "bin1"].astype(int)
    f.loc[:, "bin2"] = f.loc[:, "bin2"].astype(int)
    f["i"] = range(1, f.shape[0] + 1)  # Not sure about this yet
    f["b"] = f["scaffold"].astype(str) + ":" + (f["bin1"] // binsize2).astype(str)
    fgrouped = f.groupby("b")
    pool = mp.Pool(cores)
    ff = pd.concat(pool.starmap(_group_analyser,
                               [(fgrouped.get_group(group), binsize)
                                for group in fgrouped.groups.keys()]))
    print(ff.head())

    if ff.shape[0] > 0:
        # Sum again as we had clustered by "b"
        ff = ff.groupby(["scaffold", "bin"]).agg(n=("n", "sum")).reset_index(drop=False)
        ff = info.loc[:, ["scaffold", "length"]].merge(ff, on="scaffold", how="right")
        assert ff[~ff["bin"].isna()].shape[0] > 0
        assert ff[~ff["n"].isna()].shape[0] > 0
        assert ff[~ff["length"].isna()].shape[0] > 0, (
            info.loc[:, ["scaffold", "length"]].head(), ff.loc[:, ["scaffold", "bin", "n"]].head())
        ff.loc[:, "d"] = np.min(
            [ff["bin"], (ff["length"] - ff["bin"]) // binsize * binsize], axis=0)
        assert ff[~ff["d"].isna()].shape[0] > 0
        ff = ff.merge(ff.groupby("scaffold").agg(nbin=("length", "size")),
                      left_on="scaffold", right_index=True, how="left")
        ff = ff.merge(ff.loc[:, ["d", "n"]].groupby("d").agg(mn=("n", "mean")),
                      left_on="d", right_index=True, how="left")
        ff.loc[:, "r"] = np.log2(ff["n"] / ff["mn"])
        print(ff.head())
        z = ff.loc[ff["nbin"] >= minNbin, ["scaffold", "r"]].groupby(
            "scaffold").agg(mr=("r", "min")).sort_values("mr")
        zi = ff.loc[(ff["nbin"] > minNbin) & (ff["d"] >= innerDist),
                    ["scaffold", "r"]].groupby("scaffold").agg(mri=("r", "min"))
        ff = z.merge(ff, left_index=True, right_on="scaffold", how="right")
        ff = zi.merge(ff, left_index=True, right_on="scaffold", how="right")
        info_mr = z.merge(info, left_index=True, right_on="scaffold", how="right")
        info_mr = zi.merge(info_mr, left_index=True, right_on="scaffold", how="right")
    else:
        info_mr = info.copy()
        info_mr = info_mr.assign(mri=np.nan, mr=np.nan)

    if null is True:
        assembly["info"] = info_mr
        assembly["cov"] = ff
        assembly["binsize"] = binsize
        assembly["minNbin"] = minNbin
        assembly["innerDir"] = innerDist
        return assembly
    else:
        return {"info": info_mr, "cov": ff}
