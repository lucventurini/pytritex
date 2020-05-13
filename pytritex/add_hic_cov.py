import pandas as pd
import numpy as np
import multiprocessing as mp
pd.options.mode.chained_assignment = 'raise'


def _group_analyser(group, binsize):
    assert group["scaffold"].unique().shape[0] == 1, group["scaffold"]
    scaffold = group["scaffold"].head(1).values[0]
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

    #  info<-assembly$info
    #
    #  if("mr" %in% colnames(info) | "mri" %in% colnames(info)){
    #   stop("assembly$info already has mr and/or mri columns; aborting.")
    #  }
    info = assembly["info"]
    if "mr" in info.columns or "mri" in info.columns:
        raise KeyError("Assembly[info] already has mr and/or mri columns; aborting.")
    fpairs = assembly["fpairs"]

    if scaffolds is None:
        null = True
    else:
        info = info[info.scaffold.isin(scaffolds)]
        fpairs = fpairs[fpairs.scaffold1.isin(scaffolds)]
        null = False

    bait = (fpairs["scaffold1"] == fpairs["scaffold2"]) & (fpairs["pos1"] < fpairs["pos2"])
    assert fpairs.loc[bait, :].shape[0] > 0 and fpairs.loc[bait, "pos1"].min() >= 0, fpairs.loc[bait, "pos1"]
    assert fpairs.loc[bait, :].shape[0] > 0 and fpairs.loc[bait, "pos2"].min() >= 0, fpairs.loc[bait, "pos2"]
    temp_frame = pd.DataFrame(
        {"scaffold": fpairs.loc[bait, "scaffold1"],
         "bin1": (fpairs.loc[bait, "pos1"] // binsize) * binsize,
         "bin2": (fpairs.loc[bait, "pos2"] // binsize) * binsize,
         }
    ).loc[lambda df: df["bin2"] - df["bin1"] > 2 * binsize, :]
    temp_frame.loc[:, "i"] = range(1, temp_frame.shape[0] + 1)
    temp_frame.loc[:, "b"] = temp_frame["scaffold"] + ":" + (temp_frame["bin1"] // binsize2).astype(str)
    fgrouped = temp_frame.groupby("b")
    pool = mp.Pool(cores)
    coverage_df = pd.concat(
        pool.starmap(_group_analyser,
        [(fgrouped.get_group(group), binsize) for group in fgrouped.groups.keys()]))
    pool.close()

    if coverage_df.shape[0] > 0:
        coverage_df = coverage_df.groupby(["scaffold", "bin"]).agg(n=("n", "sum")).reset_index(drop=False)
        coverage_df = info.loc[:, ["scaffold", "length"]].merge(coverage_df, on="scaffold", how="right")
        coverage_df.loc[:, "d"] = np.minimum(
            coverage_df["bin"],
            (coverage_df["length"] - coverage_df["bin"]) // binsize * binsize)
        coverage_df = coverage_df.merge(coverage_df.groupby("scaffold").agg(nbin=("length", "size")),
                      left_on="scaffold", right_index=True, how="left")

        coverage_df = coverage_df.merge(coverage_df.loc[:, ["d", "n"]].groupby("d").agg(mn=("n", "mean")),
                      left_on="d", right_index=True, how="left")
        coverage_df.loc[:, "r"] = np.log2(coverage_df["n"] / coverage_df["mn"])
        z = coverage_df.loc[coverage_df["nbin"] >= minNbin, ["scaffold", "r"]].groupby(
            "scaffold").agg(mr=("r", "min")).sort_values("mr")
        zi = coverage_df.loc[(coverage_df["nbin"] > minNbin) & (coverage_df["d"] >= innerDist),
                    ["scaffold", "r"]].groupby("scaffold").agg(mri=("r", "min"))
        coverage_df = z.merge(coverage_df, left_index=True, right_on="scaffold", how="right")
        coverage_df = zi.merge(coverage_df, left_index=True, right_on="scaffold", how="right")
        info_mr = z.merge(info, left_index=True, right_on="scaffold", how="right")
        info_mr = zi.merge(info_mr, left_index=True, right_on="scaffold", how="right")
        info_mr.reset_index(drop=True)
    else:
        info_mr = info.copy()
        info_mr = info_mr.assign(mri=np.nan, mr=np.nan)

    if null is True:
        assembly["info"] = info_mr
        assembly["cov"] = coverage_df
        assembly["binsize"] = binsize
        assembly["minNbin"] = minNbin
        assembly["innerDist"] = innerDist
        return assembly
    else:
        if "index" in info_mr.columns:
            del info_mr["index"]
        return {"info": info_mr, "cov": coverage_df}
