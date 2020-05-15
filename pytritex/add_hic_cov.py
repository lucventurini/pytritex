import pandas as pd
import numpy as np
import pandarallel
import functools
import itertools
pd.options.mode.chained_assignment = 'raise'
from time import ctime


def _group_analyser(group, binsize):
    assert group["scaffold"].unique().shape[0] == 1, group["scaffold"]
    scaffold = group["scaffold"].head(1).values[0]
    bin_series = pd.Series(itertools.starmap(range, pd.DataFrame().assign(
        bin1=group["bin1"] + binsize, bin2=group["bin2"], binsize=binsize).astype(np.int).values),
                    index=group.index)
    assigned = group.assign(bin=bin_series).explode("bin").reset_index(drop=True).groupby(
        "bin").size().to_frame("n").reset_index().assign(scaffold=scaffold)
    return assigned


def add_hic_cov(assembly, scaffolds=None, binsize=1e3, binsize2=1e5, minNbin=50, innerDist=1e5, cores=1):
    """
    Calculate physical coverage with Hi-C links in sliding windows along the scaffolds.
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
    ).loc[lambda df: df["bin2"] - df["bin1"] > 2 * binsize, :].astype({"bin1": np.int, "bin2": np.int})
    temp_frame.loc[:, "i"] = np.arange(1, temp_frame.shape[0] + 1, dtype=np.int)
    temp_frame.loc[:, "b"] = temp_frame["scaffold"] + ":" + (temp_frame["bin1"] // binsize2).astype(str)
    pandarallel.pandarallel.initialize(nb_workers=cores)
    _gr = functools.partial(_group_analyser, binsize=binsize)
    coverage_df = temp_frame.groupby("b").parallel_apply(_gr).reset_index(level=0, drop=True)
    # coverage_df = temp_frame.groupby("b").apply(_gr).reset_index(level=0, drop=True)

    if coverage_df.shape[0] > 0:
        print(ctime(), "Merging on coverage DF (HiC)")
        coverage_df = coverage_df.groupby(["scaffold", "bin"]).agg(n=("n", "sum")).reset_index(drop=False)
        try:
            coverage_df = info.loc[:, ["scaffold", "length"]].merge(coverage_df, on="scaffold", how="right")
        except KeyError:
            raise KeyError(info.columns)
        coverage_df.loc[:, "d"] = np.minimum(
            coverage_df["bin"],
            (coverage_df["length"] - coverage_df["bin"]) // binsize * binsize)
        coverage_df.loc[:, "nbin"] = coverage_df.groupby("scaffold")["scaffold"].transform("size")
        coverage_df.loc[:, "mn"] = coverage_df.loc[:, ["d", "n"]].groupby("d")["n"].transform("mean")
        coverage_df.loc[:, "r"] = np.log2(coverage_df["n"] / coverage_df["mn"])
        z = coverage_df.loc[coverage_df["nbin"] >= minNbin, ["scaffold", "r"]].groupby(
            "scaffold").agg(mr=("r", "min")).sort_values("mr")
        zi = coverage_df.loc[(coverage_df["nbin"] > minNbin) & (coverage_df["d"] >= innerDist),
                    ["scaffold", "r"]].groupby("scaffold").agg(mri=("r", "min"))
        coverage_df = zi.merge(
            z.merge(coverage_df, left_index=True, right_on="scaffold", how="right"),
            left_index=True, right_on="scaffold", how="right").reset_index(drop=True)
        info_mr = zi.merge(z.merge(info, left_index=True, right_on="scaffold", how="right"),
                           left_index=True, right_on="scaffold", how="right").reset_index(drop=True)
        print(ctime(), "Merged on coverage DF (HiC)")
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
        info_mr.drop("index", inplace=True, errors="ignore", axis=1)
        return {"info": info_mr, "cov": coverage_df}
