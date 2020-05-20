import pandas as pd
import numpy as np
import pandarallel
import functools
import itertools
pd.options.mode.chained_assignment = 'raise'
from time import ctime


def _group_analyser(group, binsize):
    """Count how many pairs cover a given bin; return into a column called "n"."""
    assert group["scaffold_index"].unique().shape[0] == 1, group["scaffold_index"]
    scaffold = group["scaffold_index"].head(1).values[0]
    bin_series = pd.Series(itertools.starmap(range, pd.DataFrame().assign(
        bin1=group["bin1"] + binsize, bin2=group["bin2"], binsize=binsize).astype(np.int).values),
                    index=group.index)
    assigned = group.assign(bin=bin_series).explode("bin").reset_index(drop=True).groupby(
        "bin").size().to_frame("n").reset_index().assign(scaffold_index=scaffold)
    return assigned


def add_hic_cov(assembly, scaffolds=None, binsize=1e3, binsize2=1e5, minNbin=50, innerDist=1e5, cores=1,
                use_memory_fs=True):
    """
    Calculate physical coverage with Hi-C links in sliding windows along the scaffolds.
    :param assembly:
    :param scaffolds:
    :param binsize:
    :param binsize2:
    :param minNbin:
    :param innerDist:
    :param cores:
    :param use_memory_fs: for pandarallel.
    :return:
    """

    #  info<-assembly$info
    #
    #  if("mr" %in% colnames(info) | "mri" %in% colnames(info)){
    #   stop("assembly$info already has mr and/or mri columns; aborting.")
    #  }

    if binsize2 <= binsize * 3:
        raise ValueError("This function presumes that binsize is at least three times smaller than binsize2.\
Supplied values: {}, {}".format(binsize, binsize2))

    info = assembly["info"]
    if "mr" in info.columns or "mri" in info.columns:
        raise KeyError("Assembly[info] already has mr and/or mri columns; aborting.")
    fpairs = assembly["fpairs"].copy()

    if scaffolds is None:
        null = True
    else:
        info = info[info.scaffold_index.isin(scaffolds)]
        fpairs = fpairs[fpairs.scaffold_index1.isin(scaffolds)]
        null = False

    # First step: get pairs on the same scaffold
    bait = (fpairs["scaffold_index1"] == fpairs["scaffold_index2"])
    try:
        bait2 = bait & (fpairs["pos1"] > fpairs["pos2"])
    except KeyError:
        raise KeyError(fpairs.columns)
    # Switch columns for those positions
    fpairs.loc[bait2, ["pos1", "pos2"]] = fpairs.loc[bait2, ["pos2", "pos1"]].values
    # assert fpairs.loc[bait, :].shape[0] > 0 and fpairs.loc[bait, "pos1"].min() >= 0, fpairs.loc[bait, "pos1"]
    # assert fpairs.loc[bait, :].shape[0] > 0 and fpairs.loc[bait, "pos2"].min() >= 0, fpairs.loc[bait, "pos2"]

    # Bin positions of the match by BinSize; only select those bins where the distance between the two bins is
    # greater than the double of the binsize.
    temp_frame = pd.DataFrame(
        {"scaffold_index": fpairs.loc[bait, "scaffold_index1"],
         "bin1": (fpairs.loc[bait, "pos1"] // binsize) * binsize,
         "bin2": (fpairs.loc[bait, "pos2"] // binsize) * binsize,
         }
    ).loc[lambda df: df["bin2"] - df["bin1"] > 2 * binsize, :].astype({"bin1": np.int, "bin2": np.int})
    # Create a greater bin group for each bin1, based on bin1 (default: 100x bin2).
    temp_frame.loc[:, "bin_group"] = (temp_frame["bin1"] // binsize2)
    pandarallel.pandarallel.initialize(nb_workers=cores, use_memory_fs=use_memory_fs)
    _gr = functools.partial(_group_analyser, binsize=binsize)
    # Count how many pairs cover each smaller bin within the greater bin.
    coverage_df = temp_frame.groupby(["scaffold_index",
                                      "bin_group"]).parallel_apply(_gr).reset_index(level=0, drop=True)
    # coverage_df = temp_frame.groupby(["scaffold_index", "bin_group"]).apply(_gr).reset_index(level=0, drop=True)

    if coverage_df.shape[0] > 0:
        print(ctime(), "Merging on coverage DF (HiC)")
        # Group by bin, count how covered is each bin
        coverage_df = coverage_df.groupby(["scaffold_index", "bin"]).agg(n=("n", "sum")).reset_index(drop=False)
        try:
            coverage_df = info.loc[:,
                          ["scaffold_index", "length"]].merge(coverage_df, on="scaffold_index", how="right")
        except KeyError:
            raise KeyError(info.columns)
        # D is again the bin, but cutting it at the rightmost side
        coverage_df.loc[:, "d"] = np.minimum(
            coverage_df["bin"],
            (coverage_df["length"] - coverage_df["bin"]) // binsize * binsize)
        # Number of bins found in each scaffold
        coverage_df.loc[:, "nbin"] = coverage_df.groupby("scaffold_index")["scaffold_index"].transform("size")
        # Mean number of pairs covering each bin (cut at the rightmost side)
        coverage_df.loc[:, "mn"] = coverage_df.loc[:, ["d", "n"]].groupby("d")["n"].transform("mean")
        # Logarithm of number of pairs in bin divided by mean of number of pairs in bin?
        coverage_df.loc[:, "r"] = np.log2(coverage_df["n"] / coverage_df["mn"])
        # For each scaffold where the number of found bins is greater than minNbin,
        # calculate the minimum ratio (in log2) of the coverage per-bin divided by the mean coverage (in the scaffold)
        z = coverage_df.loc[coverage_df["nbin"] >= minNbin, ["scaffold_index", "r"]].groupby(
            "scaffold_index").agg(mr=("r", "min")).sort_values("mr")
        # For those scaffolds where we have at least one bin which is further away from the start than "innerDist",
        # calculate the minimum coverage *for those bins only*.
        zi = coverage_df.loc[(coverage_df["nbin"] > minNbin) & (coverage_df["d"] >= innerDist),
                    ["scaffold_index", "r"]].groupby("scaffold_index").agg(mri=("r", "min"))
        coverage_df = zi.merge(
            z.merge(coverage_df, left_index=True, right_on="scaffold_index", how="right"),
            left_index=True, right_on="scaffold_index", how="right").reset_index(drop=True).drop(
            "scaffold", axis=1, errors="ignore")
        info_mr = zi.merge(z.merge(info, left_index=True, right_on="scaffold_index", how="right"),
                           left_index=True, right_on="scaffold_index", how="right").reset_index(drop=True).drop(
            "scaffold", axis=1, errors="ignore")
        print(ctime(), "Merged on coverage DF (HiC),", coverage_df.columns)
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