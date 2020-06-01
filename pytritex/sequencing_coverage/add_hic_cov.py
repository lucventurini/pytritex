import pandas as pd
import numpy as np
# import pandarallel
import functools
import dask.dataframe as dd
pd.options.mode.chained_assignment = 'raise'
from time import ctime
import multiprocessing as mp
import numexpr as ne
import re
from collections import deque
from .collapse_bins import collapse_bins as cbn


def _group_analyser(group, binsize):
    """Count how many pairs cover a given bin; return into a column called "n"."""
    ##     """Count how many pairs cover a given bin; return into a column called "n"."""
    # (scaffold_index, bin_group), group = group
    # # assert group["scaffold_index"].unique().shape[0] == 1, group["scaffold_index"]
    # # scaffold = group["scaffold_index"].head(1).values[0]
    # bin1 = group["bin1"].values
    # bin_series = pd.Series(itertools.starmap(range, pd.DataFrame().assign(
    #     bin1=group["bin1"] + binsize, bin2=group["bin2"], binsize=binsize).astype(np.int).values),
    #                 index=group.index)
    # assigned = group.assign(bin=bin_series).explode("bin").reset_index(drop=True).groupby(
    #     "bin").size().to_frame("n").reset_index().assign(scaffold_index=scaffold_index)
    # assigned.loc[:, "bin"] = pd.to_numeric(assigned["bin"].fillna(0), downcast="unsigned")
    # assigned.loc[:, "n"] = pd.to_numeric(assigned["n"].fillna(0), downcast="unsigned")
    # return assigned

    bins = group[["bin1", "bin2"]].to_numpy() + [binsize, 0]
    counter = cbn(bins, binsize)
    counter = np.vstack([np.fromiter(counter.keys(), dtype=np.int), np.fromiter(counter.values(),
                                                                                dtype=np.int)]).T
    assigned = np.hstack([np.repeat(group.index.to_numpy()[0],
                                    counter.shape[0]).reshape(counter.shape[0], 1),
                          counter])
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

    if binsize2 <= binsize * 3:
        raise ValueError("This function presumes that binsize is at least three times smaller than binsize2.\
Supplied values: {}, {}".format(binsize, binsize2))

    info = assembly["info"]
    assert isinstance(info, dd.DataFrame), type(info)
    if "mr" in info.columns or "mri" in info.columns:
        raise KeyError("Assembly[info] already has mr and/or mri columns; aborting.")
    fpairs = assembly["fpairs"]
    fpairs = fpairs.reset_index(drop=True)

    binsize = np.int(np.floor(binsize))
    if scaffolds is None:
        null = True
    else:
        info = info[info.scaffold_index.isin(scaffolds)]
        fpairs = fpairs[fpairs.scaffold_index1.isin(scaffolds)]
        null = False

    # Switch columns for those positions
    query = "scaffold_index1 == scaffold_index2 & pos1 > pos2"
    values = fpairs[["pos1", "pos2"]].compute().to_numpy()
    mask = np.repeat(fpairs.eval(query).compute().to_numpy(), 2).reshape(values.shape)
    values = np.where(mask, values[:, [1, 0]], values)
    fpairs["pos1"] = dd.from_array(values[:, 0]).reset_index(drop=True)
    fpairs["pos2"] = dd.from_array(values[:, 1]).reset_index(drop=True)
    fpairs = fpairs.drop_duplicates()
    # Bin positions of the match by BinSize; only select those bins where the distance between the two bins is
    # greater than the double of the binsize.
    query = "scaffold_index1 == scaffold_index2"
    temp_values = fpairs.query(query)[["scaffold_index1", "pos1", "pos2"]].compute().to_numpy()
    temp_frame = pd.DataFrame(
        {"scaffold_index": temp_values[:, 0],
         "bin1": temp_values[:, 1] // binsize * binsize,
         "bin2": temp_values[:, 2] // binsize * binsize
         }).astype(dtype={"scaffold_index": np.int32,
                   "bin1": np.int32, "bin2": np.int32}).query("bin2 - bin1 > 2 * @binsize").set_index("scaffold_index")
    # Create a greater bin group for each bin1, based on bin1 (default: 100x bin2).
    temp_frame.loc[:, "bin_group"] = temp_frame["bin1"] // binsize2
    # pandarallel.pandarallel.initialize(nb_workers=cores, use_memory_fs=use_memory_fs)
    pool = mp.Pool(processes=cores)
    _gr = functools.partial(_group_analyser, binsize=binsize)
    # Count how many pairs cover each smaller bin within the greater bin.
    results = deque()
    finalised = []
    for group in iter(temp_frame.groupby(["scaffold_index", "bin_group"])):
        while len(results) >= cores:
            finalised.append(results.popleft().get())
        results.append(pool.apply_async(_gr, (group[1],)))
    finalised.extend([res.get() for res in results])
    finalised = np.vstack(finalised)
    coverage_df = pd.DataFrame().assign(
        scaffold_index=finalised[:, 0],
        bin=finalised[:, 1],
        n=finalised[:, 2]).set_index("scaffold_index")
    pool.close()
    pool.join()

    if coverage_df.shape[0] > 0:
        print(ctime(), "Merging on coverage DF (HiC)")
        # Group by bin, count how covered is each bin
        coverage_df = coverage_df.groupby(["scaffold_index", "bin"]).agg(n=("n", "sum"))
        coverage_df = coverage_df.reset_index(level=1)
        coverage_df = dd.merge(info[["length"]],
                               coverage_df, on="scaffold_index", how="right").compute()
        # D is again the bin, but cutting it at the rightmost side
        lencol = coverage_df["length"].to_numpy().view()
        lencol = lencol.astype(copy=False, casting="unsafe",
                               dtype=re.sub(r"^u", "", lencol.dtype.name))
        bincol = coverage_df["bin"].to_numpy().view()
        bincol = bincol.astype(copy=False, casting="unsafe",
                               dtype=re.sub(r"^u", "", bincol.dtype.name))
        bincol = np.floor(ne.evaluate("(lencol - bincol) / binsize")).astype(lencol.dtype)
        coverage_df.loc[:, "d"] = pd.to_numeric(np.minimum(
            coverage_df["bin"], ne.evaluate("bincol * binsize")), downcast="integer")

        # Number of bins found in each scaffold. Take "n" as the column to count.
        nbins = coverage_df.groupby(level=0).size()
        coverage_df.loc[:, "nbin"] = np.repeat(nbins, nbins)
        # Mean number of pairs covering each bin (cut at the rightmost side)
        coverage_df.loc[:, "mn"] = coverage_df[["d", "n"]].groupby("d")["n"].transform("mean")
        # Logarithm of number of pairs in bin divided by mean of number of pairs in bin?
        coverage_df = coverage_df.eval("r = log(n/mn) / log(2)")
        # For each scaffold where the number of found bins is greater than minNbin,
        # calculate the minimum ratio (in log2) of the coverage per-bin divided by the mean coverage (in the scaffold)
        min_ratio = coverage_df.query("nbin >= @minNbin")["r"].groupby(level=0).agg(
            mr=("r", "min")).sort_values("mr")
        # For those scaffolds where we have at least one bin which is further away from the start than "innerDist",
        # calculate the minimum coverage *for those bins only*.
        min_internal_ratio = coverage_df.query(
            "nbin > @minNbin & d >= @innerDist")["r"].groupby(level=0).agg(mri=("r", "min"))
        coverage_df = min_internal_ratio.merge(
            min_ratio.merge(coverage_df, on="scaffold_index", how="right"),
            on="scaffold_index", how="right").drop("scaffold", axis=1, errors="ignore")
        info_mr = dd.merge(min_ratio, info, on="scaffold_index", how="right")
        info_mr = dd.merge(min_internal_ratio, info_mr, on="scaffold_index", how="right")
        info_mr = info_mr.drop("scaffold", axis=1, errors="ignore")
        print(ctime(), "Merged on coverage DF (HiC),", coverage_df.columns)
    else:
        info_mr = info.assign(mri=np.nan, mr=np.nan)

    info_mr = info_mr.persist()
    coverage_df = dd.from_pandas(coverage_df, npartitions=np.unique(coverage_df.index.to_numpy()).shape[0])

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
