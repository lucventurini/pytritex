import pandas as pd
import numpy as np
# import pandarallel
import functools
import dask.dataframe as dd
from dask.distributed import Client
pd.options.mode.chained_assignment = 'raise'
from time import ctime
import os
from .collapse_bins import collapse_bins as cbn
from dask import delayed


def _group_analyser(group, binsize, cores=1):
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
    counter = cbn(bins, binsize, cores)
    counter = np.vstack([np.fromiter(counter.keys(), dtype=np.int), np.fromiter(counter.values(),
                                                                                dtype=np.int)]).T
    assigned = np.hstack([np.repeat(group.index.to_numpy()[0],
                                    counter.shape[0]).reshape(counter.shape[0], 1),
                          counter])
    return assigned


# Switch columns for those positions
def column_switcher(fpairs: dd.DataFrame):
    fpairs = fpairs.reset_index(drop=True)
    query = "scaffold_index1 == scaffold_index2 & pos1 > pos2"
    values = fpairs[["pos1", "pos2"]].to_dask_array(lengths=True)
    mask = np.repeat(fpairs.eval(query).compute().to_numpy(), 2).reshape(values.shape)
    values = np.where(mask, values[:, [1, 0]], values)
    fpairs["pos1"] = values[:, 0]
    fpairs["pos2"] = values[:, 1]
    fpairs = fpairs.drop_duplicates()
    return fpairs


def add_hic_cov(assembly, save_dir, client: Client,
                scaffolds=None, binsize=1e3, binsize2=1e5, minNbin=50, innerDist=1e5, cores=1,
                memory="20GB"):
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

    info = dd.read_parquet(assembly["info"])
    assert isinstance(info, dd.DataFrame), type(info)
    if "mr" in info.columns or "mri" in info.columns:
        raise KeyError("Assembly[info] already has mr and/or mri columns; aborting.")
    fpairs = dd.read_parquet(assembly["fpairs"])

    binsize = np.int(np.floor(binsize))
    if scaffolds is None:
        null = True
    else:
        info = info[info.scaffold_index.isin(scaffolds)]
        fpairs = fpairs[fpairs.scaffold_index1.isin(scaffolds)]
        null = False

    fpairs = client.submit(column_switcher, fpairs).result()
    assert isinstance(fpairs, dd.DataFrame)
    # Bin positions of the match by BinSize; only select those bins where the distance between the two bins is
    # greater than the double of the binsize.
    query = "scaffold_index1 == scaffold_index2"
    temp_values = fpairs.query(query)[["scaffold_index1", "pos1", "pos2"]].compute().to_numpy()
    temp_frame = pd.DataFrame(
        {"scaffold_index": temp_values[:, 0],
         "bin1": temp_values[:, 1] // binsize * binsize,
         "bin2": temp_values[:, 2] // binsize * binsize
         }).astype(dtype={"scaffold_index": np.int32,
                   "bin1": np.int32, "bin2": np.int32}).query(
        "bin2 - bin1 > 2 * @binsize").set_index("scaffold_index")
    # Create a greater bin group for each bin1, based on bin1 (default: 100x bin2).
    temp_frame.loc[:, "bin_group"] = temp_frame["bin1"] // binsize2
    # pandarallel.pandarallel.initialize(nb_workers=cores, use_memory_fs=use_memory_fs)
    _gr = functools.partial(_group_analyser, binsize=binsize, cores=cores)
    # Count how many pairs cover each smaller bin within the greater bin.
    finalised = []
    for group in iter(temp_frame.groupby(["scaffold_index", "bin_group"])):
        finalised.append(_gr(group[1]))
    finalised = np.vstack(finalised)
    coverage_df = pd.DataFrame().assign(
        scaffold_index=finalised[:, 0],
        bin=finalised[:, 1],
        n=finalised[:, 2]).set_index("scaffold_index")

    if coverage_df.shape[0] > 0:
        print(ctime(), "Merging on coverage DF (HiC)")
        # Group by bin, count how covered is each bin, reset the index so it is only by scaffold_index
        coverage_df = coverage_df.groupby(["scaffold_index", "bin"]).agg(n=("n", "sum")).reset_index(level=1)
        coverage_df = client.scatter(dd.from_pandas(coverage_df, chunksize=100000))
        info_length = client.scatter(info[["length"]])
        func = delayed(dd.merge)(info_length, coverage_df, on="scaffold_index", how="right")
        coverage_df = client.compute(func).result()
        assert isinstance(coverage_df, dd.DataFrame)
        # D is again the bin, but cutting it at the rightmost side
        arr = coverage_df[["bin", "length"]].to_dask_array(lengths=True)
        distance = dd.from_array(np.minimum(
            arr[:, 0],
            (arr[:, 1] - arr[:, 0]) // binsize * binsize)).to_frame("d")
        distance.index = coverage_df.index
        coverage_df = coverage_df.assign(d=distance["d"])
        # Number of bins found in each scaffold. Take "n" as the column to count.
        nbins = coverage_df.groupby(coverage_df.index.name).size().to_frame("nbin")
        coverage_df = dd.merge(coverage_df, nbins, how="left", left_index=True, right_index=True)
        # Mean number of pairs covering each bin (cut at the rightmost side)
        mn = coverage_df.groupby("d")["n"].mean().to_frame("mn")
        coverage_df = dd.merge(coverage_df.reset_index(drop=False), mn, on="d", how="left"
                               ).set_index("scaffold_index")
        # Logarithm of number of pairs in bin divided by mean of number of pairs in bin?
        coverage_df = coverage_df.eval("r = log(n/mn) / log(2)")
        assert isinstance(coverage_df, dd.DataFrame), type(coverage_df)
        # For each scaffold where the number of found bins is greater than minNbin,
        # calculate the minimum ratio (in log2) of the coverage per-bin divided by the mean coverage (in the scaffold)
        min_ratio = coverage_df[coverage_df["nbin"] >= minNbin]["r"]
        min_ratio = min_ratio.groupby("scaffold_index").min().to_frame("mr").compute().sort_values("mr")
        min_ratio = dd.from_pandas(min_ratio, chunksize=100000)
        # For those scaffolds where we have at least one bin which is further away from the start than "innerDist",
        # calculate the minimum coverage *for those bins only*.
        min_internal_ratio = coverage_df[(coverage_df["nbin"] > minNbin) & (coverage_df["d"] >= innerDist)]["r"]
        min_internal_ratio = min_internal_ratio.groupby("scaffold_index").min().to_frame("mri")
        coverage_df = dd.merge(min_ratio, coverage_df, on="scaffold_index", how="right")
        coverage_df = dd.merge(min_internal_ratio, coverage_df, on="scaffold_index", how="right")
        assert isinstance(coverage_df, dd.DataFrame), type(coverage_df)
        info_mr = dd.merge(min_ratio, info, on="scaffold_index", how="right")
        info_mr = dd.merge(min_internal_ratio, info_mr, on="scaffold_index", how="right")
        info_mr = info_mr.drop("scaffold", axis=1, errors="ignore")
        assert isinstance(info_mr, dd.DataFrame)
        print(ctime(), "Merged on coverage DF (HiC),", coverage_df.columns)
    else:
        info_mr = info.assign(mri=np.nan, mr=np.nan)

    if null is True:
        assembly["info"] = info_mr
        assembly["cov"] = coverage_df
        assembly["binsize"] = binsize
        assembly["minNbin"] = minNbin
        assembly["innerDist"] = innerDist
        for key in ["info", "cov"]:
            print(ctime(), "Storing", key, "for HiC")
            fname = os.path.join(save_dir, key + "_hic")
            dd.to_parquet(assembly[key], fname, compression="gzip", engine="pyarrow", compute=True)
            assembly[key] = fname
        print(ctime(), "Finished storing HiC")
        return assembly
    else:
        # TODO: this might need to be amended
        info_mr.drop("index", inplace=True, errors="ignore", axis=1)
        return {"info": info_mr, "cov": coverage_df}
