import pandas as pd
import numpy as np
# import pandarallel
import functools
import dask.dataframe as dd
pd.options.mode.chained_assignment = 'raise'
from time import ctime
import os
from .collapse_bins import collapse_bins as cbn
# from dask import delayed
import logging
dask_logger = logging.getLogger("dask")


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


def add_hic_cov(assembly, save_dir, scaffolds=None, binsize=1e3, binsize2=1e5, minNbin=50, innerDist=1e5):
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

    dask_logger.debug("%s Starting to add the HiC coverage", ctime())

    if binsize2 <= binsize * 3:
        raise ValueError("This function presumes that binsize is at least three times smaller than binsize2.\
Supplied values: {}, {}".format(binsize, binsize2))

    dask_logger.debug("%s Loading the info and fpairs data", ctime())
    if isinstance(assembly["info"], str):
        info = dd.read_parquet(assembly["info"])
    else:
        info = assembly["info"]
    assert isinstance(info, dd.DataFrame), type(info)
    if "mr" in info.columns or "mri" in info.columns:
        raise KeyError("Assembly[info] already has mr and/or mri columns; aborting.")
    if isinstance(assembly["fpairs"], str):
        fpairs = dd.read_parquet(assembly["fpairs"])
    else:
        assert isinstance(assembly["fpairs"], dd.DataFrame)
        fpairs = assembly["fpairs"]
    dask_logger.debug("%s Loaded the info and fpairs data", ctime())

    binsize = np.int(np.floor(binsize))
    if scaffolds is None:
        null = True
        # dask_logger.debug("%s Removing duplicates from fpairs", ctime())
        # fpairs = client.submit(column_switcher, fpairs).result()
        # dask_logger.debug("%s Removed duplicates from fpairs", ctime())
    else:
        dask_logger.debug("%s Selecting data from info and fpairs", ctime())
        info_index = info.index.compute()
        assert info.index.name == "scaffold_index"
        info = info.loc[info_index.intersection(scaffolds).values].copy()
        dask_logger.debug("%s Selected data from info", ctime())
        fpairs = fpairs[(fpairs.scaffold_index1.isin(scaffolds)) |
                        (fpairs.scaffold_index2.isin(scaffolds))].copy()
        dask_logger.debug("%s Selected data from fpairs", ctime())
        null = False

    original_info_size = info.shape[0].compute()
    assert isinstance(fpairs, dd.DataFrame)
    assert fpairs.shape[0].compute() > 0
    # Bin positions of the match by BinSize; only select those bins where the distance between the two bins is
    # greater than the double of the binsize.
    query = "scaffold_index1 == scaffold_index2"
    try:
        temp_values = fpairs.query(query)[["scaffold_index1", "pos1", "pos2"]].compute().to_numpy()
    except ValueError:
        print(fpairs.head(npartitions=-1, n=5))
        raise ValueError(fpairs.head(npartitions=-1, n=5))
    dask_logger.debug("%s Creating the temporary dataframe", ctime())
    temp_frame = pd.DataFrame(
        {"scaffold_index": temp_values[:, 0],
         "bin1": temp_values[:, 1] // binsize * binsize,
         "bin2": temp_values[:, 2] // binsize * binsize
         }).astype(dtype={"scaffold_index": np.int32,
                   "bin1": np.int32, "bin2": np.int32}).query(
        "bin2 - bin1 > 2 * @binsize").set_index("scaffold_index")
    # Create a greater bin group for each bin1, based on bin1 (default: 100x bin2).
    temp_frame.loc[:, "bin_group"] = temp_frame["bin1"] // binsize2
    dask_logger.debug("%s Created the temporary dataframe", ctime())
    temp_frame = dd.from_pandas(temp_frame, npartitions=fpairs.npartitions)
    # pandarallel.pandarallel.initialize(nb_workers=cores, use_memory_fs=use_memory_fs)
    _gr = functools.partial(_group_analyser, binsize=binsize)
    # Count how many pairs cover each smaller bin within the greater bin.
    finalised = temp_frame.groupby("scaffold_index").apply(_gr, meta=int).compute().values
    dask_logger.debug("%s Calculated the coverage bins", ctime())
    finalised = np.vstack(finalised)
    coverage_df = pd.DataFrame().assign(
        scaffold_index=finalised[:, 0],
        bin=finalised[:, 1],
        n=finalised[:, 2]).set_index("scaffold_index")

    if coverage_df.shape[0] > 0:
        dask_logger.debug("%s Merging on coverage DF (HiC)", ctime())
        # Group by bin, count how covered is each bin, reset the index so it is only by scaffold_index
        coverage_df = coverage_df.groupby(["scaffold_index", "bin"]).agg(n=("n", "sum")).reset_index(level=1)
        coverage_df = dd.merge(info[["length"]], coverage_df, on="scaffold_index", how="right")
        assert isinstance(coverage_df, dd.DataFrame)
        if coverage_df.query("length != length").shape[0].compute() > 0:
            dask_logger.critical(
                "Something went very wrong with merging the dataframes. The length column must always be populated.")
            raise AssertionError
        # D is again the bin, but cutting it at the rightmost side
        dask_logger.debug("%s Calculating the distance metric (HiC)", ctime())
        arr = coverage_df[["bin", "length"]].to_dask_array(lengths=True)
        distance = dd.from_array(np.minimum(
            arr[:, 0],
            (arr[:, 1] - arr[:, 0]) // binsize * binsize)).to_frame("d")
        distance.index = coverage_df.index
        coverage_df = coverage_df.assign(d=distance["d"])
        assert coverage_df.query("d != d").shape[0].compute() == 0
        assert isinstance(coverage_df, dd.DataFrame)
        # Number of bins found in each scaffold. Take "n" as the column to count.
        _col = coverage_df.columns[0]
        dask_logger.debug("%s Calculating the nbin metric (HiC)", ctime())
        coverage_df["nbin"] = coverage_df.groupby(
            coverage_df.index.name)[_col].transform("size", meta=int).to_dask_array(lengths=True)
        assert isinstance(coverage_df, dd.DataFrame)
        _ = coverage_df.head(npartitions=-1, n=5)
        # Mean number of pairs covering each bin (cut at the rightmost side)
        dask_logger.debug("%s Calculating the mean mn metric (HiC)", ctime())
        mn = coverage_df.reset_index(drop=True)[["d", "n"]].groupby("d")["n"].mean().to_frame("mn")
        coverage_df = coverage_df.reset_index(drop=False)
        coverage_df = coverage_df.drop("mn", axis=1, errors="ignore").set_index("d").merge(mn, how="left", on="d")
        coverage_df = coverage_df.reset_index(drop=False).set_index("scaffold_index")
        # Logarithm of number of pairs in bin divided by mean of number of pairs in bin?
        assert coverage_df.query("d != d").shape[0].compute() == 0
        coverage_df = coverage_df.eval("r = log(n/mn) / log(2)")
        assert isinstance(coverage_df, dd.DataFrame), type(coverage_df)
        assert "r" in coverage_df.columns
        # For each scaffold where the number of found bins is greater than minNbin,
        # calculate the minimum ratio (in log2) of the coverage per-bin divided by the mean coverage (in the scaffold)
        dask_logger.debug("%s Calculating the min_ratio metric (HiC)", ctime())
        bait = coverage_df["nbin"] >= minNbin
        min_ratio = coverage_df.loc[bait, ["r"]].groupby("scaffold_index").min().rename(
            columns={"r": "mr"}).drop_duplicates(ignore_index=False)
        assert isinstance(min_ratio, dd.DataFrame)
        assert coverage_df.index.name == "scaffold_index"
        assert min_ratio.index.name == coverage_df.index.name, min_ratio.index.name
        assert "r" in coverage_df.columns
        coverage_df = coverage_df.merge(min_ratio, how="left", on=coverage_df.index.name)
        # Test
        _ = coverage_df.head(npartitions=-1, n=5)
        assert "r" in coverage_df.columns
        # For those scaffolds where we have at least one bin which is further away from the start than "innerDist",
        # calculate the minimum coverage *for those bins only*.
        dask_logger.debug("%s Calculating the min internal ratio mri (HiC)", ctime())
        bait = (coverage_df["nbin"] > minNbin) & (coverage_df["d"] >= innerDist)
        min_internal_ratio = coverage_df.loc[bait, ["r"]].groupby(
            coverage_df.index.name).min().rename(columns={"r": "mri"}).drop_duplicates(ignore_index=False)
        assert coverage_df.index.name == "scaffold_index"
        assert min_internal_ratio.index.name == coverage_df.index.name, min_internal_ratio.index.name
        coverage_df = coverage_df.merge(min_internal_ratio, how="left",
                                        left_index=True, right_index=True)
        assert isinstance(coverage_df, dd.DataFrame), type(coverage_df)
        assert info.index.name == min_ratio.index.name == "scaffold_index"
        info_mr = dd.merge(min_ratio, info, left_index=True, right_index=True, how="right")
        assert info_mr.index.name == min_internal_ratio.index.name == "scaffold_index"
        info_mr = dd.merge(min_internal_ratio, info_mr, left_index=True,
                           right_index=True, how="right")
        if "scaffold" in info_mr.columns:
            info_mr = info_mr.drop("scaffold", axis=1)
        assert isinstance(info_mr, dd.DataFrame)
        assert info_mr.shape[0].compute() == original_info_size
        dask_logger.debug("%s Merged on coverage DF (HiC)", ctime())
    else:
        info_mr = info.assign(mri=np.nan, mr=np.nan)

    if null is True:
        assembly["info"] = info_mr
        assembly["cov"] = coverage_df
        assembly["binsize"] = binsize
        assembly["minNbin"] = minNbin
        assembly["innerDist"] = innerDist
        for key in ["info", "cov"]:
            if save_dir is not None:
                fname = os.path.join(save_dir, key + "_hic")
                assembly[key] = assembly[key].repartition(partition_size="100MB", force=True)
                dd.to_parquet(assembly[key], fname, compression="gzip", engine="pyarrow", compute=True)
                dask_logger.debug("%s Saved %s for HiC", ctime(), key)
                assembly[key] = fname
        dask_logger.debug("%s Finished storing HiC", ctime())
        return assembly
    else:
        dask_logger.debug("%s Preparing the HiC DF for further use", ctime())
        info_mr = info_mr.drop("index", errors="ignore", axis=1)
        assembly = {"info": info_mr, "cov": coverage_df}
        dask_logger.debug("%s Returning the HiC DF for further use", ctime())
        return assembly
