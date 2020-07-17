import pandas as pd
import numpy as np
import functools
from time import ctime
import dask.dataframe as dd
from dask.distributed import Client
import os
from .collapse_bins import collapse_bins as cbn
from dask import delayed
import time
import dask.array as da
import logging
dask_logger = logging.getLogger("dask")


def _group_analyser(group: pd.DataFrame, binsize, cores=1):
    # assert group["scaffold_index"].unique().shape[0] == 1, group["scaffold_index"]
    # scaffold = group["scaffold_index"].head(1).values[0]
    bins = group.to_numpy().astype(np.int) + np.array([binsize, 0])
    counter = cbn(bins, binsize, cores)
    counter = np.vstack([np.fromiter(counter.keys(), dtype=np.int), np.fromiter(counter.values(),
                                                                                dtype=np.int)]).T
    assigned = np.hstack([np.repeat([group.index[0]], counter.shape[0]).reshape(counter.shape[0], 1),
                          counter])
    return assigned


def add_molecule_cov(assembly: dict, save_dir, client: Client, scaffolds=None, binsize=200, cores=1):

    dask_logger.warning("%s Begin add molecule coverage (10X)", ctime())
    info = assembly["info"]
    if not isinstance(info, dd.DataFrame):
        info = dd.read_parquet(info, infer_divisions=True)
    binsize = np.int(np.floor(binsize))

    if "molecules" not in assembly:
        raise KeyError("The assembly object does not have a molecule table; aborting")
    elif assembly["molecules"] is None:
        raise KeyError("The assembly object does not have a molecule table; aborting")

    if "mr_10x" in info.columns:
        raise KeyError("Assembly['info'] already has a mr_10x column; aborting")

    if isinstance(assembly["molecules"], str):
        molecules = dd.read_parquet(assembly["molecules"], infer_divisions=True)
    elif isinstance(assembly["molecules"], (dd.DataFrame, pd.DataFrame)):
        molecules = assembly["molecules"]
    else:
        raise TypeError(("molecules", type(assembly["molecules"])))

    if scaffolds is None:
        null = True
    else:
        dask_logger.warning("%s Extracting relevant scaffolds (10X)", ctime())
        info = info.loc[scaffolds].copy()
        dask_logger.warning("%s Extracted from info (10X)", ctime())
        try:
            dask_logger.warning("%s Calculating the index of molecules", time.ctime())
            idx = molecules.index.compute()
            dask_logger.warning("%s Calculating present", time.ctime())
            present = np.unique(idx.intersection(scaffolds))
            dask_logger.warning("%s Calculated present (size: %s)", time.ctime(), present.shape[0])
        except ValueError:
            print(molecules.index.head())
            print(type(scaffolds))
            print(scaffolds[:20])
            raise
        dask_logger.warning("%s Extracting from molecules (10X)", ctime())
        molecules = molecules.loc[present].copy()
        dask_logger.warning("%s Extracted from molecules (10X)", ctime())
        null = False

    if molecules.index.name == "scaffold_index":
        index = molecules.index.compute().values
    elif "scaffold_index" in molecules.columns:
        index = molecules["scaffold_index"].compute().to_numpy()
    else:
        raise KeyError("I cannot find the scaffold_index column in molecules!")

    # dask_logger.warning("%s Creating the temp dataframe (10X)", ctime())
    dask_logger.warning("%s Calculating the temp DF", time.ctime())
    temp_dataframe = pd.DataFrame().assign(
        scaffold_index=index,
        bin1=molecules["start"].compute().to_numpy() // binsize * binsize,
        bin2=molecules["end"].compute().to_numpy() // binsize * binsize,
    ).set_index("scaffold_index")
    dask_logger.warning("%s Calculated the temp DF", time.ctime())
    # temp_dataframe.index.name = "scaffold_index"
    # dask_logger.warning("%s Created the temp dataframe, querying (10X)", ctime())
    temp_dataframe = temp_dataframe.query("bin2 - bin1 > 2 * @binsize")[:]
    temp_dataframe = dd.from_pandas(temp_dataframe, npartitions=molecules.npartitions)
    # dask_logger.warning("%s Queried the temp dataframe (10X)", ctime())
    # Now let's persist the dataframe, and submit it to the Dask cluster
    _gr = functools.partial(_group_analyser, binsize=binsize, cores=2)
    dask_logger.warning("%s Calculating coverage per-group (10X)", ctime())
    finalised = temp_dataframe.groupby("scaffold_index").apply(_gr, meta=int).compute().values
    dask_logger.warning("%s Calculated coverage per-group (10X)", ctime())
    finalised = np.vstack(finalised)
    coverage_df = pd.DataFrame().assign(
        scaffold_index=finalised[:, 0],
        bin=finalised[:, 1],
        n=finalised[:, 2]).set_index("scaffold_index")
    dask_logger.warning("%s Created the coverage per-group (10X)", ctime())
    shape = coverage_df.shape[0]
    coverage_df = dd.from_pandas(coverage_df, chunksize=100000)

    if shape > 0:
        # info[,.(scaffold, length)][ff, on = "scaffold"]->ff
        dask_logger.warning("%s Merging on coverage DF (10X)", ctime())
        assert info.index.name == coverage_df.index.name
        coverage_df = client.scatter(coverage_df)
        info_length = client.scatter(info[["length"]])
        # coverage_df = dd.merge(info_length, coverage_df,
        #          left_index=True, right_index=True,
        #          how="right", chunksize=5000)
        func = delayed(dd.merge)(info_length, coverage_df,
                                 left_index=True, right_index=True,
                                 how="right")
        coverage_df = client.compute(func).result()
        dask_logger.warning("%s Merged on coverage DF (10X)", ctime())
        assert isinstance(coverage_df, dd.DataFrame), type(coverage_df)
        arr = coverage_df[["bin", "length"]].to_dask_array(lengths=True)
        distance = dd.from_array(np.minimum(
            arr[:, 0],
            (arr[:, 1] - arr[:, 0]) // binsize * binsize)).to_frame("d")
        distance.index = coverage_df.index
        coverage_df = coverage_df.assign(d=distance["d"])
        nbins = coverage_df.groupby(coverage_df.index.name)["d"].transform("size", meta=int)
        assert nbins.shape[0].compute() == coverage_df.shape[0].compute()
        coverage_df["nbin"] = nbins
        # assert isinstance(nbins, dd.DataFrame), type(nbins)
        # assert isinstance(coverage_df, dd.DataFrame), type(coverage_df)
        # coverage_df = dd.merge(coverage_df, nbins, how="left", left_index=True, right_index=True)
        dask_logger.warning("%s Calculated the distance metric (10X)", ctime())
        assert isinstance(coverage_df, dd.DataFrame), type(coverage_df)
        # Get the average coverage ACROSS ALL SCAFFOLDS by distance to the end of the bin.
        mn = coverage_df.groupby("d")["n"].mean().to_frame("mn")
        coverage_df = dd.merge(coverage_df.reset_index(drop=False),
                               mn, on="d", how="left").set_index("scaffold_index")
        coverage_df = coverage_df.eval("r = log(n / mn) / log(2)")
        dask_logger.warning("%s Calculated the mean coverage by distance (10X)", ctime())
        coverage_df["mr_10x"] = coverage_df["r"].groupby(
            coverage_df.index.name).transform("min", meta=coverage_df.r.dtype).to_dask_array()
        assert info.index.name == "scaffold_index"
        coverage_df.index = coverage_df.index.astype(info.index.dtype)
        try:
            info_mr = dd.merge(coverage_df[["mr_10x"]].drop_duplicates(), info, how="right", on="scaffold_index")
        except ValueError:
            dask_logger.critical(coverage_df[["mr_10x"]].drop_duplicates().shape[0].compute())
            dask_logger.critical("")
            dask_logger.critical(coverage_df[["mr_10x"]].dtypes)
            dask_logger.critical("")
            dask_logger.critical(info.dtypes)
            dask_logger.critical("")
            dask_logger.critical((coverage_df[["mr_10x"]].index.dtype, info.index.dtype))
            dask_logger.critical("")
            dask_logger.critical(coverage_df[["mr_10x"]].head(npartitions=-1))
            dask_logger.critical("")
            dask_logger.critical(info.head())
            dask_logger.critical("")
            import sys
            sys.exit(1)

        dask_logger.warning("%s Calculated the coverage ratio to the average (10X)", ctime())
        assert isinstance(info_mr, dd.DataFrame), type(info_mr)
        # assert isinstance(info_mr, dd.DataFrame), type(info_mr)
        info_mr = info_mr.drop("index", axis=1, errors="ignore").drop("scaffold", axis=1, errors="ignore")
        dask_logger.warning("%s Finished calculating the 10X coverage DF", ctime())
    else:
        info_mr = info.drop("mr_10x", axis=1, errors="ignore")
        # info_mr.drop("mr_10x", inplace=True, errors="ignore", axis=1)

    if null is True:
        assembly["info"] = info_mr
        assembly["molecule_cov"] = coverage_df
        assert isinstance(assembly["molecule_cov"], dd.DataFrame), type(assembly["molecule_cov"])
        assert isinstance(assembly["info"], dd.DataFrame), type(assembly["info"])
        assembly["mol_binsize"] = binsize
        for key in ["info", "molecule_cov"]:
            # dask_logger.warning("%s Starting to rebalance %s", ctime(), key)
            assembly[key] = assembly[key].repartition(partition_size="10MB")
            # dask_logger.warning("%s Finished rebalancing %s", ctime(), key)
            if save_dir is not None:
                fname = os.path.join(save_dir, key + "_10x")
                dd.to_parquet(assembly[key], fname, compression="gzip", engine="pyarrow", compute=True)
                assembly[key] = fname
                dask_logger.warning("%s Saved %s", ctime(), key)

        dask_logger.warning("%s Finished calculating 10X coverage", ctime())
        return assembly
    else:
        # TODO: this might need to be amended if we are going to checkpoint.
        dask_logger.warning("%s Dropping index then returning", ctime())
        info_mr = info_mr.drop("index", errors="ignore", axis=1)
        assembly = {"info": info_mr, "molecule_cov": coverage_df}
        dask_logger.warning("%s Finished calculating 10X coverage", ctime())
        return assembly
