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


def add_molecule_cov(assembly: dict, save_dir, client: Client, scaffolds=None, binsize=200, cores=1,
                     memory="20GB"):
    info = dd.read_parquet(assembly["info"])
    binsize = np.int(np.floor(binsize))

    if "molecules" not in assembly:
        raise KeyError("The assembly object does not have a molecule table; aborting")
    elif assembly["molecules"] is None:
        raise KeyError("The assembly object does not have a molecule table; aborting")

    if "mr_10x" in info.columns:
        raise KeyError("Assembly['info'] already has a mr_10x column; aborting")

    if scaffolds is None:
        molecules = dd.read_parquet(assembly["molecules"])
        null = True
    else:
        info = client.scatter(info)
        info = client.submit(dd.merge, info, scaffolds, on="scaffold_index", how="left",
                             resources={"process": 1, "MEMORY": memory})
        info = client.gather(info).drop("scaffold", axis=1, errors="ignore")
        molecules = dd.read_parquet(assembly["molecules"])
        molecules = client.submit(dd.merge, molecules, scaffolds, on="scaffold_index",
                                  how="left", resources={"process": 1, "MEMORY": memory})
        molecules = client.gather(molecules).drop("scaffold",axis=1, errors="ignore")
        null = False

    if molecules.index.name == "scaffold_index":
        index = molecules.index.compute().values
    elif "scaffold_index" in molecules.columns:
        index = molecules["scaffold_index"].compute().to_numpy()
    else:
        raise KeyError("I cannot find the scaffold_index column in molecules!")

    temp_dataframe = pd.DataFrame().assign(
        scaffold_index=index,
        bin1=molecules["start"].compute().to_numpy() // binsize * binsize,
        bin2=molecules["end"].compute().to_numpy() // binsize * binsize,
    ).set_index("scaffold_index")
    # temp_dataframe.index.name = "scaffold_index"
    temp_dataframe = temp_dataframe.query("bin2 - bin1 > 2 * @binsize")[:]

    # Now let's persist the dataframe, and submit it to the Dask cluster
    _gr = functools.partial(_group_analyser, binsize=binsize, cores=cores)
    # pool = mp.Pool(processes=cores)
    # pool.close()
    # results = deque()
    finalised = []
    grouped = temp_dataframe.groupby(level=0)
    for group in iter(grouped):
        # finalised.append(_gr(group[1]))
        finalised.append(_gr(group[1]))
    finalised = np.vstack(finalised)
    coverage_df = pd.DataFrame().assign(
        scaffold_index=finalised[:, 0],
        bin=finalised[:, 1],
        n=finalised[:, 2]).set_index("scaffold_index")
    shape = coverage_df.shape[0]
    coverage_df = dd.from_pandas(coverage_df, npartitions=min(shape // 100, 1000))
    print(coverage_df.head())

    if shape > 0:
        # info[,.(scaffold, length)][ff, on = "scaffold"]->ff
        print(ctime(), "Merging on coverage DF (10X)")
        assert info.index.name == coverage_df.index.name
        coverage_df = client.scatter(coverage_df)
        info_length = client.scatter(info[["length"]])
        # coverage_df = dd.merge(info_length, coverage_df,
        #          left_index=True, right_index=True,
        #          how="right", npartitions=1000)
        func = delayed(dd.merge)(info_length, coverage_df,
                                 left_index=True, right_index=True,
                                 how="right", npartitions=1000)
        coverage_df = client.compute(func).result()
        assert isinstance(coverage_df, dd.DataFrame), type(coverage_df)
        arr = coverage_df[["bin", "length"]].to_dask_array(lengths=True)
        distance = dd.from_array(np.minimum(
            arr[:, 0],
            (arr[:, 1] - arr[:, 0]) // binsize * binsize)).to_frame("d")
        distance.index = coverage_df.index
        coverage_df = coverage_df.assign(d=distance["d"])
        nbins = coverage_df.groupby(coverage_df.index.name).size().to_frame("nbin")
        assert isinstance(nbins, dd.DataFrame), type(nbins)
        assert isinstance(coverage_df, dd.DataFrame), type(coverage_df)
        coverage_df = dd.merge(coverage_df, nbins, how="left", left_index=True, right_index=True,
                               npartitions=1000)
        assert isinstance(coverage_df, dd.DataFrame), type(coverage_df)
        # Get the average coverage ACROSS ALL SCAFFOLDS by distance to the end of the bin.
        mn = coverage_df.groupby("d")["n"].mean().to_frame("mn")
        coverage_df = dd.merge(coverage_df.reset_index(drop=False), mn, on="d", how="left"
                               ).set_index("scaffold_index")
        assert isinstance(coverage_df, dd.DataFrame), type(coverage_df)
        coverage_df = coverage_df.eval("r = log(n / mn) / log(2)")
        assert isinstance(coverage_df, dd.DataFrame), type(coverage_df)
        mr_10x = coverage_df["r"].groupby(coverage_df.index.name
                                          ).agg("min").to_frame("mr_10x")
        assert isinstance(mr_10x, dd.DataFrame)
        print(mr_10x.head())
        # assert isinstance(info, (Delayed, dd.DataFrame))
        info = client.scatter(info)
        func = delayed(dd.merge)(mr_10x, info, how="right", on="scaffold_index", npartitions=1000)
        info_mr = client.compute(func).result()
        assert isinstance(info_mr, dd.DataFrame), type(info_mr)
        # assert isinstance(info_mr, dd.DataFrame), type(info_mr)
        info_mr = info_mr.drop("index", axis=1, errors="ignore").drop("scaffold", axis=1, errors="ignore")
        print(ctime(), "Merged on coverage DF (10X)")
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
            print(time.ctime(), "Storing", key)
            fname = os.path.join(save_dir, "joblib", "pytritex", "sequencing_coverage", key + "_10x")
            dd.to_parquet(assembly[key], fname, compression="gzip", engine="pyarrow", compute=True)
            assembly[key] = fname
        print(time.ctime(), "Finished storing 10X coverage.")
        return assembly
    else:
        # TODO: this might need to be amended if we are going to checkpoint.
        return {"info": info_mr, "molecule_cov": coverage_df}
