import pandas as pd
import numpy as np
import functools
from time import ctime
import dask.dataframe as dd
import os
from .collapse_bins import collapse_bins as cbn


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


def add_molecule_cov(assembly: dict, save_dir, scaffolds=None, binsize=200, cores=1):
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
        info = dd.merge(info, scaffolds, on="scaffold_index", how="left").drop("scaffold", axis=1, errors="ignore")
        molecules = dd.read_parquet(assembly["molecules"])
        molecules = dd.merge(molecules, scaffolds, on="scaffold_index",
                             how="left").drop("scaffold",axis=1, errors="ignore")
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

    if coverage_df.shape[0] > 0:
        # info[,.(scaffold, length)][ff, on = "scaffold"]->ff
        print(ctime(), "Merging on coverage DF (10X)")
        try:
            coverage_df = dd.merge(info[["length"]], coverage_df,
                                   on="scaffold_index", how="right").compute()
        except ValueError:
            print(coverage_df.head())
            raise ValueError((coverage_df["scaffold_index"].dtype, info["scaffold_index"].dtype))
        coverage_df["d"] = pd.to_numeric(
            np.minimum(
                coverage_df["bin"],
                ((coverage_df["length"] - coverage_df["bin"]) // binsize) * binsize
            ), downcast="integer")

        nbins = coverage_df.groupby(level=0).size().to_frame("nbin")
        coverage_df = coverage_df.merge(nbins, how="left", on="scaffold_index")
        # Get the average coverage ACROSS ALL SCAFFOLDS by distance to the end of the bin.
        coverage_df["mn"] = pd.to_numeric(
            coverage_df.reset_index(drop=False).groupby("d")["n"].transform("mean"), downcast="float")
        coverage_df = coverage_df.eval("r = log(n / mn) / log(2)")
        __left = coverage_df["r"].groupby(level=0).agg(mr_10x=("r", "min"))
        __left.loc[:, "mr_10x"] = pd.to_numeric(__left["mr_10x"], downcast="float")
        info_mr = dd.merge(__left, info, how="right", on="scaffold_index").drop(
            "index", axis=1, errors="ignore").drop("scaffold", axis=1, errors="ignore")

        print(ctime(), "Merged on coverage DF (10X)")
    else:
        info_mr = info.drop("mr_10x", axis=1, errors="ignore")
        # info_mr.drop("mr_10x", inplace=True, errors="ignore", axis=1)
    info_mr = info_mr.persist()
    coverage_df = dd.from_pandas(coverage_df, npartitions=np.unique(coverage_df.index.values).shape[0])
    coverage_df = coverage_df.persist()
    print("Molecule cov (add_mol):", coverage_df.shape[0], coverage_df.columns)

    if null is True:
        assembly["info"] = info_mr
        assembly["molecule_cov"] = coverage_df
        assembly["mol_binsize"] = binsize
        for key in ["info", "molecule_cov"]:
            fname = os.path.join(save_dir, "joblib", "pytritex", "sequencing_coverage", key + "_10x")
            dd.to_parquet(assembly[key], fname, compression="gzip", engine="pyarrow", compute=True)
            assembly[key] = fname
        return assembly
    else:
        # TODO: this might need to be amended if we are going to checkpoint.
        return {"info": info_mr, "molecule_cov": coverage_df}
