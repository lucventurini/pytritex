import pandas as pd
import numpy as np
import logging
import functools
from time import ctime
import numexpr as ne
import re
import dask.dataframe as dd
import multiprocessing as mp
from collections import deque
from .collapse_bins import collapse_bins as cbn


def _group_analyser(group: pd.DataFrame, binsize):
    # assert group["scaffold_index"].unique().shape[0] == 1, group["scaffold_index"]
    # scaffold = group["scaffold_index"].head(1).values[0]
    bins = group.to_numpy().astype(np.int) + np.array([binsize, 0])
    counter = cbn(bins, binsize)
    counter = np.vstack([np.fromiter(counter.keys(), dtype=np.int), np.fromiter(counter.values(),
                                                                                dtype=np.int)]).T
    assigned = np.hstack([np.repeat([group.index[0]], counter.shape[0]).reshape(counter.shape[0], 1),
                          counter])
    return assigned


def add_molecule_cov(assembly: dict, scaffolds=None, binsize=200, cores=1):
    info = assembly["info"]
    binsize = np.int(np.floor(binsize))
    if "molecules" not in assembly:
        raise KeyError("The assembly object does not have a molecule table; aborting")
    elif assembly["molecules"] is None:
        raise KeyError("The assembly object does not have a molecule table; aborting")

    if "mr_10x" in info.columns:
        raise KeyError("Assembly['info'] already has a mr_10x column; aborting")

    if scaffolds is None:
        molecules = assembly["molecules"]
        null = True
    else:
        info = info.merge(scaffolds, on="scaffold_index", how="left").drop("scaffold", axis=1, errors="ignore")
        molecules = dd.merge(
            assembly["molecules"], scaffolds, on="scaffold_index", how="left").drop(
            "scaffold",axis=1, errors="ignore")
        null = False

    temp_dataframe = pd.DataFrame().assign(
        bin1=molecules["start"].compute() // binsize * binsize,
        bin2=molecules["end"].compute() // binsize * binsize,
    ).set_index(molecules["scaffold_index"].compute())
    temp_dataframe.index.name = "scaffold_index"
    temp_dataframe = temp_dataframe.query("bin2 - bin1 > 2 * @binsize")[:]

    # Now let's persist the dataframe, and submit it to the Dask cluster

    _gr = functools.partial(_group_analyser, binsize=binsize)
    pool = mp.Pool(processes=cores)
    # pool.close()
    results = deque()
    finalised = []
    grouped = temp_dataframe.groupby(level=0)
    for group in iter(grouped):
        # finalised.append(_gr(group[1]))
        while len(results) >= cores:
            finalised.append(results.popleft().get())
        results.append(pool.apply_async(_gr, (group[1],)))
    finalised.extend([res.get() for res in results])
    pool.close()
    finalised = np.vstack(finalised)
    coverage_df = pd.DataFrame().assign(
        scaffold_index=finalised[:, 0],
        bin=finalised[:, 1],
        n=finalised[:, 2]).set_index("scaffold_index")

    if coverage_df.shape[0] > 0:
        # info[,.(scaffold, length)][ff, on = "scaffold"]->ff
        print(ctime(), "Merging on coverage DF (10X)")
        try:
            coverage_df = dd.merge(info.loc[:, ["length"]], coverage_df,
                                   on="scaffold_index", how="right").compute().reset_index(drop=False)
        except ValueError:
            print(coverage_df.head())
            raise ValueError((coverage_df["scaffold_index"].dtype, info["scaffold_index"].dtype))
        lencol, bincol = coverage_df["length"].to_numpy().view(), coverage_df["bin"].to_numpy().view()
        lencol = lencol.astype(copy=False, casting="unsafe",
                               dtype=re.sub(r"^u", "", lencol.dtype.name))
        bincol = bincol.astype(copy=False, casting="unsafe",
                               dtype=re.sub(r"^u", "", bincol.dtype.name))
        bincol = np.floor(ne.evaluate("(lencol - bincol) / binsize")).astype(lencol.dtype)

        coverage_df.loc[:, "d"] = pd.to_numeric(np.minimum(
            coverage_df["bin"],  ne.evaluate("bincol * binsize")), downcast="integer")
        coverage_df.loc[:, "nbin"] = pd.to_numeric(coverage_df.groupby(
            "scaffold_index")["scaffold_index"].transform("size"), downcast="signed")
        coverage_df.loc[:, "mn"] = pd.to_numeric(coverage_df.groupby("d")["n"].transform("mean"),
                                                 downcast="float")
        coverage_df.loc[:, "r"] = pd.to_numeric(np.log2(coverage_df["n"] / coverage_df["mn"]),
                                                downcast="float")
        __left = coverage_df[["scaffold_index", "r"]].groupby("scaffold_index").agg(mr_10x=("r", "min"))
        __left.loc[:, "mr_10x"] = pd.to_numeric(__left["mr_10x"], downcast="float")
        info_mr = dd.merge(__left, info, how="right", on="scaffold_index").drop(
            "index", axis=1, errors="ignore").drop("scaffold", axis=1, errors="ignore")

        print(ctime(), "Merged on coverage DF (10X)")
    else:
        info_mr = info.drop("mr_10x", axis=1, errors="ignore")
        # info_mr.drop("mr_10x", inplace=True, errors="ignore", axis=1)
    info_mr.persist()
    print("Molecule cov (add_mol):", coverage_df.shape[0], coverage_df.columns)
    if null is True:
        assembly["info"] = info_mr
        assembly["molecule_cov"] = coverage_df
        assembly["mol_binsize"] = binsize
        return assembly
    else:
        return {"info": info_mr, "molecule_cov": coverage_df}
