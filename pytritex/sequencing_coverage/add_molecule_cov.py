import pandas as pd
import numpy as np
import itertools
# import pandarallel
import functools
from time import ctime
import multiprocessing as mp
import numexpr as ne
import re
from collections import deque


def _group_analyser(group: pd.DataFrame, binsize):
    # assert group["scaffold_index"].unique().shape[0] == 1, group["scaffold_index"]
    scaffold_index, group = group
    # scaffold = group["scaffold_index"].head(1).values[0]
    bin1 = group["bin1"]
    bin_series = pd.Series(itertools.starmap(range, pd.DataFrame().assign(
        bin1=group["bin1"] + binsize, bin2=group["bin2"], binsize=binsize).astype(np.int).values),
                           index=group.index,
                           name="bin")
    assigned = group.assign(bin=bin_series).explode("bin").groupby(
        ["scaffold_index", "bin"]).size().to_frame("n").reset_index(
        drop=False).assign(scaffold_index=scaffold_index)
    assigned.loc[:, "bin"] = pd.to_numeric(assigned["bin"].fillna(0), downcast="unsigned")
    assigned.loc[:, "n"] = pd.to_numeric(assigned["n"].fillna(0), downcast="unsigned")
    return assigned


def add_molecule_cov(assembly: dict, scaffolds=None, binsize=200, cores=1):
    info = assembly["info"]
    print("Starting adding molecule coverage")
    binsize = np.int(np.floor(binsize))
    if "molecules" not in assembly or assembly["molecules"] is None or assembly["molecules"].shape[0] == 0:
        raise KeyError("The assembly object does not have a molecule table; aborting")

    if "mr_10x" in info.columns:
        raise KeyError("Assembly['info'] already has a mr_10x column; aborting")

    assert not assembly["molecules"]["start"].isna().any()
    assert not assembly["molecules"]["end"].isna().any()
    if scaffolds is None:
        # scaffolds = info.loc[:, "scaffold"]
        mol = assembly["molecules"].copy()
        null = True
    else:
        info = info.merge(scaffolds, on="scaffold_index", how="left").drop("scaffold", axis=1, errors="ignore")
        mol = assembly["molecules"].merge(scaffolds, on="scaffold_index", how="left")
        null = False

    temp_dataframe = mol
    assert not temp_dataframe["start"].isna().any()
    assert not temp_dataframe["end"].isna().any()
    _start = temp_dataframe["start"].values
    _start = np.floor(ne.evaluate("_start / binsize")).astype(_start.dtype)
    temp_dataframe.loc[:, "bin1"] = ne.evaluate("_start * binsize").astype(temp_dataframe["end"].dtype)
    assert not temp_dataframe["bin1"].isna().any()
    _end = temp_dataframe["end"].values
    _end = np.floor(ne.evaluate("_end / binsize")).astype(_end.dtype)
    temp_dataframe.loc[:, "bin2"] = ne.evaluate("_end * binsize").astype(temp_dataframe["end"].dtype)
    assert not temp_dataframe["bin2"].isna().any()
    bin1, bin2 = temp_dataframe["bin1"].values, temp_dataframe["bin2"].values
    temp_dataframe = temp_dataframe.loc[
                     ne.evaluate("bin2 - bin1 > 2 * binsize"), :].copy()
    temp_dataframe.loc[:, "bin1"] = pd.to_numeric(temp_dataframe["bin1"], downcast="unsigned")
    temp_dataframe.loc[:, "bin2"] = pd.to_numeric(temp_dataframe["bin2"], downcast="unsigned")
    # pandarallel.pandarallel.initialize(nb_workers=cores, use_memory_fs=use_memory_fs)
    _gr = functools.partial(_group_analyser, binsize=binsize)
    pool = mp.Pool(processes=cores)
    results = deque()
    finalised = []
    for group in iter(temp_dataframe[["scaffold_index", "bin1", "bin2"]].groupby("scaffold_index")):
        while len(results) > 100:
            finalised.append(results.popleft().get())
        results.append(pool.apply_async(_gr, group))
    coverage_df = pd.concat(finalised).reset_index(level=0, drop=True)
    pool.close()
    pool.join()

    if coverage_df.shape[0] > 0:
        # info[,.(scaffold, length)][ff, on = "scaffold"]->ff
        print(ctime(), "Merging on coverage DF (10X)")
        try:
            coverage_df = info.loc[:, ["scaffold_index", "length"]].merge(
                coverage_df, on="scaffold_index", how="right")
        except ValueError:
            print(coverage_df.head())
            raise ValueError((coverage_df["scaffold_index"].dtype, info["scaffold_index"].dtype))
        lencol, bincol = coverage_df["length"].values.view(), coverage_df["bin"].values.view()
        lencol = lencol.astype(copy=False, casting="unsafe",
                               dtype=re.sub(r"^u", "", lencol.dtype.name))
        bincol = bincol.astype(copy=False, casting="unsafe",
                               dtype=re.sub(r"^u", "", bincol.dtype.name))
        bincol = np.floor(ne.evaluate("(lencol - bincol) / binsize")).astype(lencol.dtype)

        coverage_df.loc[:, "d"] = pd.to_numeric(np.minimum(
            coverage_df["bin"],  ne.evaluate("bincol * binsize")), downcast="integer")
        coverage_df.loc[:, "nbin"] = pd.to_numeric(coverage_df.groupby(
            "scaffold_index")["scaffold_index"].transform("size"), downcast="unsigned")
        coverage_df.loc[:, "mn"] = pd.to_numeric(coverage_df.groupby("d")["n"].transform("mean"),
                                                 downcast="float")
        coverage_df.loc[:, "r"] = pd.to_numeric(np.log2(coverage_df["n"] / coverage_df["mn"]),
                                                downcast="float")
        __left = coverage_df[["scaffold_index", "r"]].groupby("scaffold_index").agg(mr_10x=("r", "min"))
        __left.loc[:, "mr_10x"] = pd.to_numeric(__left["mr_10x"], downcast="float")
        info_mr = __left.merge(info, how="right", right_on="scaffold_index", left_index=True).drop(
            "index", axis=1, errors="ignore").drop("scaffold", axis=1, errors="ignore")
        print(ctime(), "Merged on coverage DF (10X)")
    else:
        info_mr = info.copy()
        # info_mr.drop("mr_10x", inplace=True, errors="ignore", axis=1)

    print("Molecule cov (add_mol):", coverage_df.shape[0], coverage_df.columns)
    if null is True:
        assembly["info"] = info_mr
        assembly["molecule_cov"] = coverage_df
        assembly["mol_binsize"] = binsize
        return assembly
    else:
        return {"info": info_mr, "molecule_cov": coverage_df}
