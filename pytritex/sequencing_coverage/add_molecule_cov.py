import pandas as pd
import numpy as np
import itertools
import pandarallel
import functools
from time import ctime


def _group_analyser(group: pd.DataFrame, binsize):
    assert group["scaffold_index"].unique().shape[0] == 1, group["scaffold_index"]
    scaffold = group["scaffold_index"].head(1).values[0]
    bin_series = pd.Series(itertools.starmap(range, pd.DataFrame().assign(
        bin1=group["bin1"] + binsize, bin2=group["bin2"], binsize=binsize).astype(np.int).values),
                           index=group.index,
                           name="bin")
    assigned = group.assign(bin=bin_series).explode("bin").groupby(
        ["scaffold_index", "bin"]).size().to_frame("n").reset_index(
        drop=False).assign(scaffold_index=scaffold)
    return assigned


def add_molecule_cov(assembly: dict, scaffolds=None, binsize=200, cores=1, use_memory_fs=True):
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
    temp_dataframe.loc[:, "bin1"] = temp_dataframe["start"] // binsize * binsize
    temp_dataframe.loc[:, "bin2"] = temp_dataframe["end"] // binsize * binsize
    temp_dataframe = temp_dataframe.loc[temp_dataframe.eval(
        "bin2 - bin1 > 2 *{binsize}".format(binsize=binsize)), :].copy().astype({"bin1": np.int,
                                                                                 "bin2": np.int})
    pandarallel.pandarallel.initialize(nb_workers=cores, use_memory_fs=use_memory_fs)
    _gr = functools.partial(_group_analyser, binsize=binsize)
    coverage_df = temp_dataframe.groupby("scaffold_index").parallel_apply(_gr).reset_index(level=0, drop=True)

    if coverage_df.shape[0] > 0:
        # info[,.(scaffold, length)][ff, on = "scaffold"]->ff
        print(ctime(), "Merging on coverage DF (10X)")
        try:
            coverage_df = info.loc[:, ["scaffold_index", "length"]].merge(
                coverage_df, on="scaffold_index", how="right")
        except ValueError:
            print(coverage_df.head())
            raise ValueError((coverage_df["scaffold_index"].dtype, info["scaffold_index"].dtype))
        coverage_df.loc[:, "d"] = np.minimum(
            coverage_df["bin"], ((coverage_df["length"] - coverage_df["bin"]) // binsize) * binsize)
        coverage_df.loc[:, "nbin"] = coverage_df.groupby("scaffold_index")["scaffold_index"].transform("size")
        coverage_df.loc[:, "mn"] = coverage_df.groupby("d")["n"].transform("mean")
        coverage_df.loc[:, "r"] = np.log2(coverage_df["n"] / coverage_df["mn"])
        __left = coverage_df.groupby("scaffold_index").agg(mr_10x=("r", "min"))
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
