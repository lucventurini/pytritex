import pandas as pd
import numpy as np
import itertools
import pandarallel


def add_molecule_cov(assembly: dict, scaffolds=None, binsize=200, cores=1):
    info = assembly["info"]
    print("Starting adding molecule coverage")
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
        info = info.merge(scaffolds, on="scaffold", how="left")
        mol = assembly["molecules"].merge(scaffolds, on="scaffold", how="left")
        null = False

    f = mol
    assert not f["start"].isna().any()
    assert not f["end"].isna().any()
    f.loc[:, "bin1"] = f["start"] // binsize * binsize
    f.loc[:, "bin2"] = f["end"] // binsize * binsize
    f = f.loc[f.eval("bin2 - bin1 > 2 *{binsize}".format(binsize=binsize)), :].copy()
    f.loc[:, "bin"] = pd.Series(itertools.starmap(range, pd.DataFrame().assign(
        bin1=f.loc[:, "bin1"] + binsize, bin2=f["bin2"], binsize=binsize).values))
    pandarallel.pandarell.initialize(nb_workers=cores)
    ff = f.groupby("scaffold").parallel_apply(
        lambda group: group.explode("bin").groupby(["scaffold", "bin"]).size().to_frame("n").reset_index(
            drop=False)).reset_index(level=0, drop=True)

    if ff.shape[0] > 0:
        # info[,.(scaffold, length)][ff, on = "scaffold"]->ff
        ff = info.loc[:, ["scaffold", "length"]].merge(ff, on="scaffold", how="right")
        ff.loc[:, "d"] = np.min([
            ff["bin"], ((ff["length"] - ff["bin"]) // binsize) * binsize
        ], axis=0)
        ff = ff.set_index("scaffold").merge(
            ff.groupby("scaffold").size().to_frame("nbin"),
            left_index=True, right_index=True, how="left").reset_index(drop=False)
        ff = ff.set_index("d").merge(
            ff.groupby("d").agg(mn=("n", "mean")), left_index=True, right_index=True).reset_index(drop=False)
        ff.loc[:, "r"] = np.log2(ff["n"] / ff["mn"])
        info_mr = info.merge(ff.groupby("scaffold").agg(mr_10x=("r", "min")), left_on="scaffold",
                             right_index=True, how="left").reset_index(drop=True)
        if "index" in info_mr.columns:
            del info_mr["index"]
    else:
        info_mr = info.copy()
        del info_mr["mr_10x"]

    print("Molecule cov (add_mol):", ff.shape[0], ff.columns)
    if null is True:
        assembly["info"] = info_mr
        assembly["molecule_cov"] = ff
        assembly["mol_binsize"] = binsize
        return assembly
    else:
        return {"info": info_mr, "molecule_cov": ff}
