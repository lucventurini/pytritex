import pandas as pd
import numpy as np


def find_10x_breaks(assembly: dict, scaffolds=None, interval=5e4, minNbin=20, dist=5e3, ratio=-3):
    cov = assembly["molecule_cov"].copy()
    if scaffolds is not None:
        cov = cov.merge(scaffolds, on="scaffold", how="right")
    try:
        cov.loc[:, "b"] = cov["bin"] // interval * interval
    except KeyError:
        raise KeyError(cov.head())
    bait = (cov["nbin"] >= minNbin) & (np.minimum(cov["bin"], cov["length"] - cov["bin"]) >= dist)
    bait &= (cov["r"] <= ratio)
    e = cov.loc[bait, :]
    if e.shape[0] == 0:
        return None
    e = e.sort_values("r").groupby(["scaffold", "b"]).head(1).rename(columns={"bin": "br"})
    return e.copy()[:]
