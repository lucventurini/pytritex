import pandas as pd
import numpy as np


def find_10x_breaks(cov: pd.DataFrame, scaffolds=None, interval=5e4, minNbin=20, dist=5e3, ratio=-3):
    """
    This function will take as input a coverage dataframe derived from 10X data.
    It will find those areas in scaffolds that are further away from the end point of the scaffold than the
    minimum distance AND have a coverage which is less than 8 times the average for the scaffold.
    :param cov: pd.DataFrame of the coverage (derived from 10X)
    :param scaffolds: optional; scaffold to consider for the analysis.
    :param interval:
    :param minNbin:
    :param dist:
    :param ratio:
    :return:
    """

    # cov = assembly["molecule_cov"].copy()
    if scaffolds is not None:
        cov = cov.merge(scaffolds, on="scaffold_index", how="right").drop("scaffold", axis=1, errors="ignore")
    try:
        cov.loc[:, "b"] = cov["bin"] // interval * interval
    except KeyError:
        raise KeyError(cov.head())
    bait = (cov["nbin"] >= minNbin) & (np.minimum(cov["bin"], cov["length"] - cov["bin"]) >= dist)
    bait &= (cov["r"] <= ratio)
    broken = cov.loc[bait, :]
    if broken.shape[0] == 0:
        return pd.DataFrame()
    broken = broken.sort_values("r").groupby(["scaffold_index", "b"]).head(1).rename(columns={"bin": "breakpoint"})
    return broken.copy()
