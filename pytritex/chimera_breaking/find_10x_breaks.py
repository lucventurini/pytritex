import pandas as pd
import dask.dataframe as dd
import numpy as np


def find_10x_breaks(cov: dd.DataFrame, scaffolds=None,
                    interval=5e4, minNbin=20, dist=5e3, ratio=-3) -> pd.DataFrame:
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
        # Cov is indexed by scaffold index.
        cov = cov.loc[scaffolds]
    assert isinstance(cov, dd.DataFrame), type(cov)
    cov["b"] = cov["bin"] // interval
    cov = cov.persist()
    broken = cov.query("nbin >= @minNbin & r <= @ratio", local_dict=locals())
    bait = (np.minimum(broken["length"] - broken["bin"], broken["bin"]) >= dist).compute()
    broken = broken.compute().loc[bait, :]
    if broken.shape[0] == 0:
        return pd.DataFrame()
    broken = broken.sort_values("r").groupby(
        ["scaffold_index", "b"]).head(1).rename(columns={"bin": "breakpoint"},
                                                errors="raise")
    return broken
