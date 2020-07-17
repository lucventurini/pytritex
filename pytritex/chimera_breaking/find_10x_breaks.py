import pandas as pd
import dask.dataframe as dd
import numpy as np
import time
import logging
dask_logger = logging.getLogger("dask")
from typing import Union


def find_10x_breaks(cov: dd.DataFrame, interval=5e4, minNbin=20, dist=5e3, ratio=-3) -> Union[pd.DataFrame, None]:
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
    if isinstance(cov, str):
        cov = dd.read_parquet(cov, infer_divisions=True)
    assert isinstance(cov, dd.DataFrame), type(cov)
    dask_logger.debug("%s Calculating b", time.ctime())
    cov["b"] = cov["bin"] // interval
    dask_logger.debug("%s Calculating b", time.ctime())
    broken = cov.loc[(cov["nbin"] >= minNbin) & (cov["r"] <= ratio), :].compute()
    if broken.index.name == "scaffold_index":
        broken = broken.reset_index(drop=False)
    dask_logger.debug("%s Extracted broken, creating the bait", time.ctime())
    bait = (np.minimum(broken["length"] - broken["bin"], broken["bin"]) >= dist)
    bindex = np.unique(broken.index[bait])
    if bindex.shape[0] == 0:
        return None

    broken = broken.loc[bindex, :].reset_index(drop=False).astype({"scaffold_index": int})
    dask_logger.debug("%s Extracting using the bait", time.ctime())
    # We exploit here the fact that each specific index is going to be in a different partition.

    broken = broken.sort_values("r", ascending=True).groupby(["scaffold_index", "b"]).head(1)
    broken = broken.rename(columns={"bin": "breakpoint"}, errors="raise")
    if broken.index.name is not None:
        broken = broken.reset_index(drop=False)
    dask_logger.debug("%s Extracted breakpoints, reindexing", time.ctime())
    broken = broken.drop_duplicates(subset=["scaffold_index", "breakpoint"])
    broken = broken.set_index("scaffold_index")
    dask_logger.debug("%s Found breakpoints", time.ctime())
    return broken
