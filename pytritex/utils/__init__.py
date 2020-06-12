import pandas as pd
import numpy as np
from .rolling_join import rolling_join
import re
import dask.dataframe as dd


def n50(array: np.array, p=0.5):
    array.sort()
    return array[np.where(array.cumsum() >= (1 - p) * array.sum())[0][0]]


def first(series: pd.Series):
    return series.iloc[0]


def second(series: pd.Series):
    if series.shape[0] < 2:
        return np.nan
    else:
        return series.iloc[1]


def unique_count(series: pd.Series):
    return np.unique(series).shape[0]


units = {"B": 1, "KB": 10**3, "MB": 10**6, "GB": 10**9, "TB": 10**12}
regex = re.compile(r'(\d+(?:\.\d+)?)\s*([kmgtp]?b)', re.IGNORECASE)


def parse_size(size):
    number, unit = regex.search(size.strip()).groups()
    unit = unit.upper()
    return int(float(number)*units[unit]), unit


def return_size(number, unit="GB"):
    return "{size}{unit}".format(size=max(round(number/units[unit], 2), 0.01),
                                 unit=unit)


def _rebalance_ddf(ddf: dd.DataFrame):
    """Repartition dask dataframe to ensure that partitions are roughly equal size.

    Assumes `ddf.index` is already sorted.
    """
    if not ddf.known_divisions:  # e.g. for read_parquet(..., infer_divisions=False)
        ddf = ddf.reset_index().set_index(ddf.index.name, sorted=True)
    index_counts = ddf.map_partitions(lambda _df: _df.index.value_counts().sort_index()).compute()
    index = np.repeat(index_counts.index, index_counts.values)
    divisions, _ = dd.io.io.sorted_division_locations(index, npartitions=ddf.npartitions)
    return ddf.repartition(divisions=divisions)