import pandas as pd
import numpy as np
from .rolling_join import rolling_join
import re
import dask.dataframe as dd
import dask.utils
from dask.base import compute


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


def _rebalance_ddf(ddf: dd.DataFrame, npartitions=None):
    """Repartition dask dataframe to ensure that partitions are roughly equal size.

    Assumes `ddf.index` is already sorted.
    """

    if not isinstance(ddf, dd.DataFrame):
        return ddf

    orig_shape = ddf.shape[0].compute()
    if not ddf.known_divisions:  # e.g. for read_parquet(..., infer_divisions=False)
        mins = ddf.index.map_partitions(dask.utils.M.min, meta=ddf.index)
        maxes = ddf.index.map_partitions(dask.utils.M.max, meta=ddf.index)
        mins, maxes = compute(mins, maxes)
        is_sorted = not (sorted(mins) != list(mins) or sorted(maxes) != list(maxes)
                      or any(a > b for a, b in zip(mins, maxes)))
        # print("IS SORTED:", is_sorted)
        if ddf.index.name is not None:
            ddf = ddf.reset_index().set_index(ddf.index.name, sorted=is_sorted)
        else:
            ddf = ddf.reset_index().set_index("index", sorted=is_sorted)

    index_counts = ddf.map_partitions(lambda _df: _df.index.value_counts().sort_index(),
                                      meta=int).compute()
    index = np.repeat(index_counts.index, index_counts.values)
    if npartitions is None:
        npartitions=ddf.npartitions
    divisions, _ = dd.io.io.sorted_division_locations(index, npartitions=npartitions)
    try:
        repartitioned = ddf.repartition(divisions=divisions, force=True)
    except ValueError:
        repartitioned = ddf

    assert repartitioned.shape[0].compute() == orig_shape
    return repartitioned
