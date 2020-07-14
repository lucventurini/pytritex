import pandas as pd
import numpy as np
from .rolling_join import rolling_join
import re
import dask.dataframe as dd
import dask.utils
from dask.base import compute
import itertools as it
import typing


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


def _agg_second(grouped):
    def internal(c):
        if (c != c).all():
            return [np.nan]
        f = [_ for _ in c if _ == _]
        f = [_ if isinstance(_, list) else [_] for _ in f]
        return list(it.chain.from_iterable(f))
    g = grouped.apply(internal)
    g = g.apply(lambda s: np.nan if len(s) == 1 else s[1])
    return g


second_agg = dd.Aggregation("second", lambda s: s.apply(list), _agg_second)


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


def _rebalance_ddf(ddf: dd.DataFrame, npartitions=None, target_memory=None):
    """Repartition dask dataframe to ensure that partitions are roughly equal size.

    Assumes `ddf.index` is already sorted.
    """

    if not isinstance(ddf, dd.DataFrame):
        return ddf

    ddf = ddf.persist()
    _ = ddf.head(npartitions=-1, n=1)
    import logging
    dask_logger = logging.getLogger("dask")
    orig_shape = ddf.shape[0].compute()
    dask_logger.warning("Original shape: %s", orig_shape)
    if not ddf.known_divisions:  # e.g. for read_parquet(..., infer_divisions=False)
        mins = ddf.index.map_partitions(dask.utils.M.min, meta=ddf.index)
        maxes = ddf.index.map_partitions(dask.utils.M.max, meta=ddf.index)
        mins, maxes = compute(mins, maxes)
        is_sorted = not (sorted(mins) != list(mins) or sorted(maxes) != list(maxes)
                      or any(a > b for a, b in zip(mins, maxes)))
        # dask_logger.warning("Dataframe with unknown divisions; is_sorted %s", is_sorted)
        # print("IS SORTED:", is_sorted)
        if ddf.index.name is not None:
            ddf = ddf.reset_index().set_index(ddf.index.name, sorted=is_sorted)
        else:
            ddf = ddf.reset_index().set_index("index", sorted=is_sorted)
    
    index_counts = ddf.map_partitions(lambda _df: _df.index.value_counts(dropna=False).sort_index(),
                                      meta=int).compute()
    dask_logger.warning("Index_counts: %s", sum(index_counts))
    index = np.repeat(index_counts.index, index_counts.values)
    if target_memory is not None:
        mem_usage = ddf.memory_usage(deep=True).sum().compute()
        npartitions = 1 + mem_usage // target_memory
        # dask_logger.warning("Using %s for storing %s, target: %s", npartitions, mem_usage, target_memory)
    elif npartitions is None:
        npartitions=ddf.npartitions
    divisions, _ = dd.io.io.sorted_division_locations(index, npartitions=npartitions)
    try:
        repartitioned = ddf.repartition(divisions=divisions, force=True)
    except ValueError:
        try:
            repartitioned = ddf.repartition(npartitions=npartitions)
        except ValueError:
            repartitioned = ddf

    repartitioned = repartitioned.persist()
    _ = repartitioned.head(npartitions=-1, n=1)
    new_shape = repartitioned.shape[0].compute()
    if new_shape != orig_shape:
        # Something went wrong.
        dask_logger.critical("Error in rebalancing, old shape: %s; new shape: %s", orig_shape, new_shape)
        if orig_shape < new_shape:
            dask_logger.critical("Added new rows:")
            new_rows = repartitioned.index.compute().difference(ddf.index.compute())
            dask_logger.critical(repartitioned.loc[new_rows].head(npartitions=-1, n=10))
        else:
            dask_logger.critical("Lost rows:")
            lost_rows = ddf.index.compute().difference(repartitioned.index.compute())
            dask_logger.critical(lost_rows)
            dask_logger.critical(ddf.loc[lost_rows].head(npartitions=-1, n=10))
        raise AssertionError

    assert new_shape == orig_shape, (orig_shape, new_shape)
    return repartitioned


def trim_fai(fai: typing.Union[str, dd.DataFrame]) -> dd.DataFrame:
    if isinstance(fai, str):
        new_fai = dd.read_parquet(fai, infer_divisions=True)
    else:
        assert isinstance(fai, dd.DataFrame)
        new_fai = fai[:]
    if "derived_from_split" not in new_fai.columns:  # Not initialised
        new_fai["derived_from_split"] = False
        new_fai["previous_iteration"] = new_fai.index
        new_fai = new_fai.persist()
        return new_fai

    previous = new_fai["previous_iteration"].unique().values.compute()
    original = new_fai.query("derived_from_split == False")
    original_indices = original.index.compute()
    midpoints = set.difference(set(previous), set(original_indices))
    derived = new_fai.query("derived_from_split == True")["orig_scaffold_index"].values.compute()
    all_indices = new_fai.index.compute()

    final = new_fai.loc[all_indices.difference(midpoints).difference(derived).values]
    return final.persist()
