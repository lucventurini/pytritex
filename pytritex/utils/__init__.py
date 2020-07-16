import pandas as pd
import numpy as np
from .rolling_join import rolling_join
import re
import dask.dataframe as dd
import dask.utils
from dask.base import compute
import itertools as it
import typing
import time
import logging
dask_logger = logging.getLogger("dask")


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
    # First create a list for each of the grouped chunks
    def internal(c):
        if (c != c).all():
            return [np.nan]
        f = [_ for _ in c if _ == _]
        f = [_ if isinstance(_, list) else [_] for _ in f]
        return list(it.chain.from_iterable(f))
    chunks = grouped.apply(internal)
    return chunks


def _agg_third(chunks):
    chunks = chunks.apply(lambda s: np.nan if len(s) == 1 else s[1])
    return chunks


second_agg = dd.Aggregation("second", lambda s: s.apply(list), _agg_second, _agg_third)


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


def assign_to_use_column(fai: typing.Union[str, dd.DataFrame]) -> dd.DataFrame:
    if isinstance(fai, str):
        new_fai = dd.read_parquet(fai, infer_divisions=True)
    else:
        assert isinstance(fai, dd.DataFrame)
        new_fai = fai[:]
    if "derived_from_split" not in new_fai.columns:  # Not initialised
        new_fai["derived_from_split"] = False
        new_fai["previous_iteration"] = new_fai.index
        new_fai["to_use"] = True
        return new_fai

    previous = new_fai["previous_iteration"].unique().values.compute()
    original = new_fai.query("derived_from_split == False")
    original_indices = original.index.compute()
    midpoints = set.difference(set(previous), set(original_indices))
    derived = new_fai.query("derived_from_split == True")["orig_scaffold_index"].values.compute()
    all_indices = new_fai.index.compute()
    new_fai["to_use"] = new_fai.index.isin(all_indices.difference(midpoints).difference(derived).values)
    return new_fai
