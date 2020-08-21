import pandas as pd
import numpy as np
from functools import partial
import dask.dataframe as dd
from typing import Union
import scipy.stats as stats
import dask.array as da


def _ssorted(row, by, how):
    if np.isnan(row["by"]).all():
        raise ValueError(row)
    res = np.searchsorted(row["by"], row[by], side="left")
    if res < len(row["by"]) and row["by"][res] == row[by]:
        return res
    if how == "right":
        return max(res - 1, 0)
    elif how == "left":
        return min(res, len(row["by"]) - 1)


def rank(series: Union[dd.Series, pd.Series, np.array]):
    return stats.rankdata(series) - 1
    # return np.argsort(series)


def sort_column(col):
    return np.sort(col.values).tolist()


def rolling_join(left: Union[pd.DataFrame, dd.DataFrame], right: Union[pd.DataFrame, dd.DataFrame],
                 on, by, how="right") -> Union[pd.DataFrame, dd.DataFrame]:
    """Implementation of the R data.table rolling join procedure."""

    assert on in right.columns or right.index.name == on
    if right.index.name != on:
        drop = (right.index.name is None)
        right = right.reset_index(drop=drop).set_index(on)

    left = left.copy()
    assert on in left.columns or left.index.name == on
    if left.index.name != on:
        drop = (left.index.name is None)
        left = left.reset_index(drop=drop).set_index(on)
    # right = right.sort_values([on, by])
    is_dd = isinstance(left, dd.DataFrame) or isinstance(right, dd.DataFrame)

    # Assign the rank to each "by" value
    if is_dd:
        assert left.shape[0].compute() > 0
        if not isinstance(left, dd.DataFrame):
            left = dd.from_pandas(left, npartitions=100)
        if not isinstance(right, dd.DataFrame):
            right = dd.from_pandas(right, npartitions=100)
        grouped = left[[by]].compute().groupby(on)
        chunks = tuple(left.map_partitions(len).compute())
        left["__idx_pos"] = da.from_array(grouped[by].transform(rank).astype(int),
                                          chunks=chunks)
        assert left.shape[0].compute() > 0
        s = grouped.agg(by=pd.NamedAgg(column=by, aggfunc=sort_column))[["by"]]
        s = dd.from_pandas(s, chunksize=1000)
        # assert s.index.dtype == left.index.dtype, (s.index.)
    else:
        assert left.shape[0] > 0
        grouped = left.groupby(on)
        left["__idx_pos"] = grouped[by].transform(rank)
        s = grouped[by].agg(by=(by, lambda col: np.sort(col.values).tolist() ))

    merged = dd.merge(s, right, how="outer", on=on)
    assert pd.api.types.is_numeric_dtype(merged.index.dtype) is True, np.unique(merged.index.values.compute())
    assert merged.index.name == on
    _ss = partial(_ssorted, by=by, how=how)
    merged = merged[~merged[by].isna()]
    if is_dd:
        assert "by" in merged.columns
        assert left.shape[0].compute() > 0
        shape = merged.shape[0].compute()
        assert shape > 0
        chunks = tuple(merged.map_partitions(len).compute())
        # size = sum(chunks)
        __idx_pos = da.from_array(
            merged[["by", by]].apply(_ss, axis=1, meta=int).compute(),
            chunks=chunks)
        merged = merged.assign(__idx_pos=__idx_pos)
        assert shape == merged.shape[0].compute(), (shape, merged.shape[0].compute())
    else:
        # bait = ~merged["by"].isna()
        merged["__idx_pos"] = merged[["by", by]].apply(_ss, axis=1)

    merged = merged.drop("by", axis=1)
    if is_dd:
        assert left.shape[0].compute() > 0
        assert merged.shape[0].compute() > 0
        try:
            res = dd.merge(left, merged, on=[on, "__idx_pos"], how="right", suffixes=["_y", ""])
        except ValueError:
            raise ValueError((left.index.dtype, merged.index.dtype))
        res = res.drop("__idx_pos", axis=1).drop(by + "_y", axis=1)
        res = res.reset_index(drop=False)
    else:
        res = left.drop(by, axis=1).merge(merged, on=[on, "__idx_pos"], how="right")
        res = res.drop("__idx_pos", axis=1)
    res = res[~res[by].isna()]
    return res
