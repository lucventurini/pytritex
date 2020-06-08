import pandas as pd
import numpy as np
from functools import partial
import dask.dataframe as dd
from typing import Union
import scipy.stats as stats
import dask.array as da


def _ssorted(row, by, how):
    if np.isnan(row["by"]).all():
        return np.nan
    res = np.searchsorted(row["by"], row[by], side="left")
    if res < len(row["by"]) and row["by"][res] == row[by]:
        return res
    if how == "right":
        res -= 1
        if res < 0:
            return np.nan
        else:
            return res
    elif how == "left":
        if res == len(row["by"]):
            return np.nan
        else:
            return res


def rank(series: Union[dd.Series, pd.Series, np.array]):
    return stats.rankdata(series) - 1


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
    is_dd = isinstance(left, dd.DataFrame)

    # Assign the rank to each "by" value
    if is_dd:
        grouped = left[[by]].compute().groupby(on)
        left["__idx_pos"] = dd.from_array(grouped[by].transform(rank).astype(int))
        s = grouped[by].agg(by=(by, lambda col: col.values.tolist()))
        s = dd.from_pandas(s, chunksize=1000)
    else:
        grouped = left.groupby(on)
        left["__idx_pos"] = grouped[by].transform(rank)
        s = grouped[by].agg(by=(by, lambda col: col.values.tolist()))

    merged = s.merge(right, how="outer", left_index=True, right_index=True)  # .set_index(on)
    assert merged.index.name == on
    _ss = partial(_ssorted, by=by, how=how)
    merged = merged[~merged[by].isna()]
    if is_dd:
        assert "by" in merged.columns
        shape = merged.shape[0].compute()
        chunks = tuple(merged.map_partitions(len).compute())
        # size = sum(chunks)
        __idx_pos = da.from_array(merged[["by", by]].apply(_ss, axis=1, meta=int).compute(),
                                  chunks=chunks)
        merged = merged.assign(__idx_pos=__idx_pos)
        assert shape == merged.shape[0].compute(), (shape, merged.shape[0].compute())
    else:
        # bait = ~merged["by"].isna()
        merged["__idx_pos"] = merged[["by", by]].apply(_ss, axis=1)

    merged = merged.drop("by", axis=1)
    if is_dd:
        res = dd.merge(left, merged,
                       on=[on, "__idx_pos"], how="right", suffixes=["_y", ""]).drop(
            "__idx_pos", axis=1).drop(by + "_y", axis=1)
        res = res.reset_index(drop=False)
    else:
        res = left.drop(by, axis=1).merge(merged,
                                          on=[on, "__idx_pos"], how="right").drop("__idx_pos", axis=1)
    res = res[~res[by].isna()]
    return res
