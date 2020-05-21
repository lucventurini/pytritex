import pandas as pd
import numpy as np
from functools import partial


def _ssorted(row, by, how):
    res = np.searchsorted(row["by"], row[by], side="left")
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


def rolling_join(left: pd.DataFrame, right: pd.DataFrame, on, by, how="right"):
    """Implementation of the R data.table rolling join procedure."""

    assert on in right.columns
    assert on in left.columns
    # right = right.sort_values([on, by])
    left = left.copy().sort_values([on, by])
    grouped = left.groupby(on)
    s = grouped[by].agg(by=("by", lambda col: col.values.tolist()))
    left.loc[:, "__idx_pos"] = grouped[by].transform(lambda s: np.arange(s.shape[0], dtype=np.int))
    merged = s.merge(right, how="outer",
                     left_index=True, right_on=on).sort_values([on, by])
    merged = merged[~merged[by].isna()]
    bait = ~merged["by"].isna()
    _ss = partial(_ssorted, by=by, how=how)
    merged.loc[bait, "__idx_pos"] = merged.loc[bait, ["by", by]].apply(_ss, axis=1)
    res = left.drop(by, axis=1).merge(merged.drop("by", axis=1),
                                      on=[on, "__idx_pos"], how="right").drop("__idx_pos", axis=1)
    res = res[~res[by].isna()]
    return res
