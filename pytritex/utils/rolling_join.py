import pandas as pd
import numpy as np


def rolling_join(left: pd.DataFrame, right: pd.DataFrame, on, by):
    """Implementation of the R data.table rolling join procedure."""

    assert on in right.columns
    assert on in left.columns
    right = right.sort_values([on, by])
    left = left.copy().sort_values([on, by])
    left.loc[:, "__idx_pos"] = left.groupby(on)[by].transform(lambda s: np.arange(s.shape[0], dtype=np.int))
    left = left.set_index([on])

    def ssorted(row):
        if row[on] in left.index:
            return np.maximum(0, np.searchsorted(left.loc[row[on], by], row[by]) - 1)
        else:
            return np.nan

    assert "__idx_pos" in left.columns
    assert on in right.columns

    right.loc[:, "__idx_pos"] = right.apply(ssorted, axis=1)

    res = left.reset_index(drop=False).drop(by, axis=1).merge(
        right, on=[on, "__idx_pos"], how="right").drop("__idx_pos", axis=1)
    return res
