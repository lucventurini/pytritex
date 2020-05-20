import pandas as pd
import numpy as np


def rolling_join(left: pd.DataFrame, right: pd.DataFrame, on, by):
    """Implementation of the R data.table rolling join procedure."""

    _cp_right = right.sort_values([on, by])
    left = left.sort_values([on, by])
    groups = left.groupby([on])
    left.loc[:, "__idx_pos"] = groups[by].transform(lambda s: np.arange(s.shape[0], dtype=np.int))
    assert "__idx_pos" in left.columns

    for on_key, right_indices in iter(_cp_right.groupby([on]).groups.items()):
        if on_key in groups.groups.keys():
            _cp_right.loc[right_indices, "__idx_pos"] = np.maximum(
                0, np.searchsorted(groups.get_group(on_key)[by], _cp_right.loc[right_indices, by]) - 1)
        else:
            _cp_right.loc[right_indices, "__idx_pos"] = np.nan

    res = left.drop(by, axis=1).merge(_cp_right, on=[on, "__idx_pos"], how="right").drop("__idx_pos", axis=1)
    # assert res.loc[res[by].isna()].shape[0] <= right.loc[right[by].isna()].shape[0]
    return res

