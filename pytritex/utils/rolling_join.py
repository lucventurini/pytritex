import pandas as pd
import numpy as np


def rolling_join(left: pd.DataFrame, right: pd.DataFrame, on, by, how="left"):
    """Implementation of the R data.table rolling join procedure."""

    def nearest(array, value, how="left"):
        idx = np.searchsorted(array, value, side="left")
        if idx == array.shape[0]:  # It would be beyond the data
            idx -= 1
        if how == "left":
            pass
        elif how == "right":
            idx += 1
        return idx

    _cp_right = right.sort_values([on, by])
    groups = _cp_right.groupby([on])
    _cp_right.loc[:, "__idx_pos"] = groups[by].transform(lambda s: np.arange(s.shape[0], dtype=np.int))
    hashed_groups = dict((key, groups.get_group(key)) for key in groups.groups)
    left.loc[:, "__idx_pos"] = left.apply(lambda row: np.nan if row[on] not in groups.groups else
                                          nearest(hashed_groups[row[on]][by], row[by], how=how), axis=1)
    res = left.drop(by, axis=1).merge(_cp_right, on=[on, "__idx_pos"], how="right").drop("__idx_pos", axis=1)
    assert res.loc[res[by].isna()].shape[0] <= right.loc[right[by].isna()].shape[0]
    return res
