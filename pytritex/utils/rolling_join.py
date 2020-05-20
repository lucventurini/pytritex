import pandas as pd
import numpy as np


def rolling_join(left: pd.DataFrame, right: pd.DataFrame, on, by, how="right"):
    """Implementation of the R data.table rolling join procedure."""

    assert on in right.columns
    assert on in left.columns
    # right = right.sort_values([on, by])
    left = left.copy().sort_values([on, by])
    grouped = left.groupby(on)
    left.loc[:, "__idx_pos"] = grouped[by].transform(lambda s: np.arange(s.shape[0], dtype=np.int))
    hashed = dict((_[0], _[1].values) for _ in iter(grouped[by]))

    def ssorted(row):
        if row[on] in hashed:
            l = hashed[row[on]]
            print(l.shape)
            if l.shape[0] == 1:
                if how == "right":
                    if l[0] > row[by]:
                        return np.nan
                    else:
                        return 0
                elif how == "left":
                    if l[0] < row[by]:
                        return np.nan
                    else:
                        return 0
            else:
                print(how)
                res = np.searchsorted(l, row[by], side="left")
                if how == "right":
                    res -= 1
                    if res < 0:
                        return np.nan
                    else:
                        return res
                elif how == "left":
                    if res == left.loc[row[on], "__idx_pos"].max():
                        return np.nan
                    else:
                        return res
        else:
            return np.nan

    assert "__idx_pos" in left.columns
    assert on in right.columns

    right.loc[:, "__idx_pos"] = right[[on, by]].apply(ssorted, axis=1)

    res = left.drop(by, axis=1).merge(
        right, on=[on, "__idx_pos"], how="right").drop("__idx_pos", axis=1)
    # assert res.loc[res[by].isna()].shape[0] <= right.loc[right[by].isna()].shape[0]
    return res
