import pandas as pd
import numpy as np
import re


def rolling_join(left: pd.DataFrame, right: pd.DataFrame, on, by, right_key):
    """Implementation of the R data.table rolling join procedure, when the
    "by" key is numeric."""
    restore_right_index = False
    index_names = None
    if not isinstance(right.index, pd.RangeIndex) and right.index.names != [None]:
        restore_right_index = True
        index_names = right.index.names
        right.reset_index(drop=False, inplace=True)
    initial = pd.merge(left, right, on=on, suffixes=("_left", "_right"), how="right")
    initial.loc[:, "distance"] = (initial[by + "_left"] - initial[by + "_right"])
    if not isinstance(by, list):
        by = [by]
    if not isinstance(right_key, list):
        right_key = [right_key]
    retain_columns = []
    for column in initial.columns:
        keep = True
        if column in right:
            keep = False
        elif column == "distance":
            keep = False
        if keep:
            if re.sub("_right$", "", column) in by:
                keep = False
            elif re.sub("_left$", "", column) in by:
                keep = False
            elif column in left.columns:
                keep = True
        if keep:
            retain_columns.append(column)

    retain_columns = right_key + retain_columns
    renamer = dict((_, re.sub("_right$", "", _)) for _ in retain_columns if _.endswith("_right"))
    vals = initial.sort_values("distance", ascending=True).groupby(right_key, as_index=False).head(1).loc[:,
           retain_columns].rename(columns=renamer)
    right_cols = [_ for _ in right.columns]
    right = pd.merge(vals.set_index(right_key), right.loc[:, right_cols].set_index(right_key),
                     left_index=True, right_index=True, how="inner").reset_index(drop=False)
    if restore_right_index:
        right.set_index(index_names, inplace=True)
    return right
