import pandas as pd
import numpy as np
import itertools


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


def unique_count(series: pd.Series):
    return np.unique(series).shape[0]
