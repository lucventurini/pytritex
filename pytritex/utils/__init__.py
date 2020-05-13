import pandas as pd
import numpy as np
import itertools


def n50(l: np.array, p=0.5):
    l.sort()
    return l[np.where(l.cumsum() >= (1 - p) * l.sum())[0][0]]
