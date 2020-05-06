import pandas as pd
import numpy as np


def find_10x_breaks(assembly, scaffolds=None, interval = 5e4, minNbin = 20, dist = 5e3, ratio = -3):
    cov = assembly["molecule_cov"][:]
    if scaffolds is not None:
        cov = pd.merge(cov, scaffolds, left_on="scaffold", right_on="scaffold")
    cov["b"] = cov["bin"] // interval * interval
    bait = (cov.nbin >= minNbin) & (cov.r <= ratio) & (np.minimum(cov.bin, cov.length - bin) >= dist)
    e = cov.loc[bait]
    if e.shape[0] == 0:
        return None
    # Order according to "r", group by scaffold and "b", then take the best (ie first)
    e.sort_values(by=["r"]).groupby(["scaffold", "b"]).head(1)
    e.rename(columns={"bin", "br"}, inplace=True)  # br: break point
    return e
