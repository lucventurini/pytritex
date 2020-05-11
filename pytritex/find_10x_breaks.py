import pandas as pd
import numpy as np


def find_10x_breaks(assembly: dict, scaffolds=None, interval=5e4, minNbin=20, dist=5e3, ratio=-3):
    # {
    #  cov <- copy(assembly$molecule_cov)
    #  if(!is.null(scaffolds)){
    #   cov[scaffolds, on="scaffold"]->cov
    #  }
    #  cov[, b := bin %/% interval * interval]
    #  cov[nbin >= minNbin & pmin(bin, length - bin) >= dist & r <= ratio]->e
    #  if(nrow(e) == 0){
    #   return(NULL)
    #  }
    #  e[order(r)][, idx := 1:.N, by=.(scaffold, b)][idx == 1]->e
    #  setnames(e, "bin", "br")[]
    #  e[]
    # }
    cov = assembly["molecule_cov"].copy()
    if scaffolds is not None:
        cov = cov.merge(scaffolds, on="scaffold", how="right")
    cov.loc[:, "b"] = cov["bin"] // interval * interval
    bait = (cov["nbin"] >= minNbin) & (np.minimum(cov["bin"], cov["length"] - cov["bin"]) >= dist)
    bait &= (cov["r"] <= ratio)
    e = cov.loc[bait, :]
    if e.shape[0] == 0:
        return None
    e.sort_values("r").group_by(["scaffold", "b"]).head(1).rename(columns={"bin": "br"})
    return e.copy()
