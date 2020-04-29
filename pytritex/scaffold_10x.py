import pandas as pd
import numpy as np
import subprocess as sp


def make_agp(membership: pd.DataFrame, gap_size=100):
    z = membership[["scaffold", "length", "super", "bin", "orientation"]]
    z.sort_values(axis=0, inplace=True, by=["super", "bin", "length"],
                  ascending=[True, True, False])
    z["gap"] = False
    z["index"] =
    pd.concat(z,
              pd.DataFrame(
                  {"scaffold": ["gap"], "gap": [True], "bin": [pd.NA],
                   "super": }

              )

              )


make_agp<-function(membership, gap_size=100){

  membership[, .(scaffold, length, super, bin, orientation)]->z

  setorder(z, super, bin, -length)
  z[, index := 2*1:.N-1]
  z[, gap := F]
  rbind(z, data.table(
      scaffold="gap", gap=T, super=z$super, bin=NA, length = gap_size,orientation = NA, index=z$index+1))->z
  z[order(index)][, head(.SD, .N-1), by=super]->z
  z[, n := .N, key=super]
  z[n > 1, super_start := cumsum(c(0, length[1:(.N-1)])) + 1, by = super]
  z[n == 1, super_start := 1]
  z[, super_end := cumsum(length), by = super]
  z[, n := NULL]

  z[, .(scaffold=scaffold, bed_start=0, bed_end=length, name=scaffold, score=1, strand=ifelse(is.na(orientation) | orientation == 1, "+", "-"), super=super)]->agp_bed

  list(agp=z, agp_bed=agp_bed)
 }