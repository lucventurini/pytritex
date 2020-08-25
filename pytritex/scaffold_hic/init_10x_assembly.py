import dask.dataframe as dd
import numpy as np
import pandas as pd
from ..scaffold_10x.make_agp import make_agp


def init_10x_assembly(assembly, map_10x, gap_size=100, molecules=False):

    """Function to start the process of using HiC data for 10X-super-scaffolded assemblies.
    :param assembly: the initial assembly from *before* the 10X super-scaffolding.
    :param map_10x: the super-scaffolded assembly
    :param molecules: boolean flag
    :param gap_size: length of gaps in the AGP
    """

    # super <- map_10x
    #
    #  copy(assembly$cssaln)->z
    #  z[, orig_scaffold := NULL]
    #  z[, orig_scaffold_length := NULL]
    #  z[, orig_pos := NULL]
    #  z[, scaffold_length := NULL]
    #  super$agp[, .(scaffold, super, super_start, super_end, orientation)][z, on="scaffold"]->z
    #  z[orientation == 1, pos := super_start - 1 + pos]
    #  z[orientation == -1, pos := super_end - pos + 1]
    #  z[, scaffold := NULL]
    #  setnames(z, "super", "scaffold")
    #  z[, c("super_start", "super_end", "orientation") := list(NULL, NULL, NULL)]
    #  super$info[, .(scaffold=super, scaffold_length=length)][z, on="scaffold"]->z
    #  z->s_cssaln

    map_10x["agp"] = make_agp(map_10x["membership"], info=assembly["info_10x"], gap_size=gap_size)["agp"]
    cssaln = assembly["cssaln"]
    if isinstance(cssaln, str):
        cssaln = dd.read_parquet(cssaln, infer_divisions=True)
    cssaln = cssaln.drop(["orig_scaffold", "orig_scaffold_length", "orig_pos", "scaffold_length"], axis=1)
    cssaln = dd.merge(map_10x["agp"][["super", "super_start", "super_end", "orientation"]],
                      cssaln, on="scaffold_index")
    # cssaln.
