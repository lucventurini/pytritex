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

    map_10x["agp"] = make_agp(map_10x["membership"], info=assembly["info_10x"], gap_size=gap_size)["agp"]
    cssaln = assembly["cssaln"]
    if isinstance(cssaln, str):
        cssaln = dd.read_parquet(cssaln, infer_divisions=True)
