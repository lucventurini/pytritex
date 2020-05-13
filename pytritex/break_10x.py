import pandas as pd
import numpy as np
from .find_10x_breaks import find_10x_breaks
from .break_scaffolds import break_scaffolds


def break_10x(assembly: dict, prefix="scaffold_corrected",
              breaks=None,
              species="wheat", ratio=-3, interval=5e4, minNbin=20, dist=2e3, slop=1e3,
              intermediate=False, ncores=1, maxcycle=float("inf")):

    if (dist <= 2 * slop):
        dist = 2 * slop + 1

    if breaks is None:
        breaks = find_10x_breaks(assembly, interval=interval, minNbin=minNbin, dist=dist, ratio=ratio)
    i = 0
    lbreaks = dict()
    assemblies = dict()
    if intermediate is True:
        assemblies[i] = assembly

    while breaks.shape[0] > 0 and i <= maxcycle:
        print("Cycle {i}:".format(i=i), breaks.shape[0], "break points detected")
        i += 1
        assembly = break_scaffolds(breaks, assembly=assembly, species=species, prefix="{}_".format(prefix),
                                   slop=slop, cores=ncores)
        if intermediate is True:
            assemblies[i] = assembly
        breaks = find_10x_breaks(assembly, interval=interval, minNbin=minNbin, dist=dist, ratio=ratio)
        if breaks is None or breaks.shape[0] == 0:
            break
        lbreaks[i] = breaks

    if intermediate is False:
        return {"assembly": assembly, "breaks": lbreaks}
    else:
        return {"assemblies": assemblies, "breaks": lbreaks}
