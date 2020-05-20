import pandas as pd
import numpy as np
from .break_scaffolds import break_scaffolds
from .find_10x_breaks import find_10x_breaks
from functools import partial


def break_10x(assembly, species="wheat", ratio=-3, interval=5e4, minNbin=20,
              dist=2e3, slop=1e3, intermediate=False, cores=1, maxcycle=float("inf"), use_memory_fs=False):
    """Iteratively break scaffolds using 10X physical coverage. Proceed until no more breakpoints are found.
        :param assembly: dictionary containing data so far
        :param species: default wheat, species under analysis
        :param ratio: the maximum log2(ratio) between local and average 10X coverage to detect a break. Default -3
        :param interval:
        :param minNbin: minimum number of bins in scaffold to consider finding breaks
        :param dist:
        :param slop: when breaking scaffolds, this parameter determines the "slop", ie how much we should excise
                     in both directions from the low coverage point. Default 1kbps
        :param intermediate: boolean. If true, save intermediate steps.
        :param cores:
        :param maxcycle: how many times should we try to find and resolve breaks?
    """

    dist = max(dist, 2 * slop + 1)
    break_finder = partial(find_10x_breaks, interval=interval, minNbin=minNbin, dist=dist, ratio=ratio)
    breaks = break_finder(assembly.get("molecule_cov", None))
    cycle = 0
    lbreaks = {0: breaks}  # TODO: what is this for?
    if intermediate is True:
        assemblies = {0: assembly}

    while breaks.shape[0] > 0 and cycle <= maxcycle:
        cycle += 1
        assembly = break_scaffolds(breaks=breaks, assembly=assembly, slop=slop, cores=cores, species=species,
                                   use_memory_fs=use_memory_fs)
        if intermediate is True:
            assemblies[cycle] = assembly
        breaks = break_finder(assembly.get("molecule_cov", None))
        if breaks is None:
            break
        lbreaks[cycle] = breaks
        if intermediate is False:
            return {"assembly": assembly, "breaks": lbreaks}
        else:
            return {"assemblies": assemblies, "breaks": lbreaks}