from .break_scaffolds import break_scaffolds
from .find_10x_breaks import find_10x_breaks
from functools import partial
import dask.dataframe as dd
from dask.distributed import Client
import logging
dask_logger = logging.getLogger("dask")
import time
import os
from shutil import rmtree
from joblib import Memory


def break_10x(assembly: dict, save_dir: str, client: Client,
              memory: Memory,
              species="wheat", ratio=-3, interval=5e4, minNbin=20,
              dist=2e3, slop=1e3, intermediate=False, cores=1, maxcycle=float("inf")):
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
    molecule_cov = assembly.get("molecule_cov", None)
    cycle = 0
    dask_logger.warning("%s Finding initial breaks", time.ctime())
    breaks = find_10x_breaks(molecule_cov, interval=interval, minNbin=minNbin, dist=dist, ratio=ratio)
    dask_logger.warning("%s Found initial breaks", time.ctime())
    assemblies = dict()
    if intermediate is True:
        assemblies = {0: assembly}

    orig_assembly = assembly.copy()
    base = os.path.join(save_dir, "joblib", "pytritex", "chimera_breaking")
    while breaks is not None and breaks.shape[0] > 0 and cycle <= maxcycle:
        cycle += 1
        save_dir = os.path.join(base, str(cycle))
        dask_logger.warning("%s Starting cycle %s of %s, with %s breaks",
                            time.ctime(), cycle, maxcycle, breaks.shape[0])
        orig_assembly["fai"] = assembly["fai"]
        assembly = memory.cache(break_scaffolds, ignore=["client", "cores"])(
            breaks=breaks, save_dir=save_dir,
            client=client, assembly=orig_assembly, slop=slop, cores=cores, species=species)
        if intermediate is True:
            assemblies[cycle] = assembly
        dask_logger.warning("%s Finished cycle %s of %s; finding new breaks",
                            time.ctime(), cycle, maxcycle)
        breaks = memory.cache(find_10x_breaks)(assembly.get("molecule_cov", None),
                                               interval=interval, minNbin=minNbin, dist=dist, ratio=ratio)
        dask_logger.warning("%s Finished cycle %s of %s; new breaks: %s",
                            time.ctime(), cycle, maxcycle, None if breaks is None else breaks.shape[0])

    if intermediate is False:
        return {"assembly": assembly}
    else:
        return {"assemblies": assemblies}
