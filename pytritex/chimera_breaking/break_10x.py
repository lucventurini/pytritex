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
    # if molecule_cov is not None:
    #     molecule_cov = dd.read_parquet(molecule_cov)
    breaks = find_10x_breaks(molecule_cov, interval=interval, minNbin=minNbin, dist=dist, ratio=ratio)
    cycle = 0
    # lbreaks = {0: breaks}  # TODO: what is this for?
    if intermediate is True:
        assemblies = {0: assembly}

    # for key in ['fai', 'cssaln', 'fpairs', 'molecules', 'info', 'molecule_cov', 'cov']:
    #     assembly[key] = dd.read_parquet(assembly[key], infer_divisions=True)

    base = os.path.join(save_dir, "joblib", "pytritex", "chimera_breaking")
    while breaks.shape[0].compute() > 0 and cycle <= maxcycle:
        cycle += 1
        save_dir = os.path.join(base, str(cycle))
        dask_logger.debug("%s Starting cycle %s of %s, with %s breaks",
                            time.ctime(), cycle, maxcycle, breaks.shape[0].compute())
        assembly = break_scaffolds(breaks=breaks, save_dir=save_dir,
                                   client=client,
                                   memory=memory,
                                   assembly=assembly, slop=slop, cores=cores, species=species)
        # dd.to_parquet(breaks, os.path.join(save_dir, "breaks"), compute=True,
        #               engine="pyarrow", compression="gzip")
        # for key in ['fai', 'cssaln', 'fpairs', 'molecules', 'info', 'molecule_cov', 'cov']:
        #     assembly[key] = dd.read_parquet(assembly[key], infer_divisions=True)
        if intermediate is True:
            assemblies[cycle] = assembly
        breaks = memory.cache(find_10x_breaks)(assembly.get("molecule_cov", None),
                                               interval=interval, minNbin=minNbin, dist=dist, ratio=ratio)
        if breaks is None:
            break
        # lbreaks[cycle] = os.path.join(save_dir, "breaks")

    # for key in ["fai", "cssaln", "molecules", "molecule_cov", "cov"]:
    #     name = os.path.join(base, key)
    #     dd.to_parquet(assembly[key], name,
    #                   compute=True, compression="gzip", engine="pyarrow")
    #     assembly[key] = name
    #
    # info_name = os.path.join(base, "anchored_css")
    # dd.to_parquet(assembly["info"], info_name, compute=True, compression="gzip")
    # assembly["info"] = info_name
    # fpairs_name = os.path.join(base, "anchored_hic_links")
    # dd.to_parquet(assembly["fpairs"], fpairs_name, compute=True)
    # assembly["fpairs"] = fpairs_name

    if intermediate is False:
        return {"assembly": assembly}   # , "breaks": lbreaks[cycle]}
    else:
        return {"assemblies": assemblies}   #, "breaks": lbreaks}
