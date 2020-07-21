from .find_10x_breaks import find_10x_breaks
import dask.dataframe as dd
from dask.distributed import Client
import logging
dask_logger = logging.getLogger("dask")
import time
import os
from shutil import rmtree
from ..sequencing_coverage.add_molecule_cov import add_10x_mr, add_molecule_cov
from ..sequencing_coverage.add_hic_cov import add_hic_cov
from ..anchoring import anchor_scaffolds
from .calculate_broken_scaffolds import calculate_broken_scaffolds
from .transposition.transpose_10x import _transpose_molecules
from .transposition.transpose_cssaln import _transpose_cssaln
from .transposition.transpose_hic import _transpose_fpairs
from joblib import Memory
import numpy as np


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
    fai_name = assembly["fai"][:]
    molecules = assembly["molecules"][:]
    molecule_cov = assembly["molecule_cov"][:]
    fai = dd.read_parquet(fai_name, infer_divisions=True)
    while breaks is not None and breaks.shape[0] > 0 and cycle <= maxcycle:
        cycle += 1
        save_dir = os.path.join(base, str(cycle))
        dask_logger.warning("%s Starting cycle %s of %s, with %s breaks",
                            time.ctime(), cycle, maxcycle, breaks.shape[0])
        dd.to_parquet(dd.from_pandas(breaks, chunksize=1e5), os.path.join(save_dir, "breaks"), compute=True,
                      compression="gzip", engine="pyarrow")
        fai, fai_name = calculate_broken_scaffolds(fai=fai_name, breaks=breaks, slop=slop, save_dir=save_dir, cores=cores)
        molecules = _transpose_molecules(orig_assembly["molecules"], fai, save_dir=save_dir)
        dask_logger.debug("%s Transposing the 10X alignment", time.ctime())
        molecule_cov = add_molecule_cov(
            assembly={"fai": fai, "molecules": molecules, "mol_binsize": assembly["mol_binsize"]},
            save_dir=save_dir, binsize=assembly["mol_binsize"], save_info=False)["molecule_cov"]
        breaks = memory.cache(find_10x_breaks)(molecule_cov,
                                               interval=interval, minNbin=minNbin, dist=dist, ratio=ratio)
        #
        # assembly = memory.cache(break_scaffolds, ignore=["client", "cores"])(
        #     breaks=breaks, save_dir=save_dir,
        #     client=client, assembly=orig_assembly, slop=slop, cores=cores, species=species)
        if intermediate is True:
            assemblies[cycle] = assembly
        dask_logger.warning("%s Finished cycle %s of %s; finding new breaks",
                            time.ctime(), cycle, maxcycle)

        dask_logger.warning("%s Finished cycle %s of %s; new breaks: %s",
                            time.ctime(), cycle, maxcycle, None if breaks is None else breaks.shape[0])

    # Now and only now let's do the rest
    new_assembly = dict()
    for key in ["binsize", "innerDist", "minNbin", "popseq", "mol_binsize"]:
        new_assembly[key] = assembly[key]
    new_assembly["molecules"] = molecules
    new_assembly["molecule_cov"] = molecule_cov
    new_assembly["fai"] = fai_name
    dask_logger.debug("%s Transposing the CSS alignment", time.ctime())
    new_assembly["cssaln"] = _transpose_cssaln(orig_assembly["cssaln"], fai, save_dir)
    dask_logger.debug("%s Transposing the HiC alignment", time.ctime())
    new_assembly["fpairs"] = _transpose_fpairs(orig_assembly.get("fpairs", None), fai, save_dir)
    dask_logger.debug("%s Anchoring the transposed data", time.ctime())
    new_indices = fai.query("to_use == True & derived_from_split == True").index.values.compute()
    new_assembly = anchor_scaffolds(new_assembly, scaffolds=new_indices,
                                    save=save_dir, species=species, client=client)
    info = dd.read_parquet(new_assembly["info"], infer_divisions=True)
    info = info.drop("mr_10x", axis=1, errors="ignore")
    info = info.drop("mr", axis=1, errors="ignore")
    info = info.drop("mri", axis=1, errors="ignore")
    info = add_10x_mr(info, dd.read_parquet(molecule_cov, infer_divisions=True))
    dd.to_parquet(info, os.path.join(save_dir, "info_10x"), compute=True, compression="gzip", engine="pyarrow")
    new_assembly["info"] = info
    dask_logger.debug("%s Transposing the HiC coverage data", time.ctime())
    new_assembly = add_hic_cov(new_assembly,
                               save_dir=save_dir,
                               binsize=new_assembly["binsize"], minNbin=new_assembly["minNbin"],
                               innerDist=new_assembly["innerDist"])
    info_name = new_assembly["info"][:]
    info = dd.read_parquet(new_assembly["info"], infer_divisions=True)
    for key in ["popseq_alphachr2", "popseq_alphachr", "sorted_alphachr", "sorted_alphachr2"]:
        info[key] = info[key].map({0: np.nan})
    info["derived_from_split"] = info["derived_from_split"].fillna(False)
    dd.to_parquet(info, info_name, engine="pyarrow", compute=True, compression="gzip")
    dask_logger.debug("%s Finished transposing", time.ctime())
    if "to_use" in fai.columns:
        fai = fai.query("to_use == True").drop("to_use", axis=1)
        dd.to_parquet(fai, fai_name, compute=True, engine="pyarrow", compression="gzip")
    new_assembly["fai"] = fai_name

    if intermediate is False:
        return {"assembly": new_assembly}
    else:
        assemblies[cycle] = new_assembly
        return {"assemblies": assemblies}
