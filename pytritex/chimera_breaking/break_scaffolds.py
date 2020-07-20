import numpy as np
from .transposition import _transpose_cssaln, _transpose_molecules, _transpose_fpairs
from ..sequencing_coverage import add_hic_cov, add_molecule_cov
from .calculate_broken_scaffolds import calculate_broken_scaffolds
from ..anchoring import anchor_scaffolds
from time import ctime
import dask.dataframe as dd
from dask.distributed import Client
import logging
dask_logger = logging.getLogger("dask")
import time


def break_scaffolds(breaks, client: Client, save_dir: str,
                    assembly, slop, cores=1, species="wheat") -> dict:

    dask_logger.debug("%s Calculating broken scaffolds", time.ctime())
    fai, fai_name, new_indices = calculate_broken_scaffolds(
        breaks, fai=assembly["fai"], slop=slop, save_dir=save_dir, cores=cores)
    new_assembly = {"fai": fai_name}
    for key in ["binsize", "innerDist", "minNbin", "popseq", "mol_binsize"]:
        new_assembly[key] = assembly[key]
    # Now extract the correct rows and columns from the FAI:
    # ie excluding the rows of the original scaffolds that have since been split.
    dask_logger.debug("%s Transposing the CSS alignment", ctime())
    new_assembly["cssaln"] = _transpose_cssaln(assembly["cssaln"], fai, save_dir)
    # Do the same with HiC pairs
    dask_logger.debug("%s Transposing the HiC alignment", ctime())
    new_assembly["fpairs"] = _transpose_fpairs(assembly.get("fpairs", None), fai, save_dir)

    dask_logger.debug("%s Transposing the 10X alignment", ctime())
    new_assembly["molecules"] = _transpose_molecules(assembly.get("molecules", None), fai, save_dir)

    dask_logger.debug("%s Anchoring the transposed data", ctime())
    # Extract the relevant for the anchoring.
    # fai_to_anchor = fai.loc[new_indices]
    # partial_assembly = dict((key, value) for key, value in new_assembly.items()
    #                         if key != "fai")
    # partial_assembly["fai"] = fai_to_anchor
    # partial_assembly["info"] = dd.read_parquet(assembly["info"])
    new_assembly = anchor_scaffolds(new_assembly, scaffolds=new_indices,
                                    save=save_dir, species=species, client=client)
    # Now let's recalculate the coverage
    info = dd.read_parquet(new_assembly["info"], infer_divisions=True)
    info = info.drop("mr_10x", axis=1, errors="ignore")
    info = info.drop("mr", axis=1, errors="ignore")
    info = info.drop("mri", axis=1, errors="ignore")
    new_assembly["info"] = info
    dask_logger.debug("%s Transposing the 10X coverage data", ctime())
    new_assembly = add_molecule_cov(new_assembly, save_dir=save_dir, binsize=new_assembly["mol_binsize"])
    dask_logger.debug("%s Transposing the HiC coverage data", ctime())
    new_assembly = add_hic_cov(new_assembly,
                               save_dir=save_dir,
                               binsize=new_assembly["binsize"], minNbin=new_assembly["minNbin"],
                               innerDist=new_assembly["innerDist"])
    # Change values that are 0 to nan
    dask_logger.debug("%s Transposed the 10X coverage data", ctime())
    info_name = new_assembly["info"][:]
    info = dd.read_parquet(new_assembly["info"], infer_divisions=True)
    for key in ["popseq_alphachr2", "popseq_alphachr", "sorted_alphachr", "sorted_alphachr2"]:
        info[key] = info[key].map({0: np.nan})
    info["derived_from_split"] = info["derived_from_split"].fillna(False)
    dd.to_parquet(info, info_name, engine="pyarrow", compute=True, compression="gzip")
    dask_logger.debug("%s Finished transposing", time.ctime())
    return new_assembly
