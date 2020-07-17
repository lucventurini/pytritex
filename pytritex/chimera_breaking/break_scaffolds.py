# import pandas as pd
import numpy as np
from .transposition import _transpose_cssaln, _transpose_molecules, _transpose_fpairs
from .transposition import _transpose_hic_cov, _transpose_molecule_cov
from .calculate_broken_scaffolds import calculate_broken_scaffolds
from ..anchoring import anchor_scaffolds
from time import ctime
import dask.dataframe as dd
import os
from joblib import Memory
from dask.distributed import Client
from ..utils import assign_to_use_column
import logging
dask_logger = logging.getLogger("dask")


def break_scaffolds(breaks, client: Client, save_dir: str, memory: Memory,
                    assembly, slop, cores=1, species="wheat") -> dict:

    fai = assembly["fai"]
    if isinstance(fai, str):
        fai = dd.read_parquet(fai, infer_divisions=True)
    new_assembly = calculate_broken_scaffolds(breaks, fai=fai, slop=slop, save_dir=save_dir)
    for key in ["binsize", "innerDist", "minNbin", "popseq"]:
        new_assembly[key] = assembly[key]
    # Now extract the correct rows and columns from the FAI:
    # ie excluding the rows of the original scaffolds that have since been split.
    trimmed_fai = assign_to_use_column(new_assembly["fai"])
    new_indices = trimmed_fai.index.compute().difference(fai.index.values.compute()).values
    # old_indices = trimmed_fai[trimmed_fai.to_use == True].index.compute().intersection(
    #     fai.index.values.compute())

    dask_logger.debug("%s Transposing the CSS alignment", ctime())
    new_assembly["cssaln"] = _transpose_cssaln(assembly["cssaln"], trimmed_fai, save_dir)
    # Do the same with HiC pairs
    dask_logger.debug("%s Transposing the HiC alignment", ctime())
    new_assembly["fpairs"] = _transpose_fpairs(assembly.get("fpairs", None), trimmed_fai, save_dir)

    dask_logger.debug("%s Transposing the 10X alignment", ctime())
    new_assembly["molecules"] = _transpose_molecules(assembly.get("molecules", None), trimmed_fai, save_dir)

    dask_logger.debug("%s Anchoring the transposed data", ctime())
    # Extract the relevant for the anchoring.
    # fai_to_anchor = trimmed_fai.loc[new_indices]
    # partial_assembly = dict((key, value) for key, value in new_assembly.items()
    #                         if key != "fai")
    # partial_assembly["fai"] = fai_to_anchor
    # partial_assembly["info"] = dd.read_parquet(assembly["info"])
    new_assembly = anchor_scaffolds(new_assembly, scaffolds=new_indices,
                                    save=None, species=species, client=client)
    # Now let's recalculate the coverage
    info = new_assembly["info"]
    info = info.drop("mr_10x", axis=1, errors="ignore")
    info = info.drop("mr", axis=1, errors="ignore")
    info = info.drop("mri", axis=1, errors="ignore")
    new_assembly["info"] = info
    dask_logger.debug("%s Transposing the HiC coverage data", ctime())
    new_assembly = _transpose_hic_cov(
        new_assembly,
        fai=trimmed_fai,
        old_info=assembly["info"],
        scaffolds=new_indices,
        save_dir=save_dir,
        coverage=assembly["cov"], fpairs=assembly["fpairs"])
    dask_logger.debug("%s Transposing the 10X coverage data", ctime())
    new_assembly = _transpose_molecule_cov(
        new_assembly, fai=trimmed_fai, save_dir=save_dir, scaffolds=new_indices,
        client=client, assembly=assembly, cores=cores)
    # Change values that are 0 to nan
    dask_logger.debug("%s Transposed the 10X coverage data", ctime())
    for key in ["popseq_alphachr2", "popseq_alphachr", "sorted_alphachr", "sorted_alphachr2"]:
        new_assembly["info"][key] = new_assembly["info"][key].map({0: np.nan})
    new_assembly["info"]["derived_from_split"] = new_assembly["info"]["derived_from_split"].fillna(False)
    for key in ["molecule_cov", "info", "cov"]:
        # BugFix for pyarrow not handling float16 and int16
        dask_logger.debug("%s Storing %s", ctime(), key)
        if (new_assembly[key].dtypes == np.int16).any():
            for index, col in enumerate(new_assembly[key].dtypes == np.int16):
                if col is False:
                    continue
                col = new_assembly[key].dtypes.index[index]
                new_assembly[key][col] = new_assembly[key][col].astype(np.int32)
        if (new_assembly[key].dtypes == np.float16).any():
            for index, col in enumerate(new_assembly[key].dtypes == np.float16):
                if col is False:
                    continue
                col = new_assembly[key].dtypes.index[index]
                new_assembly[key][col] = new_assembly[key][col].astype(np.float32)
        fname = os.path.join(save_dir, key)
        dd.to_parquet(new_assembly[key], fname, compression="gzip", engine="pyarrow", compute=True)
        dask_logger.debug("%s Stored %s", ctime(), key)
        new_assembly[key] = fname

    print(ctime(), "Finished transposing")
    return new_assembly
