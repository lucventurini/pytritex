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
from ..utils import _rebalance_ddf
import logging
dask_logger = logging.getLogger("dask")


def trim_fai(fai: str) -> dd.DataFrame:
    new_fai = dd.read_parquet(fai, infer_divisions=True)
    _left1 = new_fai.query("derived_from_split == True")
    _left1_indices = _left1["orig_scaffold_index"].compute()
    # Now, let us get those that are NOT split AND do not have split children
    _left2 = new_fai.query("derived_from_split == False")
    _left2 = _left2.loc[_left2.index.compute().difference(_left1_indices).values]
    # Put the unbroken scaffolds first
    new_fai = dd.concat([_left2, _left1])
    return new_fai.persist()


def break_scaffolds(breaks, client: Client, save_dir: str, memory: Memory,
                    assembly, slop, cores=1, species="wheat") -> dict:

    fai = assembly["fai"]
    new_assembly = memory.cache(calculate_broken_scaffolds)(breaks, fai=fai, slop=slop, save_dir=save_dir)
    for key in ["binsize", "innerDist", "minNbin", "popseq"]:
        new_assembly[key] = assembly[key]
    # Now extract the correct rows and columns from the FAI:
    # ie excluding the rows of the original scaffolds that have since been split.
    trimmed_fai = trim_fai(new_assembly["fai"])
    # new_shape = new_fai["orig_scaffold_index"].unique().shape[0].compute()
    # old_shape = fai["orig_scaffold_index"].unique().shape[0].compute()
    # assert new_shape == old_shape, (new_shape, old_shape)
    # Remove broken scaffolds from CSS alignment, put correct ones
    dask_logger.warning("%s Transposing the CSS alignment", ctime())
    new_assembly["cssaln"] = memory.cache(_transpose_cssaln)(
        assembly["cssaln"], trimmed_fai, save_dir)
    # Do the same with HiC pairs
    dask_logger.warning("%s Transposing the HiC alignment", ctime())
    new_assembly["fpairs"] = memory.cache(_transpose_fpairs)(
        assembly.get("fpairs", None), trimmed_fai, save_dir)

    dask_logger.warning("%s Transposing the 10X alignment", ctime())
    new_assembly["molecules"] = memory.cache(_transpose_molecules)(
        assembly.get("molecules", None), trimmed_fai, save_dir)

    dask_logger.warning("%s Anchoring the transposed data", ctime())
    new_assembly = anchor_scaffolds(new_assembly, save=save_dir, species=species, client=client)
    # Now let's recalculate the coverage
    info = new_assembly["info"]
    if isinstance(info, str):
        info = dd.read_parquet(info)
    info = info.drop("mr_10x", axis=1, errors="ignore")
    info = info.drop("mr", axis=1, errors="ignore")
    info = info.drop("mri", axis=1, errors="ignore")
    new_assembly["info"] = info
    dask_logger.warning("%s Transposing the HiC coverage data", ctime())
    new_assembly = memory.cache(_transpose_hic_cov, ignore=["cores", "client"])(
        new_assembly, assembly["info"], fai=trimmed_fai,
        save_dir=save_dir,
        coverage=assembly["cov"], fpairs=assembly["fpairs"], cores=cores, client=client)
    dask_logger.warning("%s Transposing the 10X coverage data", ctime())
    new_assembly = memory.cache(_transpose_molecule_cov, ignore=["client", "cores"])(
        new_assembly, trimmed_fai, save_dir=save_dir,
        client=client, assembly=assembly, cores=cores)
    # Change values that are 0 to nan
    for key in ["popseq_alphachr2", "popseq_alphachr", "sorted_alphachr", "sorted_alphachr2"]:
        new_assembly["info"][key] = new_assembly["info"][key].map({0: np.nan})
    new_assembly["info"]["derived_from_split"] = new_assembly["info"]["derived_from_split"].fillna(False)
    for key in ["molecule_cov", "info", "cov"]:
        # BugFix for pyarrow not handling float16 and int16
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

        new_assembly[key] = _rebalance_ddf(new_assembly[key], target_memory=5 * 10**7)
        fname = os.path.join(save_dir, key)
        dd.to_parquet(new_assembly[key], fname, compression="gzip", engine="pyarrow", compute=True)
        new_assembly[key] = fname

    print(ctime(), "Finished transposing")
    return new_assembly
