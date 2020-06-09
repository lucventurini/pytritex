import pandas as pd
import numpy as np
from .transposition import _transpose_cssaln, _transpose_molecules, _transpose_fpairs
from .transposition import _transpose_hic_cov, _transpose_molecule_cov
from .calculate_broken_scaffolds import calculate_broken_scaffolds
from ..anchoring import anchor_scaffolds
from time import ctime
import dask.dataframe as dd
import os


def trim_fai(fai, save_dir):
    new_fai = dd.read_parquet(fai)
    _left1 = new_fai.query("derived_from_split == True")
    _left1_indices = _left1["orig_scaffold_index"].compute()
    # Now, let us get those that are NOT split AND do not have split children
    _left2 = new_fai.query("derived_from_split == False")
    _left2 = new_fai.query("scaffold_index not in @_left1_indices", local_dict=locals())
    # Put the unbroken scaffolds first
    new_fai = dd.concat([_left2, _left1])
    fname = os.path.join(save_dir, "trimmed_fai")
    dd.to_parquet(new_fai, fname,
                  compute=True, compression="gzip", engine="pyarrow")
    return fname


def break_scaffolds(breaks, memory, save_dir, assembly, slop, cores=1, species="wheat") -> dict:

    fai = assembly["fai"]
    base = os.path.join(save_dir, "joblib", "pytritex", "chimera_breaking")
    new_assembly = memory.cache(calculate_broken_scaffolds)(breaks, save_dir=base, fai=fai, slop=slop)
    for key in ["binsize", "innerDist", "minNbin", "popseq"]:
        new_assembly[key] = assembly[key]
    # Now extract the correct rows and columns from the FAI:
    # ie excluding the rows of the original scaffolds that have since been split.
    trimmed_fai = memory.cache(trim_fai)(new_assembly["fai"], base)
    # new_shape = new_fai["orig_scaffold_index"].unique().shape[0].compute()
    # old_shape = fai["orig_scaffold_index"].unique().shape[0].compute()
    # assert new_shape == old_shape, (new_shape, old_shape)
    # Remove broken scaffolds from CSS alignment, put correct ones
    print(ctime(), "Transposing the CSS alignment")
    new_assembly["cssaln"] = memory.cache(_transpose_cssaln)(assembly["cssaln"], trimmed_fai, base)
    # Do the same with HiC pairs
    print(ctime(), "Transposing the HiC alignments")
    new_assembly["fpairs"] = memory.cache(_transpose_fpairs)(assembly.get("fpairs", None), trimmed_fai, base)

    print(ctime(), "Transposing the 10X alignments")
    new_assembly["molecules"] = memory.cache(_transpose_molecules)(assembly.get("molecules", None),
                                                                   trimmed_fai, base)

    print(ctime(), "Anchoring the transposed data")
    new_assembly = anchor_scaffolds(new_assembly, species=species, save=base)
    # Now let's recalculate the coverage
    info = dd.read_parquet(new_assembly["info"])
    info = info.drop("mr_10x", axis=1, errors="ignore")
    info = info.drop("mr", axis=1, errors="ignore")
    info = info.drop("mri", axis=1, errors="ignore")
    new_assembly["info"] = info
    print(ctime(), "Transposing the HiC coverage data")
    new_assembly = _transpose_hic_cov(new_assembly, assembly["info"], trimmed_fai,
                                      assembly["cov"], assembly["fpairs"], cores=cores)
    print(ctime(), "Transposing the 10X coverage data")
    new_assembly = _transpose_molecule_cov(new_assembly, trimmed_fai, assembly, cores=cores)
    # Change values that are 0 to nan
    for key in ["popseq_alphachr2", "popseq_alphachr", "sorted_alphachr", "sorted_alphachr2"]:
        new_assembly["info"][key] = new_assembly["info"][key].map({0: np.nan})
    new_assembly["info"]["derived_from_split"] = new_assembly["info"]["derived_from_split"].fillna(False)
    print(ctime(), "Finished transposing")
    return new_assembly
