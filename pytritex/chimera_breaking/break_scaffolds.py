import pandas as pd
import numpy as np
from .transposition import _transpose_cssaln, _transpose_molecules, _transpose_fpairs
from .transposition import _transpose_hic_cov, _transpose_molecule_cov
from .calculate_broken_scaffolds import calculate_broken_scaffolds
from ..anchoring import anchor_scaffolds
from time import ctime
import dask.dataframe as dd
import os


def break_scaffolds(breaks, save_dir, assembly, slop, cores=1, species="wheat") -> dict:

    fai = dd.read_parquet(assembly["fai"])
    new_assembly = calculate_broken_scaffolds(breaks, fai=fai, slop=slop)
    for key in ["binsize", "innerDist", "minNbin"]:
        new_assembly[key] = assembly[key]
    assert isinstance(new_assembly["fai"], dd.DataFrame)
    # Now extract the correct rows and columns from the FAI:
    # ie excluding the rows of the original scaffolds that have since been split.
    _left1 = new_assembly["fai"].query("derived_from_split == True")
    _left1_indices = _left1["orig_scaffold_index"].compute()
    # Now, let us get those that are NOT split AND do not have split children
    _left2 = new_assembly["fai"].query("derived_from_split == False")
    _left2 = new_assembly["fai"].query("scaffold_index not in @_left1_indices",
                                       local_dict=locals())
    # Put the unbroken scaffolds first
    new_fai = dd.concat([_left2, _left1])
    new_shape = new_fai["orig_scaffold_index"].unique().shape[0].compute()
    old_shape = fai["orig_scaffold_index"].unique().shape[0].compute()
    assert new_shape == old_shape, (new_shape, old_shape)
    base = os.path.join(save_dir, "joblib", "pytritex", "chimera_breaking")
    fname = os.path.join(base, "fai")
    dd.to_parquet(new_fai, fname, compression="gzip", engine="pyarrow")

    # Remove broken scaffolds from CSS alignment, put correct ones
    print(ctime(), "Transposing the CSS alignment")
    new_assembly["cssaln"] = _transpose_cssaln(assembly["cssaln"], new_fai)
    fname = os.path.join(base, "cssaln")
    dd.to_parquet(new_assembly["cssaln"], fname, compression="gzip", engine="pyarrow")

    # Do the same with HiC pairs
    print(ctime(), "Transposing the HiC alignments")
    fpairs = assembly.get("fpairs", None)
    if fpairs is not None:
        fpairs = dd.read_parquet(fpairs)
    new_assembly["fpairs"] = _transpose_fpairs(fpairs, new_fai)
    fname = os.path.join(base, "fpairs")
    dd.to_parquet(new_assembly["fpairs"], fname, compression="gzip", engine="pyarrow")

    print(ctime(), "Transposing the 10X alignments")
    molecules = assembly.get("molecules", None)
    if molecules is not None:
        molecules = dd.read_parquet(molecules)
    new_assembly["molecules"] = _transpose_molecules(molecules, new_fai)
    fname = os.path.join(base, "molecules")
    dd.to_parquet(new_assembly["molecules"], fname, compression="gzip", engine="pyarrow")

    print(ctime(), "Anchoring the transposed data")
    new_assembly = anchor_scaffolds(new_assembly, popseq=assembly["popseq"], species=species)
    new_assembly["info"] = new_assembly["info"].drop("mr_10x", axis=1, errors="ignore")
    new_assembly["info"] = new_assembly["info"].drop("mr", axis=1, errors="ignore")
    new_assembly["info"] = new_assembly["info"].drop("mri", axis=1, errors="ignore")
    print(ctime(), "Transposing the HiC coverage data")
    new_assembly = _transpose_hic_cov(new_assembly, assembly["info"], new_fai,
                                      assembly["cov"], assembly["fpairs"], cores=cores)
    print(ctime(), "Transposing the 10X coverage data")
    new_assembly = _transpose_molecule_cov(new_assembly, new_fai, assembly, cores=cores)
    # Change values that are 0 to nan
    for key in ["popseq_alphachr2", "popseq_alphachr", "sorted_alphachr", "sorted_alphachr2"]:
        new_assembly["info"][key] = new_assembly["info"][key].map({0: np.nan})
    new_assembly["info"]["derived_from_split"] = new_assembly["info"]["derived_from_split"].fillna(False)
    print(ctime(), "Finished transposing")
    return new_assembly
