import pandas as pd
import numpy as np
from .transposition import _transpose_cssaln, _transpose_molecules, _transpose_fpairs
from .transposition import _transpose_hic_cov, _transpose_molecule_cov
from .calculate_broken_scaffolds import calculate_broken_scaffolds
from ..anchoring import anchor_scaffolds


def break_scaffolds(breaks, assembly, slop, cores=1, species="wheat", use_memory_fs=False) -> dict:

    new_assembly = calculate_broken_scaffolds(breaks, assembly, slop)
    for key in ["binsize", "innerDist", "minNbin"]:
        new_assembly[key] = assembly[key]
    # Now extract the correct rows and columns from the FAI:
    # ie excluding the rows of the original scaffolds that have since been split.
    _left1 = new_assembly["fai"].loc[new_assembly["fai"].derived_from_split]
    # Now, let us get those that are NOT split AND do not have split children
    _left2 = new_assembly["fai"].loc[~new_assembly["fai"].derived_from_split & (
            ~new_assembly["fai"]["scaffold_index"].isin(_left1["orig_scaffold_index"]))]
    fai = pd.concat([_left1, _left2])
    new_shape = fai["orig_scaffold_index"].unique().shape[0]
    old_shape = assembly["fai"]["orig_scaffold_index"].unique().shape[0]
    assert new_shape == old_shape, (new_shape, old_shape)

    new_assembly["cssaln"] = _transpose_cssaln(assembly["cssaln"], fai)
    new_assembly["fpairs"] = _transpose_fpairs(assembly.get("fpairs", None), fai)
    new_assembly["molecules"] = _transpose_molecules(assembly.get("molecules", None), fai)
    new_assembly = anchor_scaffolds(new_assembly, popseq=assembly["popseq"], species=species)
    new_assembly["info"] = new_assembly["info"].drop("mr_10x", axis=1, errors="ignore")
    new_assembly["info"] = new_assembly["info"].drop("mr", axis=1, errors="ignore")
    new_assembly["info"] = new_assembly["info"].drop("mri", axis=1, errors="ignore")
    new_assembly = _transpose_hic_cov(new_assembly, assembly["info"], fai,
                                      assembly["cov"], assembly["fpairs"], cores=cores, use_memory_fs=use_memory_fs)
    new_assembly = _transpose_molecule_cov(new_assembly, fai, assembly, cores=cores, use_memory_fs=use_memory_fs)
    for key in ["popseq_alphachr2", "popseq_alphachr", "sorted_alphachr", "sorted_alphachr2"]:
        new_assembly["info"].loc[lambda df: df[key] == 0, key] = np.nan
    new_assembly["info"].loc[:, "derived_from_split"] = new_assembly["info"]["derived_from_split"].fillna(False)
    return new_assembly