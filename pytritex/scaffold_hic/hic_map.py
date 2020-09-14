import dask.dataframe as dd
import pandas as pd
import numpy as np
import os
from typing import Union


def hic_map(info: Union[dd.DataFrame, pd.DataFrame],
            assembly: dict, frags, species, ncores=1, min_nfrag_scaffold=50, max_cM_dist = 20,
            binsize=5e5, min_nfrag_bin=30, gap_size=100, maxiter=100, orient=True, agp_only=False,
            map=None, known_ends=True, orient_old=False, min_binsize=1e5, min_nbin=5):

    # copy(info)->hic_info
    #   hic_info[, excluded := nfrag < min_nfrag_scaffold]
    #
    #   assembly$fpairs[scaffold1 != scaffold2, .(nlinks=.N), key=.(scaffold1, scaffold2)]->hl
    #   hic_info[, .(scaffold1=scaffold, chr1=chr, cM1=cM)][hl, nomatch=0, on="scaffold1"]->hl
    #   hic_info[, .(scaffold2=scaffold, chr2=chr, cM2=cM)][hl, nomatch=0, on="scaffold2"]->hl
    #   hl[chr1 == chr2]->hl
    #   hl<-hl[abs(cM1-cM2) <= max_cM_dist | is.na(cM1) | is.na(cM2)]
    #   hl[, weight:=-log10(nlinks)]
    #
    #   cat("Scaffold map construction started.\n")
    #   make_hic_map(hic_info=hic_info, links=hl, ncores=ncores, known_ends=known_ends)->hic_map
    #   cat("Scaffold map construction finished.\n")

    hic_info = info.copy()
    hic_info["excluded"] = hic_info.eval("nfrag < @mnf", local_dict={"mnf": min_nfrag_scaffold})
    hl = assembly["fpairs"].eval("scaffold_index1 != scaffold_index2").persist()
    hl = hl.merge(hl.groupby(["scaffold_index1", "scaffold_index2"]).size().to_frame("nlinks").compute(),
                  on=["scaffold_index1", "scaffold_index2"]).persist()
    left = hic_info[["chr", "cM"]]
    left1 = left.rename(columns={"chr": "chr1", "cM": "cM1"})
    left1.index = left1.index.rename("scaffold_index1")
    left2 = hic_info.rename(columns={"chr": "chr2", "cM": "cM2"})
    left2.index = left1.index.rename("scaffold_index2")
    hl = hl.merge(left1, on="scaffold_index1", how="inner").merge(left2, on="scaffold_index2", how="inner").query(
        "chr1 == chr2")
    hl = hl[(hl["cM1"].isna()) | (hl["cM2"].isna()) | ((hl["cM1"] - hl["cM2"]).abs() <= max_cM_dist)]
    # TODO: why is the log minused?
    hl["weight"] = -np.log10(hl["nlinks"])
    hic_map = make_hic_map(hic_info=hic_info, links=hl, ncores=ncores, known_ends=known_ends)