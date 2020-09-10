import dask.dataframe as dd
import pandas as pd
import numpy as np
from ..graph_utils.make_super import make_super


def make_hic_info(cluster_info, super_global, chroms):

    """"""

    # make_hic_info<-function(cluster_info, super_global, chrs){
    #  s<-super_global$super_info
    #  s[!duplicated(s$chr),]->s
    #  s[chr %in% chrs]->s
    #
    #  super_global$membership[, .(cluster, super, bin, rank, backbone)]->tmp
    #  tmp[super %in% s$super]->tmp
    #  tmp[, super := NULL]
    #  setnames(tmp, c("cluster", "hic_bin", "hic_rank", "hic_backbone"))
    #  tmp[cluster_info, on="cluster"]->cluster_info
    #  cluster_info[order(chr, hic_bin, hic_rank, cluster)]
    # }

    sinfo = super_global["super_info"].drop_duplicates(subset=["chr"]).query(
        "chr in @chroms", local_dict={"chroms": chroms})





def make_hic_map(hic_info, links, client, ncores=1, maxiter=100, known_ends=True):
    # copy(links)->hl
    #  copy(hic_info)->info
    #
    #  setnames(hl, c("scaffold1", "scaffold2"), c("cluster1", "cluster2"))
    #  setnames(info, "scaffold", "cluster")
    #  setkey(info, "cluster")
    #  chrs <- info[!is.na(chr), unique(chr)]
    # super_global<-make_super(hl, cluster_info=info, cores=ncores, maxiter=maxiter,
    #     # 			   known_ends=known_ends, path_max=length(chrs))
    #  make_hic_info(info, super_global, chrs=chrs)->res
    #  res[order(chr, hic_bin)][, .(scaffold=cluster, chr, cM, hic_bin, hic_backbone, hic_rank)][!is.na(hic_bin)]

    hl = links.rename(columns={"scaffold_index1": "cluster1", "scaffold_index2": "cluster2"})
    if hic_info.index.name != "scaffold_index":
        info = hic_info.rename(columns={"scaffold_index": "cluster"}).set_index("cluster")
    else:
        info = hic_info
        info.index = info.index.rename("cluster")
    prev_mem = pd.DataFrame()
    chroms = info.query("chr == chr")["chr"].unique().values.compute()
    super_global = make_super(hl, cluster_info=info, cores=ncores, maxiter=maxiter,
                              known_ends=known_ends, path_max=chroms.shape[0], client=client,
                              previous_membership=prev_mem)
    result = make_hic_info(info, super_global, chroms=chroms)
    result = result.sort_values(["chr", "hic_bin"])[["chr", "cM", "hic_bin", "hic_backbone", "hic_rank"]]
    result.index = result.index.rename(columns={"cluster": "scaffold"})
    result = result.query("hic_bin == hic_bin").perist()
    return result


def hic_map(fpairs, info, popseq_dist, cores=1, min_nfrag_scaffold=50):

    hic_info = info.copy()
    hic_info["excluded"] = hic_info.eval("nfrag < @min_nfrag_scaffold",
                                         local_dict={"min_nfrag_scaffold": min_nfrag_scaffold})
    hic_info = hic_info.persist()
    base = hic_info[["chr", "cM"]]

    assert base.index.name == "scaffold_index"
    right1 = base.rename(columns={"chr": "chr1", "cM": "cM1"})
    right1.index = right1.index.rename("scaffold_index1")
    fpairs = fpairs.query("scaffold_index1 != scaffold_index2")
    fpairs = dd.merge(fpairs.set_index("scaffold_index1"), right1, on="scaffold_index1", how="left")
    fpairs = dd.merge(fpairs.set_index("scaffold_index1"), right1, on="scaffold_index1", how="left")
    right2 = base.rename(columns={"chr": "chr2", "cM": "cM2"})
    right2.index = right2.index.rename("scaffold_index2")
    fpairs = dd.merge(fpairs.reset_index(drop=False).set_index("scaffold_index2"),
                      right2, on="scaffold_index2", how="left").reset_index(drop=False)
    fpairs = fpairs.query("chr1 == chr2")
    fpairs["cMdist"] = (fpairs["cM1"] - fpairs["cM2"]).abs()
    fpairs = fpairs.query("cM1 != cM1 | cM2 != cM2 | cMdist <= @popseq_dist", local_dict={"popseq_dist": popseq_dist})
    fpairs = fpairs.persist()
    nlinks = fpairs.groupby(["scaffold_index1", "scaffold_index2"]).size().to_frame("nlinks")
    fpairs = fpairs.reset_index(drop=False).merge(nlinks, on=["scaffold_index1", "scaffold_index2"])
    fpairs = fpairs.set_index("hic_index")
    fpairs["weight"] = fpairs.eval("log(nlinks) / log(10)")
