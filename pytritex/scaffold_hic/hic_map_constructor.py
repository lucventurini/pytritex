import dask.dataframe as dd
import pandas as pd
from typing import Union
from ..graph_utils.make_super import make_super
from dask.distributed import Client


def make_hic_map(hic_info: Union[pd.DataFrame, dd.DataFrame],
                 links: Union[pd.DataFrame, dd.DataFrame],
                 client: Client,
                 ncores=1, maxiter=100, known_ends=True):
    # make_hic_map<-function(hic_info, links, ncores=1, maxiter=100, known_ends=T){
    #
    #  copy(links)->hl
    #  copy(hic_info)->info
    #
    #  setnames(hl, c("scaffold1", "scaffold2"), c("cluster1", "cluster2"))
    #  setnames(info, "scaffold", "cluster")
    #  setkey(info, "cluster")
    #  chrs <- info[!is.na(chr), unique(chr)]
    #
    #  make_hic_info(info,
    #   super_global<-make_super(hl, cluster_info=info, cores=ncores, maxiter=maxiter,
    # 			   known_ends=known_ends, path_max=length(chrs)), chrs=chrs)->res
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

    #  res[order(chr, hic_bin)][, .(scaffold=cluster, chr, cM, hic_bin, hic_backbone, hic_rank)][!is.na(hic_bin)]
    # }
    hl = links.copy().rename(columns={"scaffold_index1": "cluster1", "scaffold_index2": "cluster2"})
    info = hic_info.copy()
    info.index = info.index.rename("cluster")
    chrs = info.query("chr == chr")["chr"].unique().compute()
    sups = make_super(hl, cluster_info=info, client=client,
                      cores=ncores, maxiter=maxiter, known_ends=known_ends,
                      path_max=chrs.shape[0], previous_membership=pd.DataFrame())
    # TODO This is probably the wrong key
    super_info = sups["super_info"].drop_duplicates(["chr"]).query("chr in @chrs", local_dict={"chrs": chrs})
    if sups["membership"].index.name == "cluster":
        tmp_memb = sups["membership"][["super", "bin", "rank", "backbone"]]
    else:
        tmp_memb = sups["membership"][["cluster", "super", "bin", "rank", "backbone"]].set_index("cluster")

    tmp_memb = tmp_memb.query(
        "super in @supers", local_dict={"supers": super_info["super"].unique().compute()}).drop(["super"], axis=1)
    tmp_memb.columns = ["hic_bin", "hic_rank", "hic_backbone"]
    result = tmp_memb.merge(info, on ="cluster").compute().reset_index(drop=False).sort_values(
        ["chr", "hic_bin", "hic_rank", "cluster"])
    result = result.sort_values(["chr", "hic_bin"])[["chr", "cM", "hic_bin", "hic_backbone", "hic_rank"]]
    result = result.query("hic_bin == hic_bin")
    return result
