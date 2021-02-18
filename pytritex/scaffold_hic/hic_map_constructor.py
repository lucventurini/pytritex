import dask.dataframe as dd
import pandas as pd
from typing import Union
from ..graph_utils.make_super import make_super, add_statistics
from ..graph_utils.tip_removal.tip_remover import remove_tips
from dask.distributed import Client
import os
import dask.array as da
import numpy as np


def make_hic_map(hic_info: Union[pd.DataFrame, dd.DataFrame],
                 links: Union[pd.DataFrame, dd.DataFrame],
                 client: Client,
                 save_dir=None,
                 verbose=False,
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
    cluster_info = hic_info.copy()
    cluster_info.index = cluster_info.index.rename("cluster")
    chrs = cluster_info.query("chr == chr")["chr"].unique().compute()
    # Execute make_super until no branches of length > 1 are present

    # counter = 1
    # out, excluded, run = _iterator(counter=counter, membership=None, excluded=excluded)
    # membership = out["membership"]
    run = True
    excluded_scaffolds = cluster_info["excluded"].compute()
    # excluded = excluded_scaffolds.query(
    
    while run is True:
        cluster_info["excluded"] = excluded_scaffolds
        super_object = make_super(hl, cluster_info=cluster_info, client=client,
                      cores=ncores, maxiter=maxiter, known_ends=known_ends,
                      path_max=chrs.shape[0], previous_membership=pd.DataFrame())
        membership = super_object["membership"]
        # dd_membership = dd.read_parquet(membership, infer_divisions=True)
        # Now we have to exclude from consideration those scaffolds
        # that are the backbone of super-scaffolds where there is at least a scaffold with
        # rank > 1.
        # Ie: remove the backbones so that the branches can be reassigned to somewhere else.
        a = membership.merge(
            membership.loc[membership["rank"] > 1, ["super", "bin"]].reset_index(drop=False).drop_duplicates(),
            how="inner", on=["super", "bin"])
        assert "cluster" in a.columns
        add = a.loc[a["rank"] == 0, :]
        if add.shape[0].compute() == 0:
            run = False
        else:
            # Drop all the links between the backbone of the "fuzzy" scaffolds and the spikes.
            _to_exclude = set(add["cluster"].values.compute().tolist())
            excluded_scaffolds.loc[_to_exclude] = True
            
    excluded = set(excluded_scaffolds[excluded_scaffolds == True].index.values)
    super_object["membership"].index = super_object["membership"].index.rename("scaffold_index")
    # membership = super_object["membership"]
    super_object["membership"], super_object["super_info"] = add_statistics(super_object["membership"].reset_index(drop=False),
                                                                            client)

    membership, res, excluded = remove_tips(
        links=links, excluded=excluded,
        out=super_object, info=hic_info,
        client=client,
        save_dir=os.path.join(save_dir, "tip_removal"),
        ncores=ncores, verbose=verbose,
        min_dist=1e4)

    super_object = {"membership": membership, "super_info": res}

    # TODO This is probably the wrong key
    # super_info = super_object["super_info"].drop_duplicates(["chr"]).query("chr in @chrs", local_dict={"chrs": chrs})
    super_object["super_info"] = dd.read_parquet(super_object["super_info"], infer_divisions=True)
    super_info = super_object["super_info"]
    # assert super_object["super_info"].shape[0].compute() == chrs.shape[0]
    super_object["membership"] = dd.read_parquet(super_object["membership"], infer_divisions=True)
    # Get a temporary membership table
    #  super_global$membership[, .(cluster, super, bin, rank, backbone)]->tmp
    if super_object["membership"].index.name == "scaffold_index":
        tmp_memb = super_object["membership"][["super", "bin", "rank", "backbone"]]
    else:
        tmp_memb = super_object["membership"][
            ["scaffold_index", "super", "bin", "rank", "backbone"]].set_index("scaffold_index")
        
    #  tmp[super %in% s$super]->tmp
    #  tmp[, super := NULL]
    tmp_memb = tmp_memb.query(
        "super in @supers", local_dict={"supers": np.unique(super_info.index.values.compute())}).rename(
        columns={"bin": "hic_bin", "rank": "hic_rank", "backbone": "hic_backbone"})
    #  setnames(tmp, c("cluster", "hic_bin", "hic_rank", "hic_backbone"))
    #  tmp[cluster_info, on="cluster"]->cluster_info
    #  cluster_info[order(chr, hic_bin, hic_rank, cluster)]
    tmp_memb.columns = ["super", "hic_bin", "hic_rank", "hic_backbone"]
    chr_result = tmp_memb.merge(hic_info, on="scaffold_index").compute().reset_index(drop=False).sort_values(
        ["chr", "hic_bin", "hic_rank", "scaffold_index"])
    chr_result = chr_result.sort_values(["chr", "cM", "hic_bin"])[
        ["scaffold_index", "chr", "cM", "hic_bin", "hic_backbone", "hic_rank"]
    ].query("hic_bin == hic_bin").set_index("scaffold_index")
    assert isinstance(chr_result, pd.DataFrame), (type(chr_result), chr_result.head())
    chr_result = dd.from_pandas(chr_result, chunksize=int(1e6))
    if save_dir is not None:
        dd.to_parquet(chr_result, os.path.join(save_dir, "chr_result"))
        dd.to_parquet(super_object["super_info"], os.path.join(save_dir, "result"))
        dd.to_parquet(super_object["membership"], os.path.join(save_dir, "membership"))
        print(chr_result, os.path.join(save_dir, "chr_result"))

    return chr_result
