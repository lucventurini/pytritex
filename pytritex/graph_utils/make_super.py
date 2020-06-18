import pandas as pd
import numpy as np
# import graph_tool
# import graph_tool as gt
# import graph_tool.clustering
import time
import networkit as nk
import dask.dataframe as dd
from pytritex.graph_utils.make_super_path import make_super_path


# make_super<-function(hl, cluster_info, prefix="super", cores=1, paths=T, path_max=0, known_ends=F,
# 		     maxiter=100, verbose=T){
#
#  hl[cluster1 %in% cluster_info[excluded == F]$cluster & cluster2 %in% cluster_info[excluded == F]$cluster]->hl
#
#  hl[cluster1 < cluster2]->e
#  graph.edgelist(as.matrix(e[, .(cluster1, cluster2)]), directed=F)->g
#  E(g)$weight<-e$weight
#
#  data.table(cluster=V(g)$name, super=paste(prefix, sep="_", clusters(g)$membership))->mem
#  cluster_info[mem, on="cluster"]->mem
#
#  mem[, .(super_size=.N, length=.N, chr=unique(na.omit(chr))[1], cM=mean(na.omit(cM))), keyby=super]->info
#  mem[, .(cluster1=cluster, super)][hl, on="cluster1"]->e
#
#  list(super_info=info, membership=mem, graph=g, edges=e)->s
#
#  if(paths){
#   if(path_max > 0){
#    idx<-head(s$super_info[order(-length)], n=path_max)$super
#   } else {
#    idx<-s$super_info$super
#   }
#   rbindlist(mclapply(mc.cores=cores, idx, function(i) {
#    start <- end <- NULL
#    # Take terminal nodes from genetic map
#    if(known_ends){
#     s$mem[super == i & !is.na(cM)][order(cM)]$cluster->x
#     start=x[1]
#     end=tail(x,1)
#    }
#    make_super_path(s, idx=i, start=start, end=end, maxiter=maxiter, verbose=verbose)->x
#    if(verbose){
#     cat(paste0("Chromosome ", head(s$mem[super == i]$chr, n=1), " finished.\n"))
#    }
#    x
#   }))[s$membership, on="cluster"]->s$membership
#  }
#
#  s
# }


def _concatenator(edges, membership, known_ends=False, maxiter=100, verbose=False):
    start = end = None
    if known_ends:
        x = membership.query("cM == cM").sort_values("cM").index
        start, end = x.head(1), x.tail(1)
    final = make_super_path(edges, membership,
                            start=start, end=end,
                            maxiter=maxiter, verbose=verbose)
    return final


def make_super(hl: dd.DataFrame, cluster_info: dd.DataFrame,
               cores=1, paths=True, path_max=0, known_ends=False,
               maxiter=100, verbose=True):

    #  hl[cluster1 %in% cluster_info[excluded == F]$cluster & cluster2 %in% cluster_info[excluded == F]$cluster]->hl
    hl = hl.copy()
    non_excluded = cluster_info.loc[cluster_info["excluded"] == False].index.values.compute()
    bait1 = hl["cluster1"].isin(non_excluded)
    bait2 = hl["cluster2"].isin(non_excluded)
    hl = hl.loc[bait1 & bait2, :]
    # hl = hl.loc[(hl["cluster1"].isin(
    #     cluster_info.loc[~cluster_info["excluded"]].index.values))
    #     & (hl["cluster2"].isin(cluster_info.loc[~cluster_info["excluded"]].index.values))]
    edge_list = hl.query("cluster1 < cluster2")[["cluster1", "cluster2", "weight"]].compute()
    # Unique cluster IDs
    cidx = np.unique(edge_list[["cluster1", "cluster2"]].values.flatten())
    cidx = np.vstack([cidx, np.arange(cidx.shape[0])])
    cidx = pd.DataFrame().assign(cluster=cidx[0, :], cidx=cidx[1, :]).set_index("cluster")

    # Now merge into the edge table
    cidx1 = cidx[:]
    cidx1.index.names = ["cluster1"]
    cidx1.columns = ["cidx1"]
    cidx2 = cidx[:]
    cidx2.index.names = ["cluster2"]
    cidx2.columns = ["cidx2"]

    edge_list = edge_list.merge(cidx1, how="inner", on="cluster1").merge(cidx2, how="inner", on="cluster2")
    edge_list = edge_list.reset_index(drop=False)
    # Now we are ready to create the graph using the indices.
    #  hl[cluster1 < cluster2]->e
    #  graph.edgelist(as.matrix(e[, .(cluster1, cluster2)]), directed=F)->g
    #  E(g)$weight<-e$weight
    nk.setNumberOfThreads(cores)
    graph = nk.Graph(n=cidx.shape[0], weighted=True, directed=False)
    edge_list.apply(lambda row: graph.addEdge(row["cidx1"], row["cidx2"], row["weight"]), axis=1)
    components = nk.components.ConnectedComponents(graph)
    # Get the connected components
    components.run()
    cidx_list, ssuper = [], []
    for idx, comp in enumerate(components.getComponents()):
        cidx_list.extend(comp)
        ssuper.extend([idx] * len(comp))
        # ssuper.extend(["{prefix}_{idx}".format(prefix=prefix, idx=idx)] * len(comp))

    membership = pd.DataFrame().assign(cidx=cidx_list, super=ssuper).set_index("cidx")
    membership = membership.merge(cidx.reset_index(drop=False), on=["cidx"], how="outer")
    membership = cluster_info.merge(membership, on="cluster", how="right")
    # info <- mem[, .(super_size=.N, length=.N, chr=unique(na.omit(chr))[1], cM=mean(na.omit(cM))),
    # keyby=super]
    # chr=("chr", lambda s: np.nan if s.dropna().shape[0] == 0 else s.dropna().unique()[0])

    grouped = membership.groupby("super")
    info = grouped.agg(
        {"super": "size",
         "cM": "mean"}
    )
    # print(info.head(npartitions=-1, n=5))
    info["chr"] = grouped["chr"].unique().apply(lambda s: [_ for _ in s if not np.isnan(_)][0],
                                                meta=np.float)
    print(info.head(npartitions=-1, n=5))
    assert "cluster" in membership, membership.columns
    edge_list = membership.rename(columns={"cluster": "cluster1"})[["cluster1", "super"]].merge(
        hl, on="cluster1", how="right")

    super_object = {"super_info": info, "membership": membership, "graph": graph, "edges": edge_list}

    if paths is True:
        cms = super_object["membership"].loc[
              :, ["cluster", "super", "cM"]].drop_duplicates().set_index("cluster")

        if path_max > 0:
            idx = super_object["super_info"].sort_values("length", ascending=False).head(path_max)["super"].unique()
        else:
            idx = super_object["super_info"]["super"].unique()
        grouped_edges = super_object["edges"].groupby("super").groups
        grouped_membership = cms.groupby("super").groups
        order = super_object["membership"].query("super in @idx")[[
            "super", "chr"]].sort_values(["chr"])
        print(time.ctime(), "Starting super, with", idx.shape[0], "positions to analyse.")
        results = []
        for row in order.itertuples(name=None):
            index, ssuper, popseq_chr = row
            my_edges = super_object["edges"].loc[grouped_edges[ssuper]]
            my_membership = cms.loc[grouped_membership[ssuper]]
            results.append(_concatenator(my_edges, my_membership, known_ends, maxiter, verbose))

        results = np.vstack(results)
        results = pd.DataFrame().assign(
            cluster=results[:, 0], bin=results[:, 1], rank=results[:, 2], backone=results[:, 3])
        super_object["membership"] = results.merge(super_object["membership"], on="cluster")
        assert "cluster" in super_object["membership"].columns

    print(time.ctime(), "Finished make_super run")
    # import sys
    # sys.exit(1)
    return super_object
