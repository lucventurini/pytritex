import pandas as pd
# import graph_tool
# import graph_tool as gt
# import graph_tool.clustering
import time
import networkit as nk
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


def _concatenator(super_object, ssuper, known_ends=False, maxiter=100, verbose=False, prev_chrom=None):
    start = end = None
    if known_ends:
        x = super_object["membership"].query("(super == @ssuper) & (cM == cM)").sort_values("cM")["cluster"]
        start, end = x.head(1), x.tail(1)
    final = make_super_path(super_object,
                            idx=ssuper, start=start, end=end,
                            maxiter=maxiter, verbose=verbose)
    return final


def make_super(hl, cluster_info,
               prefix="super", cores=1, paths=True, path_max=0, known_ends=False,
               maxiter=100, verbose=True):

    #  hl[cluster1 %in% cluster_info[excluded == F]$cluster & cluster2 %in% cluster_info[excluded == F]$cluster]->hl
    hl = hl.copy()
    hl = hl.loc[(hl["cluster1"].isin(cluster_info.loc[~cluster_info["excluded"], "cluster"]))
                & (hl["cluster1"].isin(cluster_info.loc[~cluster_info["excluded"], "cluster"]))]
    edge_list = hl.loc[hl.eval("cluster1 < cluster2"), ["cluster1", "cluster2", "weight"]]
    # Unique cluster IDs
    cidx = pd.Series(
        pd.concat([edge_list.cluster1, edge_list.cluster2]).unique()
        ).rename("cluster").reset_index(drop=False)[["cluster", "index"]].rename(columns={"index": "cidx"})
    # Now merge into the edge table
    edge_list = edge_list.merge(
        cidx.rename(columns={"cluster": "cluster1", "cidx": "cidx1"}), how="inner", on="cluster1").merge(
        cidx.rename(columns={"cluster": "cluster2", "cidx": "cidx2"}), how="inner", on="cluster2")
    # Now we are ready to create the graph using the indices.
    #  hl[cluster1 < cluster2]->e
    #  graph.edgelist(as.matrix(e[, .(cluster1, cluster2)]), directed=F)->g
    #  E(g)$weight<-e$weight
    graph = nk.Graph(n=cidx.shape[0], weighted=True, directed=False)
    edge_list.apply(lambda row: graph.addEdge(row["cidx1"], row["cidx2"], row["weight"]), axis=1)
    components = nk.components.ConnectedComponents(graph)
    # Get the connected components
    components.run()
    membership = pd.concat([pd.DataFrame().assign(
        cidx=comp,
        super="{prefix}_{idx}".format(prefix=prefix, idx=idx))
        for idx, comp in enumerate(components.getComponents())])
    membership = membership.merge(cidx, on=["cidx"])
    membership = cluster_info.merge(membership, on="cluster", how="right")
    # info <- mem[, .(super_size=.N, length=.N, chr=unique(na.omit(chr))[1], cM=mean(na.omit(cM))),
    # keyby=super]
    info = membership.groupby("super").agg(
        super_size=("super", "size"), length=("super", "size"), cM=("cM", "mean")).reset_index(drop=False)
    assert "cluster" in membership, membership.columns
    edge_list = membership.rename(columns={"cluster": "cluster1"})[["cluster1", "super"]].merge(
        hl, on="cluster1", how="right")

    super_object = {"super_info": info, "membership": membership, "graph": graph, "edges": edge_list}

    if paths is True:
        if path_max > 0:
            idx = super_object["super_info"].sort_values("length", ascending=False).head(path_max)["super"].unique()
        else:
            idx = super_object["super_info"]["super"].unique()

        order = super_object["membership"].query("super in @idx")[[
            "super", "popseq_chr"]].sort_values(["popseq_chr"])
        print(time.ctime(), "Starting super, with", idx.shape[0], "positions to analyse.")
        results = []
        prev_chrom = None
        print(order["popseq_chr"].value_counts())
        for row in order.itertuples(name=None):
            index, ssuper, popseq_chr = row
            popseq_chr = int(popseq_chr)
            check = (popseq_chr == prev_chrom)
            assert check in (False, True)
            # print(type(prev_chrom), prev_chrom, type(popseq_chr), popseq_chr, popseq_chr == prev_chrom)
            if check is False:
                if prev_chrom is not None:
                    print(time.ctime(), "Finished chromosome", prev_chrom)
                print(time.ctime(), "Starting chromosome", popseq_chr)
            prev_chrom = int(popseq_chr)
            results.append(_concatenator(super_object, ssuper, known_ends, maxiter, verbose))

        super_object["membership"] = pd.concat(results).merge(
            super_object["membership"], on="cluster")

    print(time.ctime(), "Finished make_super run")
    import sys
    sys.exit(1)
    return super_object
