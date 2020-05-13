import pandas as pd
# import graph_tool
import multiprocessing as mp
import numpy as np
import itertools
import graph_tool as gt
import graph_tool.clustering
from .make_super_path import make_super_path


def _concatenator(s, i, group, known_ends=False, maxiter=100, verbose=False):
    start, end = None, None
    if known_ends is True:
        mem = s["mem"]
        x = mem[(mem["super"] == i) & (~mem["cM"].isna())].sort_values("cM")["cluster"]
        start, end = x.head(1), x.tail(1)
    x = make_super_path(s, idx=i, start=start, end=end, maxiter=maxiter, verbose=verbose)
    return x

def make_super(hl: pd.DataFrame,
               cluster_info, prefix="super", cores=1, paths=True, path_max=0,
               known_ends=False, maxiter=100, verbose=True):

    bait1 = hl["cluster1"].isin(cluster_info[cluster_info["excluded"] == False]["cluster"])
    bait2 = hl["cluster2"].isin(cluster_info[cluster_info["excluded"] == False]["cluster"])
    hl = hl.loc[bait1 & bait2,:][:]
    g = gt.Graph(directed=False)
    weight = g.new_edge_property("int")
    e = hl[hl["cluster1"] < hl["cluster2"]]  # Edges?
    # Create undirected graph
    idx = cluster_info[["cluster"]].reset_index()["cluster"]
    idx["index"] = idx.index
    e["index1"] = pd.merge(idx, e["cluster1"], left_on="cluster", right_on="cluster1")["index"]
    e["index2"] = pd.merge(idx, e["cluster2"], left_on="cluster", right_on="cluster2")["index"]
    g.add_edge_list(e[["index1", "index2", "weight"]].to_numpy(), eprops=[weight])

    # Now get the component for each vertex
    mem = pd.merge(
        pd.DataFrame({"cidx": g.get_vertices(),
                      "super": graph_tool.topology.label_components(g)[0].get_array()}),
        idx, left_on="cidx", right_on="index")[["cluster", "super"]]

    if prefix is not None:
        mem["super"] = mem["super"].astype(str).str.replace("^(.)", "{}_\\1".format(prefix), regex=True)
    if mem.shape[0] == 0:
        return None  # This is the error

    def unique_value(col):
        if col.shape[0] == 0 or np.where(col.isna())[0].shape[0] == col.shape[0]:
            return np.nan
        else:
            return np.unique(col)[0]

    mem = pd.merge(cluster_info, mem, left_on="cluster", right_on="cluster")
    info = mem.dropna(axis=0, subset=["chr"]).groupby("super").agg({
        "super": ["count"], "chr": [unique_value], "cM": ["mean"]
    })
    # Flatten the multiindex, then rename
    info.columns = [col for col in info.columns.values]
    info = info.rename(columns={("super", "count"): "super_size",
                                ("cM", "mean"): "cM",
                                ("chr", "unique_value"): "chr"})
    e = pd.merge(mem[["cluster", "super"]].rename(columns={"cluster": "cluster1"}), hl,
                 left_on="cluster1", right_on="cluster1")
    s = {"super_info": info, "membership": mem, "graph": g, "edges": e}

    if paths is True:
        if path_max > 0:
            idx = s["super_info"].sort_values("length", ascending=False).head(path_max)["super"]
        else:
            idx = s["super_info"]["super"]
        pool = mp.Pool(cores)

        pd.concat()

    return s
 # if(paths){
 #  rbindlist(mclapply(mc.cores=cores, idx, function(i) {
 #   start <- end <- NULL
 #   # Take terminal nodes from genetic map
 #   if(known_ends){
 #    s$mem[super == i & !is.na(cM)][order(cM)]$cluster->x
 #    start=x[1]
 #    end=tail(x,1)
 #   }
 #   make_super_path(s, idx=i, start=start, end=end, maxiter=maxiter, verbose=verbose)->x
 #   if(verbose){
 #    cat(paste0("Chromosome ", head(s$mem[super == i]$chr, n=1), " finished.\n"))
 #   }
 #   x
 #  }))[s$membership, on="cluster"]->s$membership
 # }

