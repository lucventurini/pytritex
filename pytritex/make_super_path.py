import pandas as pd
import numpy as np
# import networkit as nk
# import networkx as nx
import graph_tool.topology
import graph_tool as gt
from .kopt2 import kopt2


def insert_node(df: pd.DataFrame, el: pd.DataFrame):

    def checker(df, el):
        inner_bait = el[(~el["cluster1"].isin(df["cluster"])) & (el["cluster2"].isin(df["cluster"]))]["cluster1"].unique()
        bait1 = el["cluster2"].isin(df["cluster"])
        bait2 = el["cluster1"].isin(inner_bait)

        temp = pd.merge(el.loc[bait1 & bait2].set_index("cluster2"),
                 df[["cluster", "bin"]].rename(columns={"cluster": "cluster2", "bin": "bin"}).set_index("cluster2"),
                 left_index=True, right_index=True, how="outer")
        y = temp[~temp["cluster1"].isna()].sort_values(["cluster1", "bin"]).agg(
            {"cluster1": }
        )
        y = y.sort_values(["cluster1", "bin"])
        idx = y.loc[y["dist"] == 1, :].index
        return y, idx

    y, idx = checker(df, el)
    el_col_map = {"cluster1": "path_node1", "cluster2": "path_node2", "weight": "old_path"}
    while idx.shape[0] > 0:
        z = pd.merge(el.rename(columns=el_col_map).set_index(["path_node1", "path_node2"]),
                     pd.DataFrame(
                         {
                             "cluster": y.loc[idx, "cluster1"],
                             "path_node1": y.loc[idx, "cluster2"],
                             "path_node2": y.loc[idx + 1, "cluster2"],
                             "weight1": y.loc[idx, "weight"],
                             "weight2": y.loc[idx + 1, "weight"]
                }
                     ).set_index(["path_node1", "path_node2"]),
                     left_index=True, right_index=True)
        z = z.assign(diff=z["weight1"] + z["weight1"] - z["old_path"]).sort_values("diff").head(1)
        m = z["bin"]
        df = pd.DataFrame(
            pd.concat([
                df.loc[range(0, m), :][:],
                pd.DataFrame({"cluster": z["cluster"], "bin": m}),
                df.loc[range(m + 1, df.shape[0] + 1)].assign(bin=df["bin"] + 1)[:]
            ])
        )
        y, idx = checker(df, el)
    return df


def make_super_path(super_info, idx=None, start=None, end=None, maxiter=100, verbose=True):

    submem = super_info["mem"][
        super_info["mem"].eval("super == {}".format(idx) if idx is not None else "super.isna()")]
    el = super_info["edges"][super_info["edges"].eval(
        "super == {}".format(idx) if idx is not None else "super.isna()"
    )][["cluster1", "cluster2", "weight"]]
    graph_view = gt.GraphView(super_info["graph"],
                              # Get the indices in the graph vertex list that are in the submem["cluster"] column
                              # TODO However, we might actually have to use the indices, not the "cluster" column
                              vfilt=np.in1d(super_info["graph"].get_vertices(), submem["cluster"])
                              )
    mst = gt.topology.min_spanning_tree(graph_view)
    m_graph_view = gt.GraphView(graph_view, efilt=mst)
    if start is None or end is None:
        # TODO this is wrong, will fail
        start, end = graph_tool.topology.pseudo_diameter(m_graph_view)[1]
    vertices = graph_tool.topology.shortest_path(m_graph_view, start, end)
    # Now back to the names. WE NEED TO KEEP THE ORDER.
    vertices = submem.set_index("index").loc[vertices]["cluster"]
    df = pd.DataFrame({"cluster": vertices, "bin": range(1, vertices.shape[0])})
    df = insert_node(df, el)
    df = kopt2(df, el, super_info["cluster"])

    df <- node_relocation(df, el, maxiter=maxiter, verbose=verbose)

    # Now we have to get the names of these vertices
        # Now: find the longest geodesic, ie the first (in a breadth-first search) path with the length equal to the diameter
        # Then go back to the names

#  df<-insert_node(df, el)
#  # Traveling salesman heuristics
#  df<-kopt2(df, el)
#  df<-node_relocation(df, el, maxiter=maxiter, verbose=verbose)
#
#  data.frame(df)->df
#  data.frame(el)->el
#  data.frame(cluster=df$cluster, rank = 0)->ranks
#  r=0
#
#  while(length(n<-unique(subset(el, !cluster1 %in% df$cluster & cluster2 %in% df$cluster)$cluster1)) > 0) {
#   r = r+1
#   subset(el, cluster2 %in% df$cluster & cluster1 %in% n)->tmp
#   tmp[!duplicated(tmp$cluster1),]->tmp
#   rbind(ranks, data.frame(cluster=tmp$cluster1, rank=r))->ranks
#   merge(tmp, df[c("cluster", "bin")], by.x="cluster2", by.y="cluster")->x
#   rbind(df, data.frame(cluster=x$cluster1, bin=x$bin))->df
#  }
#
#  merge(df, submem)->df
#  df$bin<-as.integer(df$bin)
#  flip<-with(df, suppressWarnings(cor(bin, cM))) < 0
#  if((!is.na(flip)) & flip) {
#   with(df, max(bin) - bin + 1)->df$bin
#  }
#  ranks$rank<-as.numeric(ranks$rank)
#  merge(df, ranks)->df
#  df$backbone <- df$cluster %in% dia
#
#  data.table(df)[, .(cluster, bin, rank, backbone)]
# }