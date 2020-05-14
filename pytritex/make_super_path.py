import pandas as pd
import numpy as np
import networkit as nk
import networkx as nx
from .kopt2 import kopt2
import scipy.sparse.csgraph
import scipy.sparse
from .insert_node import insert_node
from .node_relocation import node_relocation


def make_super_path(super_object, idx=None, start=None, end=None, maxiter=100, verbose=True, ncores=1):
    """Order scaffolds using Hi-C links for one chromosome."""

    # Get backbone from minimum spanning tree (MST)
    sub_membership = super_object["membership"].loc[lambda df: df["super"] == idx]
    edge_list = super_object["edges"].loc[lambda df: df["super"] == "idx",
                                          ["cluster1", "cluster2", "cidx1", "cidx2", "weight"]]
    cidx = pd.DataFrame(
        {"cidx": np.concatenate([edge_list["cidx1"], edge_list["cidx2"]]),
         "cluster": np.concatenate([edge_list["cluster1"], edge_list["cluster2"]])}
    ).drop_duplicates().set_index("cidx")
    _el_list = edge_list[["cidx1", "cidx2", "weight"]].values
    assert isinstance(_el_list, np.array)
    weights, coords = _el_list[:, 2], (_el_list[:, 0], _el_list[:, 1])
    shape = [(_el_list.max(axis=0)[:2] + 1).max()] * 2
    matrix = scipy.sparse.csc_matrix((weights, coords), shape=shape, dtype=_el_list.dtype)

    # the sparse matrix from nk will already contain the weights
    # matrix = nk.graphtools.subgraphFromNodes(graph, pd.concat([edge_list["cidx1"], edge_list["cidx2"]]).unique())
    mst = scipy.sparse.csgraph.minimum_spanning_tree(matrix)
    # Now back to NK
    nk.setNumberOfThreads(ncores)
    mgraph = nk.Graph(mst.shape[0], weighted=True, directed=False)
    [mgraph.addEdge(x, y, mst[x, y]) for x, y in zip(mst.nonzero()[0], mst.nonzero()[1])]
    # Now that we have the new minimum-spanning-tree graph we can get the diameter or the distance
    if start is None or end is None:
        diameter = nk.distance.Diameter(mgraph).run().getDiameter()[0]
        best = np.where(np.array(nk.distance.APSP(mgraph).run().getDistances()) == diameter)
        start, end = best[0][0], best[1][0]
    path = nk.distance.BFS(mgraph, start, target=end).run().getPath(end)
    # Bins here are ordered from 1 to X
    df = pd.DataFrame().assign(cluster=cidx.loc[path, "cluster"],
                               bin=range(len(path)))
    df = insert_node(df, edge_list)
    # # Traveling salesman heuristics
    df = kopt2(df, edge_list)
    df = node_relocation(df, edge_list, maxiter=maxiter, verbose=verbose)
    ranks = pd.DataFrame().assign(cluster=df["cluster"], rank=0)
    r = 0

    n = edge_list.loc[(~edge_list["cluster1"].isin(df["cluster"])) &
                      (edge_list["cluster2"].isin(df["cluster"])), "cluster1"].unique()

    while n.shape[0] > 0:
        r += 1
        tmp = edge_list.loc[(edge_list["cluster2"].isin(df["cluster"])) &
                            (edge_list["cluster1"].isin(n))].loc[lambda _tmp: ~_tmp["cluster1"].duplicated()]
        # rbind(ranks, data.frame(cluster=tmp$cluster1, rank = r))->ranks
        ranks = pd.concat([ranks, pd.DataFrame().assign(cluster=tmp["cluster1"], rank=r)])
        #   merge(tmp, df[c("cluster", "bin")], by.x="cluster2", by.y="cluster")->x
        # Do a merge and drop superfluous column
        x = pd.merge(tmp, df[["cluster", "bin"]], left_on="cluster2", right_on="cluster").drop("cluster", axis=1)
        df = pd.concat([df, pd.DataFrame().assign(cluster=x["cluster1"], bin=x["bin"])])
        n = edge_list.loc[(~edge_list["cluster1"].isin(df["cluster"])) &
                          (edge_list["cluster2"].isin(df["cluster"])), "cluster1"].unique()
        continue

    df = pd.merge(df, sub_membership)
    df.loc[:, "bin"] = df.loc[:, "bin"].astype(np.int)
    # Now calculate the correlation between bins and cM. If the correlation is negative
    # we have to flip the order!
    flip = df["bin"].corr(df["cM"], method="pearson")
    if not np.isnan(flip) and flip < 0:
        # Flip the order
        df.loc[:, "bin"] = df["bin"].max() - df["bin"] + 1
    ranks.loc[:, "rank"] = pd.to_numeric(ranks["rank"])
    df = pd.merge(df, ranks)
    df.loc[:, "backbone"] = df["cluster"].isin(path)
    return df[["cluster", "bin", "rank", "backbone"]]
