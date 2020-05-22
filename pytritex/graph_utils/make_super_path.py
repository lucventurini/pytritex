import pandas as pd
import numpy as np
import networkit as nk
from sys import float_info
import scipy.sparse.csgraph
import scipy.sparse
from .k_opt_tsp import tsp_2_opt


# This is necessary for networkit, as currently the distances can be "infinite"
# See https://github.com/networkit/networkit/issues/541
max_double = float_info.max


# TODO: remove all unnecessary pandas code from here. It only slows down the function.
def make_super_path(super_object, idx, start=None, end=None, maxiter=100, verbose=True, ncores=1):
    """Order scaffolds using Hi-C links for one chromosome."""

    # Get backbone from minimum spanning tree (MST)
    sub_membership = super_object["membership"].query("super == @idx")
    edge_list = super_object["edges"].query(
        "(super == @idx) & (cluster1 != cluster2)")[["cluster1", "cluster2", "weight"]]

    # Get the new index and make sure it does not use much memory
    cidx = pd.DataFrame().assign(
        cluster=np.unique(edge_list[["cluster1", "cluster2"]].values.flatten())
    ).reset_index(drop=False).rename(columns={"index": "cidx"})
    cidx.loc[:, "cidx"] = pd.to_numeric(cidx["cidx"], downcast="unsigned")

    # Assign the new index
    edge_list = cidx.rename(columns={"cidx": "cidx1", "cluster": "cluster1"}).merge(edge_list, on="cluster1")
    edge_list = cidx.rename(columns={"cidx": "cidx2", "cluster": "cluster2"}).merge(edge_list, on="cluster2")
    coords, weights = edge_list[["cidx1", "cidx2"]].values, edge_list["weight"].values
    shape = [int(coords.max() + 1)] * 2
    try:
        matrix = scipy.sparse.coo_matrix((weights, (coords[:, 0], coords[:, 1])), dtype=weights.dtype,
                                         shape=shape)
    except TypeError:
        raise TypeError((weights, coords))
    except ValueError:
        raise ValueError((weights, coords))

    # the sparse matrix from nk will already contain the weights
    # matrix = nk.graphtools.subgraphFromNodes(graph, pd.concat([edge_list["cidx1"], edge_list["cidx2"]]).unique())
    mst = scipy.sparse.csgraph.minimum_spanning_tree(matrix)
    # Now back to NK
    nk.setNumberOfThreads(ncores)
    mgraph = nk.Graph(mst.shape[0], weighted=True, directed=False)
    [mgraph.addEdge(x, y, 1) for x, y in zip(mst.nonzero()[0], mst.nonzero()[1])]
    # Now that we have the new minimum-spanning-tree graph we can get the diameter or the distance
    nodes = np.fromiter(mgraph.iterNodes(), np.int)
    if start is None or end is None:
        apsp = nk.distance.APSP(mgraph).run()
        distances = np.array(apsp.getDistances())
        # Set unavailable distances to -1
        distances[(distances == max_double) | np.isnan(distances)] = -1
        best = np.where(distances == distances.max())
        start, end = nodes[best[0][0]], nodes[best[1][0]]
    else:
        # Get the local coordinates
        start = cidx.query("cluster == @start").index.values[0]
        end = cidx.query("cluster == @end").index.values[0]

    path = np.array(nk.distance.BFS(mgraph, start, target=end).run().getPath(end), dtype=np.int)
    # Now convert to coordinates
    assert np.in1d(path, cidx).all()
    # Now do the 2-OPT optimisation to find the correct path
    path = tsp_2_opt(mst.toarray(), path)
    # Get a database of the shape "cluster", "bin"
    df = cidx.set_index("cidx").merge(
        pd.DataFrame().assign(cidx=path, bin=np.arange(path.shape[0], dtype=np.int)).set_index("cidx"),
        left_index=True, right_index=True).reset_index(drop=True)
    ranks = pd.DataFrame().assign(cluster=df["cluster"], rank=0)
    r = 0
    n = edge_list.loc[(~edge_list["cluster1"].isin(df["cluster"])) &
                      (edge_list["cluster2"].isin(df["cluster"])), "cluster1"].unique()

    while n.shape[0] > 0:
        r += 1
        tmp = edge_list.loc[(edge_list["cluster2"].isin(df["cluster"])) &
                            (edge_list["cluster1"].isin(n))].loc[
            lambda _tmp: ~_tmp["cluster1"].duplicated(keep="first")]
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
    # print("Final DF", df, sep="\n")
    # print("Path", path, sep="\n")
    # Remember we added the index earlier
    df.eval("backbone = cidx in @path", inplace=True)
    return df[["cluster", "bin", "rank", "backbone"]]
