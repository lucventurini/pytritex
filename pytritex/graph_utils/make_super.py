import pandas as pd
import numpy as np
# import graph_tool
# import graph_tool as gt
# import graph_tool.clustering
import time
import networkit as nk
import dask.dataframe as dd
from pytritex.graph_utils.make_super_path import make_super_path
from dask.distributed import Client
from dask.delayed import delayed
import logging
logger = logging.getLogger("distributed.worker")


def _concatenator(edges, membership, known_ends=False, maxiter=100, verbose=False):
    start = end = None
    if known_ends:
        x = membership.query("cM == cM").sort_values("cM").index
        start, end = x.head(1), x.tail(1)
    final = make_super_path(edges.compute(), membership.compute(),
                            start=start, end=end,
                            maxiter=maxiter, verbose=verbose)
    return final


def make_super(hl: dd.DataFrame,
               cluster_info: dd.DataFrame,
               previous_membership: pd.DataFrame,
               client: Client,
               cores=1, paths=True, path_max=0, known_ends=False,
               maxiter=100, verbose=True):

    #  hl[cluster1 %in% cluster_info[excluded == F]$cluster & cluster2 %in% cluster_info[excluded == F]$cluster]->hl
    hl = hl.copy()
    non_excluded = cluster_info.loc[cluster_info["excluded"] == False].index.values.compute()
    logger.info("Retaining %s scaffolds out of %s",
                non_excluded.shape[0],
                cluster_info.shape[0].compute())
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
    nk.setNumberOfThreads(cores)
    graph = nk.Graph(n=cidx.shape[0], weighted=True, directed=False)
    edge_list.apply(lambda row: graph.addEdge(row["cidx1"], row["cidx2"], row["weight"]), axis=1)
    components = nk.components.ConnectedComponents(graph)
    # Get the connected components
    components.run()
    cidx_list, ssuper = [], []
    # Aggregate into super-scaffolds
    for idx, comp in enumerate(components.getComponents()):
        cidx_list.extend(comp)
        ssuper.extend([idx] * len(comp))

    # Create a dataframe of cluster/super-scaffolds relationship
    membership = pd.DataFrame().assign(cidx=cidx_list, super=ssuper).set_index("cidx")
    membership = membership.merge(cidx.reset_index(drop=False), on=["cidx"], how="outer")
    membership = cluster_info.merge(membership, on="cluster", how="right")
    membership = membership.persist()

    # Where is each super-scaffold located?
    grouped = membership.groupby("super")
    info = grouped.agg(
        {"super": "size",
         "cM": "mean"}
    )
    # print(info.head(npartitions=-1, n=5))
    info["chr"] = grouped["chr"].unique().apply(lambda s: [_ for _ in s if not np.isnan(_)][0],
                                                meta=np.float)
    # print(info.head(npartitions=-1, n=5))
    assert "cluster" in membership, membership.columns
    edge_list = membership.rename(columns={"cluster": "cluster1"})[["cluster1", "super"]].merge(
        hl, on="cluster1", how="right")

    super_object = {"super_info": info, "membership": membership, "graph": graph, "edges": edge_list}

    if paths is True:
        # TODO: first off, we need to format the previous_membership
        # previous_membership = previous_membership.merge(cidx, right_index=True, left_index=True)
        # TODO now we have to set up a dictionary
        previous_membership = previous_membership.reset_index(drop=False).set_index("super")
        _mem = previous_membership.groupby(level=0).apply(
            lambda group: tuple(sorted(group["scaffold_index"].values.tolist()))
        ).to_frame("key").reset_index(drop=False)
        lookup = dict((_.key, _.super) for _ in _mem.itertuples(name="pandas", index=False))

        cms = membership.loc[
              :, ["cluster", "super", "cM"]].drop_duplicates().set_index("cluster")

        if path_max > 0:
            idx = super_object["super_info"].sort_values(
                "length", ascending=False).head(path_max)["super"].unique()
        else:
            idx = super_object["super_info"]["super"].unique()
        idx = idx.compute()
        grouped_edges = super_object["edges"].groupby("super")
        grouped_membership = cms.groupby("super")
        order = membership
        order = order.loc[order["super"].isin(idx),
                          ["super", "chr"]].drop_duplicates().compute().sort_values(["chr"])
        print(time.ctime(), "Starting super, with", idx.shape[0], "positions to analyse.")
        results = []
        previous_results = []
        # pool = mp.Pool(cores)
        # membership = client.scatter(membership)
        for row in order.itertuples(name=None):
            index, ssuper, popseq_chr = row
            print(time.ctime(), "Analysing", index, "super", ssuper, "on chromosome", popseq_chr)
            indices = grouped_edges.get_group(ssuper).index
            # TODO: this is the point. IF we have the result already in the previous membership
            # TODO: then we should get those out here.
            # TODO: we also have to make sure that the index associated with the super scaffold
            # TODO: is changed, otherwise we risk ID collisions.
            key = tuple(sorted(indices.compute().values.tolist()))
            if key in lookup:
                old_super = lookup[key]
                old_result = previous_membership.loc[
                    old_super, ["scaffold_index", "bin", "rank", "backbone"]].values
                assert old_result.shape[0] == len(key)
                previous_results.append(old_result)
            else:
                my_edges = client.scatter(edge_list.loc[indices])
                my_membership = client.scatter(
                    cms.loc[grouped_membership.get_group(ssuper).index])
                results.append(client.submit(_concatenator,
                    my_edges, my_membership, known_ends, maxiter, verbose))

        results = client.gather(results)
        results = np.vstack(results + previous_results)
        results = pd.DataFrame().assign(
            cluster=results[:, 0], bin=results[:, 1], rank=results[:, 2], backbone=results[:, 3])
        super_object["membership"] = dd.merge(results, super_object["membership"], on="cluster")
        assert "cluster" in super_object["membership"].columns

    print(time.ctime(), "Finished make_super run")
    # import sys
    # sys.exit(1)
    return super_object
