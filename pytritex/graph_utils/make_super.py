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
from ..utils import _rebalance_ddf
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


def find_previous_results(raw_membership, previous_membership) -> (set, list):
    # TODO: first off, we need to format the previous_membership
    # previous_membership = previous_membership.merge(cidx, right_index=True, left_index=True)
    # TODO now we have to set up a dictionary
    if raw_membership.shape[0] == 0 or previous_membership.shape[0] == 0:
        return {}, []

    previous_membership = previous_membership.reset_index(drop=False).set_index("super")
    assert previous_membership.scaffold_index.unique().shape[0] == previous_membership.shape[0]
    _previous = previous_membership.groupby(level=0).apply(
        lambda group: tuple(sorted(group["scaffold_index"].values.tolist()))
    ).to_frame("key").reset_index(drop=False).set_index("key")
    _previous = _previous.rename(columns={"super": "previous_super"})

    _raw_mem = raw_membership.reset_index(drop=False).set_index("super")
    _raw = _raw_mem.groupby(level=0).apply(
        lambda group: tuple(sorted(group["cluster"].values.tolist()))
    ).to_frame("key").reset_index(drop=False).set_index("key")

    merged = _previous.merge(_raw, left_index=True, right_index=True, how="inner")
    to_skip = set(merged["super"].values.tolist())
    scaffolds_to_skip = set(_raw_mem.loc[to_skip]["cluster"].values)
    previous = previous_membership.loc[
        set(merged["previous_super"].values.tolist()),
        ["scaffold_index", "bin", "rank", "backbone"]].drop_duplicates(subset=["scaffold_index"])
    assert previous["rank"].max() <= 1
    assert scaffolds_to_skip == set(previous["scaffold_index"].values)
    logger.warning("Retaining %s groups, with %s scaffolds, from the previous iteration",
                   len(to_skip), len(scaffolds_to_skip))
    previous = previous.values
    return to_skip, [previous]


def make_super(hl: dd.DataFrame,
               cluster_info: dd.DataFrame,
               previous_membership: pd.DataFrame,
               client: Client,
               cores=1, paths=True, path_max=0, known_ends=False,
               maxiter=100, verbose=True):

    #  hl[cluster1 %in% cluster_info[excluded == F]$cluster & cluster2 %in% cluster_info[excluded == F]$cluster]->hl
    hl = hl.copy()
    non_excluded = cluster_info.loc[cluster_info["excluded"] == False].index.values.compute()
    logger.warning("Retaining %s scaffolds out of %s",
                   non_excluded.shape[0],
                   cluster_info.shape[0].compute())
    bait1 = hl["cluster1"].isin(non_excluded)
    bait2 = hl["cluster2"].isin(non_excluded)
    hl = hl.loc[bait1 & bait2, :].persist()
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

    if len(cidx_list) > len(set(cidx_list)):
        logger.error("The clustering into components is NOT unique.")
        import sys
        sys.exit(1)
    # Create a dataframe of cluster/super-scaffolds relationship
    num_supers = len(ssuper)
    length = len(cidx_list)
    raw_membership = pd.DataFrame().assign(cidx=cidx_list, super=ssuper).set_index("cidx")
    raw_membership = raw_membership.merge(
        cidx.reset_index(drop=False), on=["cidx"], how="outer").set_index("cluster")
    if raw_membership.shape[0] != length:
        logger.error("Duplicated indices after merging")
        import sys
        sys.exit(1)
    membership = cluster_info.merge(raw_membership, left_index=True,
                                    right_index=True, how="right")
    membership = membership.persist()
    if membership.shape[0].compute() != length:
        logger.error("Duplicated indices after merging with cluster_info")
        import sys
        sys.exit(1)

    # Where is each super-scaffold located?
    grouped = membership.groupby("super")
    info = grouped.agg(
        {"super": "size",
         "cM": "mean"}
    ).rename(columns={"super": "size"}).persist()
    # print(info.head(npartitions=-1, n=5))
    info["chr"] = grouped["chr"].unique().apply(lambda s: [_ for _ in s if not np.isnan(_)][0],
                                                meta=np.float)
    # print(info.head(npartitions=-1, n=5))
    assert membership.index.name == "cluster"
    edge_list = membership[["super"]]
    edge_list.index = edge_list.index.rename("cluster1")

    edge_list = edge_list.merge(hl, on="cluster1", how="right")
    edge_list = _rebalance_ddf(edge_list, npartitions=100)

    super_object = {"super_info": info,
                    "membership": membership, "graph": graph, "edges": edge_list}

    if paths is True:
        cms = membership.loc[:, ["super", "cM"]
              ].reset_index(drop=False).drop_duplicates().set_index("super").persist()
        cms = cms.repartition(npartitions=num_supers)
        # cms = _rebalance_ddf(cms, npartitions=100)
        if path_max > 0:
            idx = np.unique(super_object["super_info"][["super", "length"]].compute().sort_values(
                "length", ascending=False).head(path_max).index.values)
        else:
            idx = np.unique(super_object["super_info"].index.values.compute())
        edges = super_object["edges"].set_index("super").repartition(npartitions=num_supers)
        # grouped_edges = super_object["edges"].groupby("super")
        # grouped_membership = cms.groupby("super")
        order = membership
        # Removing the useless index (=cluster)
        order = order.loc[order["super"].isin(idx),
                          ["super", "chr"]].drop_duplicates().compute().sort_values(["chr"]).reset_index(drop=True)
        print(time.ctime(), "Starting super, with", idx.shape[0], "positions to analyse.")
        results = []
        to_skip, previous_results = find_previous_results(raw_membership, previous_membership)
        # First thing: let's check whether we have previous hits. We will consider these separately.
        # pool = mp.Pool(cores)
        # membership = client.scatter(membership)
        total = order.shape[0]
        analysed, cached = 0, 0
        for chrom in sorted(np.unique(order["chr"].values)):
            chrom_analysed, chrom_cached = 0, 0
            chrom_results = []
            subset = order[order["chr"] == chrom]
            chrom_total = subset.shape[0]
            logger.warning("%s Starting chromosome %s", time.ctime(), chrom)
            for row in subset.itertuples(name=None):
                index, ssuper, popseq_chr = row
                if ssuper in to_skip:
                    logger.debug("%s Using cached result for super %s on chromosome %s",
                                 time.ctime(), ssuper, popseq_chr)
                    chrom_cached += 1
                    continue
                chrom_analysed += 1
                logger.debug("%s Analysing super %s on chromosome %s", time.ctime(), ssuper, popseq_chr)
                my_edges = client.scatter(edges.loc[ssuper])
                my_membership = client.scatter(cms.loc[ssuper])
                chrom_results.append(client.submit(_concatenator,
                                                   my_edges, my_membership, known_ends, maxiter, verbose))
            logger.warning("%s Finished chr. %s (%s, %s%%), analysed %s, cached %s",
                           time.ctime(), chrom, chrom_total,
                           round(100 * chrom_total / total, 2),
                           chrom_analysed, chrom_cached
                           )
            analysed += chrom_analysed
            cached += chrom_cached
            chrom_results = client.gather(chrom_results)
            results.extend(chrom_results)

        logger.warning("%s Finished make_super_scaffolds", time.ctime())
        results = np.vstack(results + previous_results)
        results = pd.DataFrame().assign(
            cluster=results[:, 0], bin=results[:, 1], rank=results[:, 2], backbone=results[:, 3])
        assert results.cluster.unique().shape[0] == results.shape[0]
        results = results.set_index("cluster")
        super_object["membership"] = dd.merge(results, super_object["membership"], on="cluster")
        assert super_object["membership"].index.name == "cluster"

    super_object["membership"] = super_object["membership"].drop("cidx", axis=1)
    logger.warning("%s Finished make_super run", time.ctime())
    return super_object
