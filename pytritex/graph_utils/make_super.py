import pandas as pd
import numpy as np
from dask.delayed import delayed
import time
import networkit as nk
import dask.dataframe as dd
from pytritex.graph_utils.make_super_path import make_super_path
from dask.distributed import Client
import logging
from pytritex.utils import first
logger = logging.getLogger("distributed.worker")
import sys


def _concatenator(edges, membership, known_ends=False, maxiter=100, verbose=False):
    start = end = None
    if known_ends:
        assert "cluster" in membership
        x = membership.query("cM == cM").sort_values("cM")["cluster"].values
        start, end = x[0], x[-1]
    final = make_super_path(edges, membership, start=start, end=end, maxiter=maxiter, verbose=verbose)
    return final


def add_statistics(membership, client):
    # m[,.(n=.N, nbin=max(bin), max_rank = max(rank), length = sum(length)), key = super]->res
    # res[,.(super, super_size=n, super_nbin=nbin)][m, on = "super"]->m

    # Note that membership has *no index* here.
    membership = membership.drop(["super_size", "super_nbin"], axis=1, errors="ignore")
    mem_sup_group = membership.groupby("super")
    res_size = mem_sup_group.size().to_frame("n")
    res_cols = mem_sup_group.agg(
        {"bin": "max", "rank": "max", "length": "sum"})
    res_cols.columns = ["nbin", "max_rank", "length"]
    res = dd.merge(res_size, res_cols, on="super")
    left = res.loc[:, ["n", "nbin"]].rename(
        columns={"n": "super_size", "nbin": "super_nbin"})
    if left.index.dtype != membership.index.dtype:
        left = left.reset_index(drop=False)
        left["super"] = left["super"].astype(membership.index.dtype)
        left = left.set_index("super")
    func = delayed(dd.merge)(left, membership, left_index=True, right_on="super",
                             how="right")
    membership = client.compute(func).result()
    assert "scaffold_index" in membership.columns
    membership = membership.set_index("scaffold_index")
    if not isinstance(membership, dd.DataFrame):
        membership = dd.from_pandas(membership, chunksize=10**5)
    if not isinstance(res, dd.DataFrame):
        res = dd.from_pandas(res, chunksize=10**5)
    return membership, res



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

    hl = hl.copy()
    non_excluded = cluster_info.loc[cluster_info["excluded"] == False].index.values.compute()
    logger.warning("Retaining %s scaffolds out of %s",
                   non_excluded.shape[0],
                   cluster_info.shape[0].compute())
    bait1 = hl["cluster1"].isin(non_excluded)
    bait2 = hl["cluster2"].isin(non_excluded)
    hl = hl.loc[bait1 & bait2, :]
    try:
        edge_list = hl.query("cluster1 < cluster2")[["cluster1", "cluster2", "weight"]].compute().drop_duplicates()
    except ValueError:
        logger.critical(hl.dtypes)
        logger.critical(hl.head(npartitions=-1))
        raise

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
    edge_list = edge_list.reset_index(drop=True)
    # Now we are ready to create the graph using the indices.
    nk.setNumberOfThreads(cores)
    graph = nk.Graph(n=cidx.shape[0], weighted=True, directed=False)
    edge_list.apply(lambda row: graph.addEdge(row["cidx1"], row["cidx2"], row["weight"]), axis=1)
    components = nk.components.ConnectedComponents(graph)
    # Get the connected componentsx
    components.run()
    cidx_list, ssuper = [], []
    # Aggregate into super-scaffolds
    for idx, comp in enumerate(components.getComponents()):
        cidx_list.extend(comp)
        ssuper.extend([idx] * len(comp))

    if len(cidx_list) > len(set(cidx_list)):
        logger.error("The clustering into components is NOT unique.")
        sys.exit(1)
    # Create a dataframe of cluster/super-scaffolds relationship
    num_supers = len(ssuper)
    length = len(cidx_list)
    raw_membership = pd.DataFrame().assign(cidx=cidx_list, super=ssuper).astype(
        {"cidx": int, "super": int}).set_index("cidx")
    raw_membership = raw_membership.merge(
        cidx.reset_index(drop=False), on=["cidx"], how="outer").set_index("cluster")
    if raw_membership.shape[0] != length:
        logger.error("Duplicated indices after merging")
        sys.exit(1)
    membership = cluster_info.merge(raw_membership, left_index=True,
                                    right_index=True, how="right", npartitions=cluster_info.npartitions)
    membership = membership.reset_index(drop=False).astype({"cluster": int}).set_index("cluster").persist()
    if membership.shape[0].compute() != length:
        logger.error("Duplicated indices after merging with cluster_info")
        sys.exit(1)

    # Where is each super-scaffold located?
    grouped = membership.groupby("super")
    info = grouped.agg(
        {"super": "size",
         "cM": "mean",
         "length": "sum",
         "chr": first}
    ).rename(columns={"super": "size"})
    assert membership.index.name == "cluster"
    edge_list = membership[["super"]]
    edge_list.index = edge_list.index.rename("cluster1")
    hl = hl.astype({"cluster1": int, "cluster2": int}).set_index("cluster1").persist()
    edge_list = edge_list.merge(hl, on="cluster1", how="right").reset_index(drop=False).persist()
    assert edge_list.query("super != super").shape[0].compute() == 0  # Check no undefined row is present
    # Check we have the same supers
    assert (edge_list["super"].unique().compute() == membership["super"].unique().compute()).all()
    
    super_object = {"super_info": info,
                    "membership": membership, "graph": graph, "edges": edge_list}

    if paths is True:
        cms = membership.loc[:, ["super", "chr", "cM"]].copy().reset_index(drop=False).set_index("super").rename(
                  columns={"scaffold_index": "cluster"})
        if "cluster" not in cms.columns:
            logger.critical("\n\n%s\n\n", cms.head())
            sys.exit(1)
        assert "cluster" in cms.columns, cms.head()
        if path_max > 0:
            assert info.index.name == "super"
            info["max_length"] = info.groupby("chr")["length"].transform("max", meta=int)
            info["is_max"] = info["max_length"] == info["length"]
            idx = info.query("is_max == True & chr == chr").index.values.compute()
            assert idx.shape[0] == path_max
        else:
            idx = np.unique(super_object["super_info"].index.values.compute())
        edges = super_object["edges"].set_index("super")

        order = membership
        # Removing the useless index (=cluster)
        order = order.loc[order["super"].isin(idx),
                          ["super", "chr"]].drop_duplicates().compute().sort_values(["chr"]).reset_index(drop=True)
        to_skip, previous_results = find_previous_results(raw_membership, previous_membership)
        # First thing: let's check whether we have previous hits. We will consider these separately.
        results = []
        submitted = []
        nrows = order[~order.super.isin(to_skip)].shape[0]
        if nrows > 0:
            logger.warning("%s Starting to analyse %s rows", time.ctime(), order[~order.super.isin(to_skip)].shape[0])
            for row in order[~order.super.isin(to_skip)].itertuples():
                index, ssuper, popseq_chr = row
                my_edges = edges.loc[ssuper].compute()
                my_membership = cms.loc[ssuper].compute()
                submitted.append(client.submit(_concatenator, my_edges, my_membership, known_ends, maxiter, verbose))
            logger.warning("%s Retrieving final results", time.ctime())
            submitted = client.gather(submitted)
            assert isinstance(submitted[0], np.ndarray), (type(submitted[0]), submitted[0])
            results.extend(submitted)
            del submitted
            logger.warning("%s Finished make_super_scaffolds", time.ctime())
            results = np.vstack(results + previous_results)
            results = pd.DataFrame().assign(
                cluster=results[:, 0], bin=results[:, 1], rank=results[:, 2], backbone=results[:, 3])
            assert results.cluster.unique().shape[0] == results.shape[0]
            results = results.set_index("cluster")
            super_object["membership"] = dd.merge(results, super_object["membership"], on="cluster")
            if not isinstance(super_object["membership"], dd.DataFrame):
                super_object["membership"] = dd.from_pandas(super_object["membership"], chunksize=10 ** 5)
            assert super_object["membership"].index.name == "cluster"
            super_object["membership"] = super_object["membership"].drop("cidx", axis=1).persist()
        else:
            super_object["membership"] = dd.merge(
                pd.DataFrame().assign(cluster=[], bin=[], rank=[], backbone=[]),
                super_object["membership"], on="cluster")
            super_object["membership"] = super_object["membership"].set_index("cluster")

    logger.warning("%s Finished make_super run", time.ctime())
    logger.warning("%s Membership size and dtypes: %s, %s", time.ctime(), super_object["membership"].shape[0].compute(),
                   super_object["membership"].dtypes)
    return super_object
