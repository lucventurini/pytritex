import pandas as pd
from pytritex.graph_utils.make_super import make_super
import dask.dataframe as dd
import os
from dask.distributed import Client
import numpy as np
import dask.array as da
from dask.delayed import delayed
import logging
from typing import Union
logger = logging.getLogger("distributed.worker")


def get_previous_groups(membership):
    """This function will retrieve the previous groups to which scaffolds have been assigned.
    We only consider groups where the maximum rank is 1.
    This has the purpose of avoiding recomputing groups for nothing during the TSP stage."""

    if membership is not None:
        if isinstance(membership, str):
            membership = dd.read_parquet(membership, infer_divisions=True)
        else:
            assert isinstance(membership, dd.DataFrame)
        membership = membership.loc[
            membership.super_size > 1, ["super", "bin", "rank", "backbone"]].compute()
        membership["max_rank"] = membership.groupby("super")["rank"].transform("max")
        membership = membership.loc[membership["max_rank"] <= 1].drop("max_rank", axis=1)
    else:
        membership = pd.DataFrame(
            columns=["scaffold_index", "super", "bin", "rank", "backbone"]).set_index("scaffold_index")
    return membership


def prepare_tables(links, info, membership, excluded):
    if isinstance(links, str):
        links = dd.read_parquet(links, infer_divisions=True)
    else:
        assert isinstance(links, dd.DataFrame)
    if excluded is None:
        logger.warning("No excluded scaffolds.")
    else:
        logger.warning("No. of excluded scaffolds: %s", len(excluded))
    if isinstance(info, str):
        info = dd.read_parquet(info, infer_divisions=True)
    else:
        assert isinstance(info, dd.DataFrame)
    membership = get_previous_groups(membership)
    info2 = info.loc[:, ["popseq_chr", "popseq_cM", "length"]].rename(
        columns={"popseq_chr": "chr", "popseq_cM": "cM"})
    assert "popseq_chr" not in info2.columns and "chr" in info2.columns
    assert "popseq_cM" not in info2.columns and "cM" in info2.columns
    if excluded is not None:
        excluded_scaffolds = pd.Series(list(excluded), name="scaffold_index")
    else:
        excluded = pd.Series([], name="scaffold_index")
        excluded_scaffolds = excluded.copy()
    cluster_info = info2.assign(excluded=info.index.isin(excluded_scaffolds))
    cluster_info.index = cluster_info.index.rename("cluster")
    hl = links.copy().rename(columns={"scaffold_index1": "cluster1", "scaffold_index2": "cluster2"})
    return links, info, membership, excluded_scaffolds, cluster_info, hl


def add_missing_scaffolds(info, membership, maxidx, excluded_scaffolds, client):
    # _to_concatenate <- info[!m$scaffold, on = "scaffold"]
    # _to_concatenate <- _to_concatenate[,.(
    # scaffold, bin=1, rank=0, backbone=T, chr=popseq_chr, cM=popseq_cM, length=length,
    # excluded=scaffold % in % excluded_scaffolds, super = paste0(
    #     prefix, "_", maxidx + 1:.N))])
    # rbind(m, [
    #   ->m

    indices = membership.index.values.compute()
    if indices.shape[0] > len(set(indices)):
        logger.error("Something went wrong in the calculation, we have duplicated indices.")
        import sys
        sys.exit(1)

    assert membership.index.name == "scaffold_index"
    # TODO: this part of the code is probably causing the double counting
    info_index = info.index.values.compute()
    bait_index = sorted(set.difference(set(info_index), set(indices)))
    # bait = info_index[bait_index]
    # # Now of course the sum of bait and indices should be equal to info_index.
    # if bait.shape[0] + indices.shape[0] != info_index.shape[0]:
    #     if bait.shape[0] + indices.shape[0] > info_index.shape[0]:
    #         logger.error("Something went wrong in getting the bait, we have duplicated indices.")
    #     else:
    #         logger.error("Something went wrong in getting the bait, we have missing indices.")
    #     import sys
    #     sys.exit(1)

    _to_concatenate = info.loc[bait_index, ["popseq_chr", "popseq_cM", "length"]].rename(
        columns={"popseq_chr": "chr", "popseq_cM": "cM"})
    assert _to_concatenate.shape[0].compute() == len(bait_index)

    chunks = client.compute(delayed(
        lambda df: df.map_partitions(len))(client.scatter(_to_concatenate))).result()
    chunks = tuple(chunks.compute().values.tolist())

    sup_column = da.from_array(
        maxidx + pd.Series(range(1, len(bait_index) + 1)), chunks=chunks)

    excl_vector = np.in1d(bait_index, excluded_scaffolds)
    if excl_vector.shape[0] != len(bait_index):
        logger.error("Wrong length of excl_vector")
        import sys
        sys.exit(1)
    excl_column = da.from_array(excl_vector, chunks=chunks)

    _to_concatenate = _to_concatenate.assign(
        bin=1, rank=0, backbone=True,
        excluded=excl_column,
        super=sup_column)
    func = delayed(dd.concat)([client.scatter(membership), client.scatter(_to_concatenate)])
    new_membership = client.compute(func).result()
    new_membership = new_membership.reset_index(drop=False)
    assert "scaffold_index" in new_membership.columns
    sis = new_membership["scaffold_index"].compute()
    if sis.shape[0] > len(set(sis.values)):
        logger.error("""Duplicated indices after concatenating.
Original size: %s
Membership: %s
To add: %s
Concatenated: %s""", info_index.shape[0], indices.shape[0], _to_concatenate.shape[0].compute(),
                     new_membership.shape[0].compute())
        import sys
        sys.exit(1)

    new_membership = new_membership.persist()
    return new_membership


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
    left = client.scatter(left)
    membership = client.scatter(membership)
    func = delayed(dd.merge)(left, membership, left_index=True, right_on="super",
                             how="right")
    membership = client.compute(func).result()
    assert "scaffold_index" in membership.columns
    membership = membership.set_index("scaffold_index")
    return membership, res


def make_super_scaffolds(links: Union[str, dd.DataFrame],
                         info: Union[str, dd.DataFrame],
                         membership: Union[str, dd.DataFrame],
                         save_dir: str,
                         client: Client,
                         excluded=None,
                         to_parquet=True,
                         ncores=1):

    # Links, info are as they were.
    # "Membership", if given, provides the pre-calculated groupings.
    # Excluded scaffolds is the list of scaffolds to ignore.
    # cluster_info is a dataframe of the form:
    # chr, cM, length, excluded
    # index is cluster == scaffold_index
    # HL: is the list of links to consider

    links, info, membership, excluded_scaffolds, cluster_info, hl = prepare_tables(
        links, info, excluded=excluded, membership=membership)

    super_scaffolds = make_super(
        hl=hl,
        cluster_info=cluster_info,
        previous_membership=membership,
        client=client,
        verbose=False, cores=ncores,
        paths=True, path_max=0, known_ends=False, maxiter=100)

    # Take the membership table. Extract the scaffolds IDs which are *not* present in the table.
    membership = super_scaffolds["membership"]
    assert membership.index.name == "cluster"
    membership.index = membership.index.rename("scaffold_index")
    maxidx = super_scaffolds["super_info"]["super"].max().compute()
    membership = add_missing_scaffolds(info, membership, maxidx, excluded_scaffolds, client)
    membership, res = add_statistics(membership, client)
    # mem_copy = dd.from_pandas(mem_copy, chunksize=1000)
    if to_parquet is True:
        dd.to_parquet(membership, os.path.join(save_dir, "membership"))
        # res = dd.from_pandas(res, chunksize=1000)
        dd.to_parquet(res, os.path.join(save_dir, "result"))
        return {"membership": os.path.join(save_dir, "membership"),
                "info": os.path.join(save_dir, "result")}
    else:
        return {"membership": membership, "info": res}
