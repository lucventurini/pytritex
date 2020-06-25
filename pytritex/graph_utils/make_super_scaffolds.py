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


def make_super_scaffolds(links: Union[str, dd.DataFrame],
                         info: Union[str, dd.DataFrame],
                         membership: Union[str, dd.DataFrame],
                         save_dir: str,
                         client: Client,
                         excluded=None,
                         to_parquet=True,
                         ncores=1):
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

    if membership is not None:
        if isinstance(membership, str):
            membership = dd.read_parquet(membership, infer_divisions=True)
        else:
            assert isinstance(membership, dd.DataFrame)
        # cluster=results[:, 0], bin=results[:, 1], rank=results[:, 2], backbone=results[:, 3])
        membership = membership.loc[
            membership.super_size > 1, ["super", "bin", "rank", "backbone"]].compute()
        membership["max_rank"] = membership.groupby("super")["rank"].transform("max")
        membership = membership.loc[membership["max_rank"] <= 1].drop("max_rank", axis=1)
    else:
        membership = pd.DataFrame(
            columns=["scaffold_index", "super", "bin", "rank", "backbone"]).set_index("scaffold_index")

    info2 = info.loc[:, ["popseq_chr", "popseq_cM", "length"]].rename(
        columns={"popseq_chr": "chr", "popseq_cM": "cM"})
    assert "popseq_chr" not in info2.columns and "chr" in info2.columns
    assert "popseq_cM" not in info2.columns and "cM" in info2.columns
    if excluded is not None:
        excluded_scaffolds = pd.Series(list(excluded), name="scaffold_index")
    else:
        excluded = pd.Series([], name="scaffold_index")
        excluded_scaffolds = excluded.copy()
    input_df = info2.assign(excluded=info.index.isin(excluded_scaffolds))
    input_df.index = input_df.index.rename("cluster")
    hl = links.copy().rename(columns={"scaffold_index1": "cluster1", "scaffold_index2": "cluster2"})
    super_scaffolds = make_super(
        hl=hl,
        cluster_info=input_df,
        previous_membership=membership,
        client=client,
        verbose=False, cores=ncores,
        paths=True, path_max=0, known_ends=False, maxiter=100)
    mem_copy = super_scaffolds["membership"]
    assert mem_copy.index.name == "cluster"
    mem_copy.index = mem_copy.index.rename("scaffold_index")
    maxidx = super_scaffolds["super_info"]["super"].max().compute()
    assert mem_copy.index.name == "scaffold_index"
    iindex = info.index.values.compute()
    bait_index = ~np.in1d(iindex, mem_copy.index.values.compute())
    bait = iindex[bait_index]

    _to_concatenate = info.loc[bait, ["popseq_chr", "popseq_cM", "length"]].rename(
        columns={"popseq_chr": "chr", "popseq_cM": "cM"})

    _to_concatenate = client.scatter(_to_concatenate)
    chunks = delayed(lambda df: df.map_partitions(len))(_to_concatenate)
    chunks = tuple(client.compute(chunks).result().compute().values.tolist())
    _to_concatenate = client.gather(_to_concatenate)

    sup_column = da.from_array(
        maxidx + pd.Series(range(1, info.loc[bait].shape[0].compute() + 1)), chunks=chunks)
    excl_column = da.from_array(np.in1d(bait, excluded_scaffolds), chunks=chunks)

    _to_concatenate = _to_concatenate.assign(
        bin=1, rank=0, backbone=True,
        excluded=excl_column,
        super=sup_column)
    _to_concatenate = client.scatter(_to_concatenate)
    mem_copy = client.scatter(mem_copy)
    func = delayed(dd.concat)([mem_copy, _to_concatenate])
    mem_copy = client.compute(func).result()
    mem_sup_group = mem_copy.groupby("super")
    res_size = mem_sup_group.size().to_frame("n")
    res_cols = mem_sup_group.agg(
        {"bin": "max", "rank": "max", "length": "sum"})
    res_cols.columns = ["nbin", "max_rank", "length"]
    res = dd.merge(res_size, res_cols, on="super")
    left = res.loc[:, ["n", "nbin"]].rename(
        columns={"n": "super_size", "nbin": "super_nbin"})
    left = client.scatter(left)
    mem_copy = client.scatter(mem_copy)
    func = delayed(dd.merge)(left, mem_copy, left_index=True, right_on="super",
                             how="right")
    mem_copy = client.compute(func).result()
    mem_copy_iname = mem_copy.index.name
    mem_copy = mem_copy.reset_index(drop=False).set_index(mem_copy_iname)
    # mem_copy = dd.from_pandas(mem_copy, chunksize=1000)
    if to_parquet is True:
        dd.to_parquet(mem_copy, os.path.join(save_dir, "membership"))
        # res = dd.from_pandas(res, chunksize=1000)
        dd.to_parquet(res, os.path.join(save_dir, "result"))
        return {"membership": os.path.join(save_dir, "membership"),
                "info": os.path.join(save_dir, "result")}
    else:
        return {"membership": mem_copy, "info": res}
