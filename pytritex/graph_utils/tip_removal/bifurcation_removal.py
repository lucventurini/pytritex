from pytritex.graph_utils.make_super_scaffolds import make_super_scaffolds
from dask.distributed import Client
from . import minimum_distance, _calculate_degree
from .bulge_removal import _remove_bulges
import dask.dataframe as dd
import numpy as np
from dask.delayed import delayed
import pandas as pd


def _remove_bifurcations(links: dd.DataFrame,
                         excluded,
                         membership: dd.DataFrame,
                         info: dd.DataFrame,
                         client: Client, save_dir: str,
                         min_dist=minimum_distance, ncores=1):
    # resolve length-one-bifurcations at the ends of paths
    #    m[rank > 0][bin == 2 | super_nbin - 1 == bin ][,
    #    .(super, super_nbin, type = bin == 2, scaffold, length, bin0=bin)]->x
    #    unique(rbind(
    #    m[x[type == T, .(super, bin0, bin=1)], on=c("super", "bin")],
    #    m[x[type == F, .(super, bin0, bin=super_nbin)], on=c("super", "bin")]
    #    ))->a
    #    a[, .(super, bin0, scaffold2=scaffold, length2=length)][
    #    x, on=c("super", "bin0")][, ex := ifelse(length >= length2, scaffold2, scaffold)]$ex -> add
    #
    #    if(length(add) > 0){
    #     ex <- c(ex, add)
    #     make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
    #     out$membership -> m
    #    }


    if membership.shape[0].compute() == 0:
        return membership, info, excluded
    keys = ["super", "super_nbin", "length", "bin"]
    bifurcated = membership.query("(rank > 0) & ((bin == 2) | (super_nbin - 1 == bin) )")
    if bifurcated.shape[0].compute() == 0:
        return bifurcated, info, excluded
    bifurcated = bifurcated[keys]
    bifurcated = bifurcated.eval("type = (bin == 2)").rename(columns={"bin": "bin0"})
    bifurcated = bifurcated[["super", "super_nbin", "type", "length", "bin0"]]
    bifurcated = bifurcated.reset_index(drop=False).drop_duplicates()
    key = ["super", "bin"]
    # indexed = dd_membership.set_index(key)
    # m[x[type == T,.(super, bin0, bin=1)], on = c("super", "bin")],
    membership = client.scatter(membership.reset_index(drop=False))
    upper = client.scatter(bifurcated.loc[bifurcated.type, ["super", "bin0"]].assign(bin=1))
    func = delayed(dd.merge)(membership, upper, on=key)
    upper = client.compute(func).result()
    upper = upper.persist()
    assert "scaffold_index" in upper.columns
    upper = client.scatter(upper)
    # m[x[type == F,.(super, bin0, bin=super_nbin)], on = c("super", "bin")]
    lower = client.scatter(bifurcated.loc[~bifurcated.type, ["super", "bin0", "super_nbin"]].rename(
        columns={"super_nbin": "bin"}))
    func = delayed(dd.merge)(membership, lower, on=key)
    lower = client.compute(func).result()
    lower.persist()
    lower = client.scatter(lower)
    assert "scaffold_index" in lower.columns
    # assert lower.index.name == "scaffold_index"
    func = delayed(dd.concat)([upper, lower])
    a = client.compute(func).result()
    a = a.drop_duplicates().rename(columns={"length": "length2", "scaffold_index": "scaffold_index2"})
    a = client.scatter(a)
    # Now merge with "x" to find places to exclude
    bifurcated = client.scatter(bifurcated)
    func = delayed(dd.merge)(a, bifurcated, on=["super", "bin0"], how="right")
    a = client.compute(func).result().compute()
    add = np.where((a.length >= a.length2), a.scaffold_index2, a.scaffold_index)
    out = {"membership": membership, "info": None}
    if add.shape[0] > 0:
        excluded.update(add.tolist())
        out = make_super_scaffolds(links=links, info=info,
                                   membership=membership,
                                   excluded=excluded, ncores=ncores,
                                   client=client, save_dir=save_dir,
                                   to_parquet=False)
        membership = out["membership"]
    membership, res, excluded = _remove_bulges(links=links, excluded=excluded,
                                               membership=membership,
                                               save_dir=save_dir,
                                               client=client,
                                               info=info, min_dist=min_dist, ncores=ncores)
    if res is not None:
        out["info"] = res
    return membership, out["info"], excluded
