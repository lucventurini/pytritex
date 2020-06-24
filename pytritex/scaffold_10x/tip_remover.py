import pandas as pd
from pytritex.graph_utils.make_super_scaffolds import make_super_scaffolds
import dask.dataframe as dd
import numpy as np
from time import ctime
from dask.distributed import Client


minimum_distance = 1e4


def _calculate_degree(links: str, excluded) -> pd.DataFrame:
    links = dd.read_parquet(links, infer_divisions=True)
    degree = links[~links["scaffold_index2"].isin(excluded)].groupby(
        "scaffold_index1").size().to_frame("degree").compute()
    degree.index = degree.index.rename("scaffold_index")
    return degree


def _remove_bulges(links: str, excluded,
                   membership: str,
                   client: Client,
                   save_dir: str,
                   info: str, min_dist=minimum_distance, ncores=1):
    membership = dd.read_parquet(membership, infer_divisions=True)
    add = membership[(membership["rank"] == 1) & (membership["length"] <= min_dist)].index.values.compute()
    if add.shape[0] > 0:
        excluded.update(add.tolist())
        out = make_super_scaffolds(links=links, client=client, save_dir=save_dir,
                                   info=info, excluded=excluded, ncores=ncores)
        membership = out["membership"]
        res = out["info"]
    else:
        res = None
    return membership, res, excluded


def _remove_short_tips(links: str, excluded, membership: str, info: str,
                       client: Client, save_dir: str,
                       min_dist=minimum_distance, ncores=1):
    # #   links[!scaffold2 %in% ex][, .(degree=.N), key=.(scaffold_index=scaffold1)]-> degree
    # # inner <- m[rank == 1][, .(super=c(super, super, super), bin=c(bin, bin-1, bin+1))]
    # # middle <- m[inner, on=c("super", "bin")]
    # # a <- degree[middle, on="scaffold_index"]
    # #  add <-  a[degree == 1 & length <= 1e4]$scaffold
    # #   if(length(add) > 0){
    # #    ex <- c(ex, add)
    # #    make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
    # #    out$m -> m
    # #   }
    # #   # remove short tips/bulges of rank 1
    # #   m[rank == 1 & length <= 1e4]$scaffold -> add
    # #   ex <- c(ex, add)
    # #   if(length(add) > 0){
    # #    ex <- c(ex, add)
    # #    make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
    # #    out$m -> m
    # #   }

    membership = dd.read_parquet(membership, infer_divisions=True)
    inner = membership.loc[membership["rank"] == 1, ["super", "bin"]].values.compute()
    inner0 = np.tile(inner[:, 0], 3)
    inner1 = np.tile(inner[:, 1], 3).reshape(3, inner.shape[0]) + np.array([0, -1, 1]).reshape(3, 1)
    inner1 = inner1.flatten()
    inner = pd.DataFrame().assign(super=inner0, bin=inner1)
    middle = dd.merge(membership, inner, how="right",
                      on=["super", "bin"]).reset_index(drop=False).set_index("scaffold_index")
    links = dd.read_parquet(links, infer_divisions=True)
    degree = _calculate_degree(links, excluded)
    a = degree.merge(middle, on="scaffold_index").reset_index(drop=False)
    add = a.query("(degree == 1) & (length <= @min_dist)")["scaffold_index"].values
    out = {"membership": membership, "info": None}
    if add.shape[0] > 0:
        excluded.update(add.tolist())
        out = make_super_scaffolds(links=links, info=info, client=client,
                                   save_dir=save_dir,
                                   excluded=excluded, ncores=ncores)
        membership = out["membership"]

    info = dd.read_parquet(info, infer_divisions=True)
    membership, res, excluded = _remove_bulges(links=links, excluded=excluded,
                                               membership=membership,
                                               client=client, save_dir=save_dir,
                                               info=info, min_dist=min_dist, ncores=ncores)
    if res is not None:
        out["info"] = res

    return membership, out["info"], excluded


def _remove_bifurcations(links, excluded, membership, info,
                         client: Client, save_dir: str,
                         min_dist=minimum_distance, ncores=1):
    # resolve length-one-bifurcations at the ends of paths
    if membership.shape[0] == 0:
        return membership, info, excluded
    keys = ["super", "super_nbin", "scaffold_index", "length", "bin"]
    x = membership.query("(rank > 0) & ( (bin == 2) | (super_nbin - 1 == bin) )")
    if x.shape[0] == 0:
        return x, info, excluded
    x = x[keys]
    x = x.eval("type = (bin == 2)").rename(columns={"bin": "bin0"})
    x = x[["super", "super_nbin", "type", "scaffold_index", "length", "bin0"]]
    key = ["super", "bin"]
    indexed = membership.set_index(key)
    # m[x[type == T,.(super, bin0, bin=1)], on = c("super", "bin")],
    upper = indexed.merge(
        x.loc[x.type, ["super", "bin0"]].assign(bin=1).set_index(key),
        on=key)
    # m[x[type == F,.(super, bin0, bin=super_nbin)], on = c("super", "bin")]
    lower = indexed.merge(
        x.loc[~x.type, ["super", "bin0", "super_nbin"]].rename(columns={"super_nbin": "bin"}).set_index(key),
        on=key)
    a = pd.concat([upper, lower]).reset_index(drop=False).drop_duplicates().rename(
        columns={"scaffold_index": "scaffold_index2", "length": "length2"},
        errors="raise").set_index(["super", "bin0"])
    # Now merge with "x" to find places to exclude
    try:
        a = a.merge(x.set_index(["super", "bin0"]), on=["super", "bin0"], how="right")
    except ValueError:
        print(a.head())
        print(a.shape[0])
        print(a.index.get_level_values(0).dtype, a.index.get_level_values(1).dtype)
        print(x.head())
        print(x.shape[0])
        print(x[["super", "bin0"]].dtypes)
    add = np.where((a.length >= a.length2), a.scaffold_index2, a.scaffold_index)
    out = {"membership": membership, "info": None}
    if add.shape[0] > 0:
        excluded.update(add.tolist())
        out = make_super_scaffolds(links=links, info=info, excluded=excluded, ncores=ncores,
                                   client=client, save_dir=save_dir)
        membership = out["membership"]
    membership, res, excluded = _remove_bulges(links=links, excluded=excluded,
                                               membership=membership,
                                               save_dir=save_dir,
                                               client=client,
                                               info=info, min_dist=min_dist, ncores=ncores)
    if res is not None:
        out["info"] = res
    return membership, out["info"], excluded


def remove_tips(links: str, excluded, out: dict, info: str,
                client: Client,
                save_dir: str,
                ncores=1,
                verbose=False, min_dist=minimum_distance):

    if verbose:
        print(ctime(), "Starting tip removal")
    # out = {"membership": membership, "info": info}
    membership = out["membership"]
    if membership.shape[0] == 0:
        return out

    membership, res, excluded = _remove_short_tips(links, excluded, membership, info,
                                                   client=client, save_dir=save_dir,
                                                   min_dist=min_dist, ncores=ncores)
    if membership.shape[0] == 0:
        print("This set of parameters leads to lose everything.")
        return membership, info, excluded
    if res is not None:
        out["info"] = res

    membership, res, excluded = _remove_bifurcations(links, excluded, membership, info,
                                                     client=client, save_dir=save_dir,
                                                     min_dist=min_dist, ncores=ncores)
    if membership.shape[0] == 0:
        print("This set of parameters leads to lose everything.")
        return membership, info, excluded
    if res is not None:
        out["info"] = res

    degree = _calculate_degree(links, excluded)
    add = degree.merge(
        membership.query("rank == 1"), on="scaffold_index").query("degree == 1")
    add = add.index.values.compute()
    if add.shape[0] > 0:
        excluded.update(add.tolist())
        out = make_super_scaffolds(links=links, info=info, excluded=excluded, ncores=ncores,
                                   client=client, save_dir=save_dir)
        membership = dd.read_parquet(out["membership"], infer_divisions=True)

    add = membership.query("rank > 0")["scaffold_index"]
    if add.shape[0] > 0:
        excluded.update(add.tolist())
        out = make_super_scaffolds(links=links, info=info, excluded=excluded, ncores=ncores,
                                   client=client, save_dir=save_dir)

    return membership, out["info"], excluded
