import pandas as pd
from pytritex.graph_utils.make_super_scaffolds import make_super_scaffolds
from functools import partial
from .scaffold_unanchored import _scaffold_unanchored
import numpy as np
from time import ctime


minimum_distance = 1e4


def _calculate_degree(links, excluded):
    degree = links.query("scaffold_index2 not in @excluded").groupby("scaffold_index1").size().to_frame("degree")
    degree.index.names = ["scaffold_index"]
    return degree


def _remove_bulges(links, excluded, membership, info, min_dist=minimum_distance, ncores=1, prefix=None):
    add = membership.query("(rank == 1) & (length <= @min_dist)")["scaffold_index"]
    if add.shape[0] > 0:
        excluded.update(add.tolist())
        out = make_super_scaffolds(links=links, info=info, excluded=excluded, ncores=ncores, prefix=prefix)
        membership = out["membership"]
        res = out["info"]
    else:
        res = None
    return membership, res, excluded


def _remove_short_tips(links, excluded, membership, info, min_dist=minimum_distance, ncores=1, prefix=None):
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

    inner = membership.query("rank == 1")[["super", "bin"]].values
    inner0 = np.tile(inner[:, 0], 3)
    inner1 = np.tile(inner[:, 1], 3).reshape(3, inner.shape[0]) + np.array([0, -1, 1]).reshape(3, 1)
    inner1 = inner1.flatten()
    inner = pd.DataFrame().assign(super=inner0, bin=inner1)
    middle = membership.merge(inner, how="right", on=["super", "bin"]).set_index("scaffold_index")
    degree = _calculate_degree(links, excluded)
    a = degree.merge(middle, on="scaffold_index").reset_index(drop=False)
    add = a.query("(degree == 1) & (length <= @min_dist)")["scaffold_index"].values
    out = {"membership": membership, "info": None}
    if add.shape[0] > 0:
        excluded.update(add.tolist())
        out = make_super_scaffolds(links=links, info=info, excluded=excluded, ncores=ncores, prefix=prefix)
        membership = out["membership"]

    membership, res, excluded = _remove_bulges(links=links, excluded=excluded,
                                               membership=membership,
                                               info=info, min_dist=min_dist, ncores=ncores, prefix=prefix)
    if res is not None:
        out["info"] = res

    return membership, out["info"], excluded


def _remove_bifurcations(links, excluded, membership, info, min_dist=minimum_distance, ncores=1, prefix=None):
    # resolve length-one-bifurcations at the ends of paths
    # m[rank > 0][bin == 2 | super_nbin - 1 == bin][,.(super, super_nbin, type = bin == 2, scaffold, length, bin0=bin)]->x
    # upper <- m[x[type == T,.(super, bin0, bin=1)], on=c("super", "bin")]
    # lower <- m[x[type == F,.(super, bin0, bin=super_nbin)], on = c("super", "bin")]
    # a <- unique(rbind(upper, lower))
    # a[,.(super, bin0, scaffold2=scaffold, length2=length)][x, on = c("super", "bin0")][,
    #     ex := ifelse(length >= length2, scaffold2, scaffold)]$ex -> add
    #
    # if (length(add) > 0){
    # ex < - c(ex, add)
    # make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
    # out$membership -> m
    # }
    #
    # # remove short tips/bulges of rank 1
    # m[rank == 1 & length <= 1e4]$scaffold -> add
    # if (length(add) > 0){
    # ex < - c(ex, add)
    # make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
    # out$m -> m
    # }
    keys = ["super", "super_nbin", "scaffold_index", "length", "bin"]
    x = membership.query("(rank > 0) & ( (bin == 2) | (super_nbin - 1 == bin) )")[keys]
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
        columns={"scaffold_index": "scaffold_index2", "length": "length2"}).set_index(["super", "bin0"])
    # Now merge with "x" to find places to exclude
    a = a.merge(x.set_index(["super", "bin0"]), on=["super", "bin0"], how="right")
    add = np.where((a.length >= a.length2), a.scaffold_index2, a.scaffold_index)
    out = {"membership": membership, "info": None}
    if add.shape[0] > 0:
        excluded.update(add.tolist())
        out = make_super_scaffolds(links=links, info=info, excluded=excluded, ncores=ncores, prefix=prefix)
        membership = out["membership"]
    membership, res, excluded = _remove_bulges(links=links, excluded=excluded,
                                               membership=membership,
                                               info=info, min_dist=min_dist, ncores=ncores, prefix=prefix)
    if res is not None:
        out["info"] = res
    return membership, out["info"], excluded


def remove_tips(links, excluded, out, info, ncores=1, prefix=None,
                verbose=False, min_dist=minimum_distance):

    if verbose:
        print(ctime(), "Starting tip removal")
    # out = {"membership": membership, "info": info}
    membership = out["membership"]
    membership, res, excluded = _remove_short_tips(links, excluded, membership, info,
                                                   min_dist=min_dist, ncores=ncores, prefix=prefix)
    if res is not None:
        out["info"] = res

    membership, res, excluded = _remove_bifurcations(links, excluded, membership, info,
                                                     min_dist=min_dist, ncores=ncores, prefix=prefix)
    if res is not None:
        out["info"] = res

    # remove tips of rank 1
    # links[!scaffold2 % in % ex][,.(degree=.N), key =.(scaffold=scaffold1)]->b
    # b[m[rank == 1], on = "scaffold"][degree == 1]$scaffold -> add
    # if (length(add) > 0){
    # ex < - c(ex, add)
    # make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
    # out$m -> m
    # }
    degree = _calculate_degree(links, excluded)
    add = degree.merge(membership.query("rank == 1"), on="scaffold_index").query("degree == 1")["scaffold_index"]
    if add.shape[0] > 0:
        excluded.update(add.tolist())
        out = make_super_scaffolds(links=links, info=info, excluded=excluded, ncores=ncores, prefix=prefix)
        membership = out["membership"]

    add = membership.query("rank > 0")["scaffold_index"]
    if add.shape[0] > 0:
        excluded.update(add.tolist())
        out = make_super_scaffolds(links=links, info=info, excluded=excluded, ncores=ncores, prefix=prefix)

    return membership, out["info"], excluded
