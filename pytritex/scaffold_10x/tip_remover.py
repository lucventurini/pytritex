import pandas as pd
from pytritex.graph_utils.make_super_scaffolds import make_super_scaffolds
from functools import partial
from .scaffold_unanchored import _scaffold_unanchored
import numpy as np


def _relaunch_scaffolding(f):
    """Wrapper to re-run the super-scaffolding algorithm every time we have to add excluded scaffolds."""
    def _wrapped(*args, **kwargs):
        membership, excluded, add = f(*args, **kwargs)
        if add.shape[0] > 0:
            excluded.update(add.tolist())
            out = make_super_scaffolds(links=kwargs["links"],
                                       prefix=kwargs["prefix"],
                                       info=kwargs["info"],
                                       excluded=excluded,
                                       ncores=kwargs["ncores"])
            membership = out["membership"]
        return membership, excluded
    return _wrapped


def _remove_bulges(f):

# if(!raw){
#   if(verbose){
#    cat("Tip removal.\n")
#   }
#   # remove short tips of rank 1
#   links[!scaffold2 %in% ex][, .(degree=.N), key=.(scaffold=scaffold1)]->b
#   a <- b[
#     m[m[rank == 1][, .(super=c(super, super, super), bin=c(bin, bin-1, bin+1))], on=c("super", "bin")],
#          on="scaffold_index"]
#   a[degree == 1 & length <= 1e4]$scaffold->add
#   if(length(add) > 0){
#    ex <- c(ex, add)
#    make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
#    out$m -> m
#   }


@_relaunch_scaffolding
def _remove_short_tips(links: pd.DataFrame, excluded, membership: pd.DataFrame, info: pd.DataFrame,
                       prefix=None, ncores=1, popseq_dist=1e4):
    """First operation: remove scaffolds where the rank is 1 and the degree (=> no. links) is 1."""
    degree = links.query("scaffold_index1 not in @excluded").groupby("scaffold_index1").size().to_frame("degree")
    degree.index.names = ["scaffold_index"]
    inner = membership.query("rank == 1")[["super", "bin"]].values
    inner = pd.DataFrame().assign(
        super=np.tile(inner[:, 0], 3),
        bin=np.tile(inner[:, 1], 3).reshape(inner[:, 1].shape[0], 3) + np.array([0, -1, 1]).reshape(3, 1))
    a = membership.merge(inner, on=["super", "bin"], how="right")
    a = degree.merge(a, on="scaffold_index").query("(degree == 1) & (length <= @popseq_dist)")
    add = a["scaffold_index"].values
    return membership, excluded, add


@_relaunch_scaffolding
def _remove_bulges(links: pd.DataFrame, excluded, membership: pd.DataFrame, info: pd.DataFrame,
                   prefix=None, ncores=1, popseq_dist=1e4):
    add = membership.query("rank == 1 & length <= @popseq_dist")
#
#   # remove short tips/bulges of rank 1
#   m[rank == 1 & length <= 1e4]$scaffold -> add
#   ex <- c(ex, add)
#   if(length(add) > 0){
#    ex <- c(ex, add)
#    make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
#    out$m -> m
#   }
#
#   # resolve length-one-bifurcations at the ends of paths
#   m[rank > 0][bin == 2 | super_nbin - 1 == bin ][, .(super, super_nbin, type = bin == 2, scaffold, length, bin0=bin)]->x
#   unique(rbind(
#   m[x[type == T, .(super, bin0, bin=1)], on=c("super", "bin")],
#   m[x[type == F, .(super, bin0, bin=super_nbin)], on=c("super", "bin")]
#   ))->a
#   a[, .(super, bin0, scaffold2=scaffold, length2=length)][x, on=c("super", "bin0")][, ex := ifelse(length >= length2, scaffold2, scaffold)]$ex -> add
#
#   if(length(add) > 0){
#    ex <- c(ex, add)
#    make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
#    out$membership -> m
#   }
#
#   # remove short tips/bulges of rank 1
#   m[rank == 1 & length <= 1e4]$scaffold -> add
#   if(length(add) > 0){
#    ex <- c(ex, add)
#    make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
#    out$m -> m
#   }
#
#   # remove tips of rank 1
#   links[!scaffold2 %in% ex][, .(degree=.N), key=.(scaffold=scaffold1)]->b
#   b[m[rank == 1], on="scaffold_index"][degree == 1]$scaffold -> add
#   if(length(add) > 0){
#    ex <- c(ex, add)
#    make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
#    out$m -> m
#   }
#
#   # remove remaining nodes of rank > 0
#   m[rank > 0]$scaffold -> add
#   if(length(add) > 0){
#    ex <- c(ex, add)
#    make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
#    out$m -> m
#   }
#
#   if(popseq_dist > 0 & unanchored == T){
#    if(verbose){
#     cat("Including unanchored scaffolds.\n")
#    }
#    # use unanchored scaffolds to link super-scaffolds
#    ww2[is.na(popseq_chr1), .(scaffold_link=scaffold1, link_length=length1, scaffold1=scaffold2)]->x
#    ww2[is.na(popseq_chr1), .(scaffold_link=scaffold1, scaffold2=scaffold2)]->y
#    x[y, on="scaffold_link", allow.cartesian=T][scaffold1 != scaffold2]->xy
#
#    m[, .(scaffold1=scaffold, super1=super, chr1=chr, cM1=cM, size1=super_nbin, d1 = pmin(bin - 1, super_nbin - bin))][xy, on="scaffold_index1"]->xy
#    m[, .(scaffold2=scaffold, super2=super, chr2=chr, cM2=cM, size2=super_nbin, d2 = pmin(bin - 1, super_nbin - bin))][xy, on="scaffold_index2"]->xy
#    xy[super2 != super1 & d1 == 0 & d2 == 0 & size1 > 1 & size2 > 1 & chr1 == chr2]->xy
#    xy[scaffold1 < scaffold2, .(nscl=.N), scaffold_link][xy, on="scaffold_link"]->xy
#    xy[nscl == 1] -> xy
#    xy[super1 < super2][, c("n", "g"):=list(.N, .GRP), by=.(super1, super2)][order(-link_length)][!duplicated(g)]->zz
#
#    sel <- zz[, .(scaffold1=c(scaffold_link, scaffold_link, scaffold1, scaffold2),
# 	  scaffold2=c(scaffold1, scaffold2, scaffold_link, scaffold_link))]
#    rbind(links, ww2[sel, on=c("scaffold_index1", "scaffold_index2")])->links2
#
#    make_super_scaffolds(links=links2, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
#    out$m -> m
#
#    #resolve branches
#    m[rank > 0][bin == 2 | super_nbin - 1 == bin ][, .(super, super_nbin, type = bin == 2, scaffold, length, bin0=bin)]->x
#    unique(rbind(
#    m[x[type == T, .(super, bin0, bin=1)], on=c("super", "bin")],
#    m[x[type == F, .(super, bin0, bin=super_nbin)], on=c("super", "bin")]
#    ))->a
#    a[, .(super, bin0, scaffold2=scaffold, length2=length)][x, on=c("super", "bin0")][, ex := ifelse(length >= length2, scaffold2, scaffold)]$ex -> add
#
#    if(length(add) > 0){
#     ex <- c(ex, add)
#     make_super_scaffolds(links=links2, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
#     out$membership -> m
#    }
#
#    # remove remaining nodes of rank > 0
#    m[rank > 0]$scaffold -> add
#    if(length(add) > 0){
#     ex <- c(ex, add)
#     make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
#    }
#   }
#
#   out$m -> m
#   out$info -> res
#  }


@_relaunch_scaffolding
def _remove_bifurcations(links, excluded, membership, info, prefix=None, ncores=1):
    keys = ["super", "super_nbin", "bin", "scaffold_index", "length"]
    new_keys = ["super", "super_nbin", "type", "scaffold_index", "length", "bin0"]
    x = membership.query("(rank > 0) & ((bin == 2) | (super_nbin - 1 == bin))", inplace=False)[keys]
    x = x.eval("type = (bin == 2)").eval("bin0 = bin")[new_keys]
    upper = membership.merge(x.loc[x["type"], ["super", "bin0"]].eval("bin=1"), on=["super", "bin"], how="right")
    lower = x.loc[~x["type"], ["super", "bin0", "super_nbin"]]
    lower.columns = ["super", "bin0", "bin"]
    lower = membership.merge(lower, on=["super", "bin"], how="right")
    a = pd.concat([upper, lower]).reset_index(drop=True)[["super", "bin0", "scaffold_index", "length"]]
    a.columns = ["super", "bin0", "scaffold_index2", "length2"]
    a = a.merge(x, on=["super", "bin0"]).eval("excluded = scaffold_index2 if length >= length2 else scaffold_index")
    add = a["excluded"].values
    return membership, excluded, add


def non_raw_analyser(links: pd.DataFrame, excluded, membership: pd.DataFrame, info: pd.DataFrame,
                     prefix=None, ncores=1, verbose=True, unanchored=True, popseq_dist=1e4):

    """This part of scaffold10x is for when we are not in "raw" mode"""

    if verbose:
        print("Tip removal")
    membership = _remove_short_tips(links, excluded, membership, info=info, prefix=prefix,
                                    ncores=ncores, popseq_dist=popseq_dist)


    return
