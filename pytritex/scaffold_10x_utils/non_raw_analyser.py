import pandas as pd
import numpy as np
from ..make_super_scaffolds import make_super_scaffolds
from .scaffold_unanchored import _scaffold_unanchored

# if(!raw){
#   if(verbose){
#    cat("Tip removal.\n")
#   }
#   # remove short tips of rank 1
#   links[!scaffold2 %in% ex][, .(degree=.N), key=.(scaffold=scaffold1)]->b
#   a <- b[
#     m[m[rank == 1][, .(super=c(super, super, super), bin=c(bin, bin-1, bin+1))], on=c("super", "bin")],
#          on="scaffold"]
#   a[degree == 1 & length <= 1e4]$scaffold->add
#   if(length(add) > 0){
#    ex <- c(ex, add)
#    make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
#    out$m -> m
#   }
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
#   b[m[rank == 1], on="scaffold"][degree == 1]$scaffold -> add
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
#    m[, .(scaffold1=scaffold, super1=super, chr1=chr, cM1=cM, size1=super_nbin, d1 = pmin(bin - 1, super_nbin - bin))][xy, on="scaffold1"]->xy
#    m[, .(scaffold2=scaffold, super2=super, chr2=chr, cM2=cM, size2=super_nbin, d2 = pmin(bin - 1, super_nbin - bin))][xy, on="scaffold2"]->xy
#    xy[super2 != super1 & d1 == 0 & d2 == 0 & size1 > 1 & size2 > 1 & chr1 == chr2]->xy
#    xy[scaffold1 < scaffold2, .(nscl=.N), scaffold_link][xy, on="scaffold_link"]->xy
#    xy[nscl == 1] -> xy
#    xy[super1 < super2][, c("n", "g"):=list(.N, .GRP), by=.(super1, super2)][order(-link_length)][!duplicated(g)]->zz
#
#    sel <- zz[, .(scaffold1=c(scaffold_link, scaffold_link, scaffold1, scaffold2),
# 	  scaffold2=c(scaffold1, scaffold2, scaffold_link, scaffold_link))]
#    rbind(links, ww2[sel, on=c("scaffold1", "scaffold2")])->links2
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


def non_raw_analyser(links: pd.DataFrame, excluded, membership: pd.DataFrame, info: pd.DataFrame,
                     prefix=None, ncores=1, verbose=True, unanchored=True, popseq_dist=1e4):

    """This part of scaffold10x is for when we are not in "raw" mode"""

    if verbose:
        print("Tip removal")

    # remove short tips of rank 1
    #   links[!scaffold2 %in% ex][, .(degree=.N), key=.(scaffold=scaffold1)]->b
    b = links.loc[~links["scaffold2"].isin(excluded)].rename(
        columns={"scaffold1": "scaffold"}).groupby(["scaffold"]).size().to_frame("degree")
    bait = membership["rank"] == 1
    a = b.merge(
        membership.merge(membership,
                         pd.DataFrame().assign(super=pd.concat([membership.loc[bait, "super"]]),
                              bin=pd.concat([membership.loc[bait]["bin"],
                                             membership.loc[bait]["bin"] - 1,
                                             membership.loc[bait]["bin"] + 1])), on=["super", "bin"]),
        left_index=True, right_on=["scaffold"], how="right")
    # TODO: why the hard limit on the length?
    add = a.loc[(a["degree"] == 1) & (a["length"] <= 1e4)]["scaffold"]
    if add.shape[0] > 0:
        excluded = excluded.extend(add)
        out = make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=excluded, ncores=ncores)
        membership = out["membership"]
    # Now do the same but on the new membership
    add = membership.loc[(membership["rank"] == 1) & (membership["length"] <= 1e4), "scaffold"]
    if add.shape[0] > 0:
        excluded = excluded.extend(add)
        out = make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=excluded, ncores=ncores)
        membership = out["membership"]

    # resolve length-one-bifurcations at the ends of paths
    x = membership.loc[(membership["rank"] > 0) &
                       ((membership["bin"] == 2) | (membership["super_nbin"] - 1 == membership["bin"])),
    ["super", "super_nbin", "bin", "scaffold", "length"]].assign(type=lambda df: df["bin"] == 2,
                                                                 bin0=lambda df: df["bin"]).loc[:,
        ["super", "super_nbin", "type", "scaffold", "length", "bin0"]]
    a = pd.concat(
        [
         membership.merge(x.loc[~x["type"], ["super", "bin0"]].assign(bin=1), on=["super", "bin"], how="right"),
         membership.merge(x.loc[x["type"], ["super", "bin0", "super_nbin"]].rename(
             columns={"super_nbin": "bin"}), on=["super", "bin"], how="right")
        ]
    )

    __left = a.rename(columns={"scaffold": "scaffold2", "length": "length2"})[
        ["super", "bin0", "scaffold2", "length2"]]
    add = __left.merge(x, on=["super", "bin0"], how="right").apply(
        lambda row: row.scaffold2 if row.length >= row.length2 else row.scaffold)
    if add.shape[0] > 0:
        excluded = excluded.extend(add)
        out = make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=excluded, ncores=ncores)
        membership = out["membership"]

    # Now re-do the same but on the new membership
    add = membership.loc[(membership["rank"] == 1) & (membership["length"] <= 1e4), "scaffold"]
    if add.shape[0] > 0:
        excluded = excluded.extend(add)
        out = make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=excluded, ncores=ncores)
        membership = out["membership"]

    # remove tips of rank 1
    b = links[~links["scaffold2"].isin(excluded)].rename(columns={"scaffold1": "scaffold"}).groupby(
        ["scaffold"]).size().to_frame("degree")
    add = b.merge(membership.loc[membership["rank"] == 1, :], left_index=True, right_on=["scaffold"], how="inner").loc[
        lambda df: df["degree"] == 1, "scaffold"]
    if add.shape[0] > 0:
        excluded = excluded.extend(add.to_list())
        out = make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=excluded, ncores=ncores)
        membership = out["membership"]

    # remove remaining nodes of rank > 0
    add = membership.loc[membership["rank"] > 0, "scaffold"]
    if add.shape[0] > 0:
        excluded = excluded.extend(add.to_list())
        out = make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=excluded, ncores=ncores)
        membership = out["membership"]

    if popseq_dist > 0 and unanchored is True:
        out, excluded = _scaffold_unanchored()

    membership, res = out["membership"], out["info"]
    return membership, res, excluded
