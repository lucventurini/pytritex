import pandas as pd
import numpy as np
from ..make_super_scaffolds import make_super_scaffolds


# if(verbose){
#         cat("Including unanchored scaffolds.\n")
#        }
#        # use unanchored scaffolds to link super-scaffolds
#        ww2[is.na(popseq_chr1), .(scaffold_link=scaffold1, link_length=length1, scaffold1=scaffold2)]->x
#        ww2[is.na(popseq_chr1), .(scaffold_link=scaffold1, scaffold2=scaffold2)]->y
#        x[y, on="scaffold_link", allow.cartesian=T][scaffold1 != scaffold2]->xy
#
#        m[, .(scaffold1=scaffold, super1=super, chr1=chr, cM1=cM, size1=super_nbin, d1 = pmin(bin - 1, super_nbin - bin))][xy, on="scaffold1"]->xy
#        m[, .(scaffold2=scaffold, super2=super, chr2=chr, cM2=cM, size2=super_nbin, d2 = pmin(bin - 1, super_nbin - bin))][xy, on="scaffold2"]->xy
#        xy[super2 != super1 & d1 == 0 & d2 == 0 & size1 > 1 & size2 > 1 & chr1 == chr2]->xy
#        xy[scaffold1 < scaffold2, .(nscl=.N), scaffold_link][xy, on="scaffold_link"]->xy
#        xy[nscl == 1] -> xy
#        xy[super1 < super2][, c("n", "g"):=list(.N, .GRP), by=.(super1, super2)][order(-link_length)][!duplicated(g)]->zz
#
#        sel <- zz[, .(scaffold1=c(scaffold_link, scaffold_link, scaffold1, scaffold2),
#     	  scaffold2=c(scaffold1, scaffold2, scaffold_link, scaffold_link))]
#        rbind(links, ww2[sel, on=c("scaffold1", "scaffold2")])->links2
#
#        make_super_scaffolds(links=links2, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
#        out$m -> m
#
#        #resolve branches
#        m[rank > 0][bin == 2 | super_nbin - 1 == bin ][, .(super, super_nbin, type = bin == 2, scaffold, length, bin0=bin)]->x
#        unique(rbind(
#        m[x[type == T, .(super, bin0, bin=1)], on=c("super", "bin")],
#        m[x[type == F, .(super, bin0, bin=super_nbin)], on=c("super", "bin")]
#        ))->a
#        a[, .(super, bin0, scaffold2=scaffold, length2=length)][x, on=c("super", "bin0")][, ex := ifelse(length >= length2, scaffold2, scaffold)]$ex -> add
#
#        if(length(add) > 0){
#         ex <- c(ex, add)
#         make_super_scaffolds(links=links2, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
#         out$membership -> m
#        }
#
#        # remove remaining nodes of rank > 0
#        m[rank > 0]$scaffold -> add
#        if(length(add) > 0){
#         ex <- c(ex, add)
#         make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
#        }


def _scaffold_unanchored(links, excluded, membership, info, ww2, ncores=1,
                         prefix=None, verbose=False):

    if verbose:
        print("Including unanchored scaffolds")

    # use unanchored scaffolds to link super-scaffolds
    x = ww2.loc[(ww2["popseq_chr1"].isna()), ["scaffold1", "length1", "scaffold2"]].rename(
        columns={"scaffold1": "scaffold_link", "length1": "link_length", "scaffold2": "scaffold1"}
    )
    y = ww2.loc[(ww2["popseq_chr1"].isna()), ["scaffold1", "scaffold2"]].rename(columns={"scaffold1": "scaffold_link"})
    xy = x.merge(y, on=["scaffold_link"], how="outer").loc[lambda df: df["scaffold1"] != df["scaffold2"]]

    # m[, .(scaffold1=scaffold, super1=super, chr1=chr, cM1=cM, size1=super_nbin, d1 = pmin(bin - 1, super_nbin - bin))
    # ][xy, on="scaffold1"]->xy

    xy = membership.assign(
        d1=np.minimum(membership["bin"] - 1, membership["super_nbin"] - membership["bin"])).rename(
        columns={"scaffold": "scaffold1", "super": "super1", "chr": "chr1", "cM": "cM1", "super_nbin": "size1"}
    )[["scaffold1", "super1", "chr1", "cM1", "size1", "d1"]].merge(xy, on="scaffold1", how="right")

    # m[, .(scaffold2=scaffold, super2=super, chr2=chr, cM2=cM, size2=super_nbin, d2 = pmin(bin - 1, super_nbin - bin))
    # ][xy, on="scaffold2"]->xy

    xy = membership.assign(
        d2=np.minimum(membership["bin"] - 1, membership["super_nbin"] - membership["bin"])).rename(
        columns={"scaffold": "scaffold2", "super": "super2", "chr": "chr2", "cM": "cM2", "super_nbin": "size2"}
    )[["scaffold2", "super2", "chr2", "cM2", "size2", "d2"]].merge(xy, on="scaffold2", how="right")

    xy = xy.loc[xy.eval("super2 != super1 & d1 == 0 & d2 == 0 & size1 > 1 & size2 > 1 & chr1 == chr2")]
    xy = xy.loc[xy.eval("scaffold1 < scaffold2")].groupby("scaffold_link").size().to_frame("nscl").merge(
        xy, how="right", left_index=True, right_on="scaffold_link")
    xy = xy.loc[xy["nscl"] == 1]
    n_g = xy.loc[xy.eval("super1 < super2")].groupby(
        ["super1", "super2"]).size().to_frame("n").assign(g=lambda df: np.arange(df.shape[0], dtype=int))
    zz = n_g.merge(xy, left_index=True, right_on=n_g.index.names).sort_values("link_length", ascending=False).loc[
        lambda df: ~df["g"].duplicated()]

    sel = pd.DataFrame().assign(
        scaffold1=pd.concat([zz["scaffold_link"], zz["scaffold_link"], zz["scaffold1"], zz["scaffold2"]]),
        scaffold2=pd.concat([zz["scaffold1"], zz["scaffold2"], zz["scaffold_link"], zz["scaffold_link"]])
    )
    links2 = pd.concat(links, ww2.merge(sel, on=["scaffold1", "scaffold2"], how="right"))
    out = make_super_scaffolds(links=links2, prefix=prefix, info=info, excluded=excluded, ncores=ncores)
    membership = out["membership"]
    #resolve branches
    x = membership.loc[(membership["rank"] > 0) &
                       ((membership["bin"] == 2) | (membership["super_nbin"] - 1 == membership["bin"])),
        ["super", "super_nbin", "bin", "scaffold", "length"]].assign(type=lambda df: df["bin"] == 2).rename(
        columns={"bin": "bin0"})[["super", "super_nbin", "type", "scaffold", "length", "bin0"]]
    a = pd.concat([
        membership.merge(x.loc[x["type"], ["super", "bin0"]].assign(bin=1), on=["super", "bin"], how="right"),
        membership.merge(x.loc[~x["type"], ["super", "bin0", "super_nbin"]].rename(columns={
            "super_nbin": "bin"}), on=["super", "bin"], how="right")
    ])

    add = a.rename(columns={"scaffold": "scaffold2", "length": "length2"})[
        ["super", "bin0", "scaffold2", "length2"]].merge(x, on=["super", "bin0"]).apply(
        lambda row: row.scaffold2 if row.length >= row.length2 else row.scaffold, axis=1)

    if add.shape[0] > 0:
        excluded = excluded.extend(add.to_list())
        out = make_super_scaffolds(links=links2, prefix=prefix, info=info, excluded=excluded, ncores=ncores)
        membership = out["membership"]

    add = membership.loc[membership["rank"] > 0, "scaffold"]
    if add.shape[0] > 0:
        excluded = excluded.extend(add.to_list())
        # TODO: why are we going back to links instead of using again links2 here?
        out = make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=excluded, ncores=ncores)

    return out, excluded
