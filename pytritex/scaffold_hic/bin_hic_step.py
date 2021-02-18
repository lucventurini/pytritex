from functools import partial

import pandas as pd
import dask.dataframe as dd
import dask.array as da
import numpy as np
from typing import Dict


def weighted_avg(group, by, weights):
    """Weighted average for groups.

    See https://stackoverflow.com/a/33392117/5188860"""
    d = group[by]
    w = group[weights]

    return (d * w).sum() / w.sum()


def _same_step(chr_hic, frags, binsize, chrs):
    chr_hic["bin1"] = chr_hic["start1"] // binsize * binsize
    chr_hic["bin2"] = chr_hic["start2"] // binsize * binsize
    chr_hic["size"] = 1
    chr_hic["nlinks"] = chr_hic.groupby(["chr1", "chr2", "bin1", "bin2"])["size"].transform(sum,
                                                                                            meta=int).compute()
    mat = chr_hic[["chr1", "chr2", "bin1", "bin2", "nlinks"]].drop_duplicates()
    mat["idx"] = da.from_array(np.arange(1, mat.shape[0].compute() + 1, dtype=int))
    mat = mat.set_index("idx")

    chr_frags = frags.query("chr in @chrs", local_dict={"chrs": chrs})
    chr_frags["bin"] = chr_frags["start"] // binsize * binsize
    chr_frags["bin"] = chr_frags["bin"].astype(int)
    grouped = chr_frags.groupby(["chr", "bin"])
    nfrags = grouped.size().rename("nfrags")
    wavg = partial(weighted_avg, by="cov", weights="length")
    cov = grouped.apply(wavg, meta=float).rename("cov")
    gc = grouped.apply(lambda group: (group.nC.sum() + group.nG.sum()) / (
            group.nA.sum() + group.nC.sum() + group.nG.sum() + group.nT.sum()), meta=float).rename("gc")
    chr_frags = dd.merge(cov, gc, on=["chr", "bin"]).merge(nfrags, on=["chr", "bin"])
    chr_frags["frag_index"] = 1
    chr_frags["frag_index"] = chr_frags["frag_index"].cumsum()
    chr_frags = chr_frags.set_index("frag_index")
    return mat, chr_frags


def _different_step(chrlen, chr_hic, frags, binsize, step, chrs):
    bins = chrlen.query("chr in @chrs", local_dict={"chrs": chrs})
    bins["bin"] = bins["length"].apply(lambda x: range(0, x, step))
    bins = bins.explode("bin")[["chr", "bin"]].drop_duplicates()
    bins["win"] = bins.apply(lambda row: np.arange(row["bin"], row["bin"] + binsize - 1, 25, dtype=int), axis=1)
    bins = bins.explode("win")

    # bins["start"] = bins["bin"]
    # bins["end"] = bins["bin"] + binsize - 1
    chr_hic = chr_hic.query("chr1 == chr2")
    chr_hic["win1"] = (chr_hic["start1"] // step * step).astype(int)
    chr_hic["win2"] = (chr_hic["start2"] // step * step).astype(int)
    mat = dd.merge(chr_hic, bins.rename(columns={"chr": "chr1", "win": "win1"}), on=["chr1", "win1"]).merge(
        bins.rename(columns={"chr": "chr1", "win": "win1"}), on=["chr2", "win2"])
    mat = mat.groupby(["chr1", "bin1", "chr2", "bin2"]).size().rename("nlinks").to_frame().reset_index(drop=False)
    mat["idx"] = 1
    mat["idx"] = mat["idx"].cumsum()
    mat = mat.set_index("idx")
    mat["dist"] = mat.eval("bin2 - bin1").abs()

    chr_frags = frags.query("chr in @chrs", local_dict={"chrs": chrs})
    chr_frags["win"] = (chr_frags["start"] // step * step).astype(int)
    dd.merge(bins, chr_frags, on=["chr", "win"])
    grouped = chr_frags.groupby(["chr", "bin"])
    nfrags = grouped.size().rename("nfrags")
    wavg = partial(weighted_avg, by="cov", weights="length")
    cov = grouped.apply(wavg, meta=float).rename("cov")
    gc = grouped.apply(lambda group: (group.nC.sum() + group.nG.sum()) / (
            group.nA.sum() + group.nC.sum() + group.nG.sum() + group.nT.sum()), meta=float).rename("gc")
    chr_frags = dd.merge(cov, gc, on=["chr", "bin"]).merge(nfrags, on=["chr", "bin"])
    return mat, chr_frags


def bin_hic_step(hic: dd.DataFrame, frags: dd.DataFrame,
                 binsize: int, chrlen: pd.DataFrame, step=None, chrs=None, cores=1) -> Dict[(str, dd.DataFrame)]:
    """Bin HiC links in windows of fixed size"""

    if chrs is None:
        chrs = np.unique(chrlen["chr"].values)

    if step is None:
        step = binsize

    chr_hic = hic.query("chr1 in @chrs & chr2 in @chrs", local_dict={"chrs": chrs})
    if step == binsize:
        mat, chr_frags = _same_step(chr_hic, frags, binsize, chrs)
    else:
        mat, chr_frags = _different_step(chr_hic, frags, binsize=binsize, step=step, chrs=chrs)
    return {"bins": chr_frags, "mat": mat}


# # Bin Hi-C links in windows of a fixed size
# bin_hic_step<-function(hic, frags, binsize, step=NULL, chrlen, chrs=NULL, cores=1){
#  if(is.null(chrs)){
#   chrs <- unique(chrlen$chr)
#  }
#
#  if(is.null(step)){
#   step <- binsize
#  }
#
#  options(scipen=20)
#  bins <- chrlen[chr %in% chrs, .(bin=seq(0, length, step)), key=chr]
#  bins[, end := bin + binsize - 1][, start := bin]
#  bins[, id := paste0(chr, ":", bin)]
#
#  hic[chr1 %in% chrs & chr2 %in% chrs] -> z
#
#  if(step == binsize){

#   z[, .(nlinks=.N), keyby=.(chr1,bin1,chr2,bin2)]->z
#   z[, id1 := paste(sep=":", chr1, bin1)]
#   z[, id2 := paste(sep=":", chr2, bin2)]
#   z->mat
#
#   frags[chr %in% chrs]->f
#   f[, bin := as.integer(start %/% binsize * binsize)]
#   f[, .(nfrags=.N, eff_length=sum(length), cov=weighted.mean(cov, length),
# 	gc=(sum(nC)+sum(nG))/(sum(nA)+sum(nC)+sum(nG)+sum(nT))), keyby=.(chr, bin)]->f
#   f[, id := paste(sep=":", chr, bin)]
#  } else {
#   copy(bins)->b
#   b[, idx := 1:.N]
#   b[, .(win=seq(bin, bin+binsize-1, step)),.(chr, bin)]->b
#
#   z[chr1 == chr2]->z
#   z[, win1:=as.integer(start1 %/% step * step)]
#   z[, win2:=as.integer(start2 %/% step * step)]
#
#   rbindlist(mclapply(mc.cores=cores, mc.preschedule=F, chrs, function(i){
#    z[chr1 == i]->y
#    y[, .(nlinks=.N), keyby=.(chr1,win1,chr2,win2)]->y
#    setnames(copy(b), paste0(names(b), 1))[y, on=c("chr1", "win1"), allow.cartesian=T]->y
#    setnames(copy(b), paste0(names(b), 2))[y, on=c("chr2", "win2"), allow.cartesian=T]->y
#    y[, .(nlinks=.N), keyby=.(chr1,bin1,chr2,bin2)]
#   }))->y
#   y[, id1:=stri_c(sep=":", chr1, bin1)]
#   y[, id2:=stri_c(sep=":", chr2, bin2)]
#   y[, dist := abs(bin1 - bin2)]
#   y->mat
#
#   copy(frags)->f
#   f[, win := as.integer(start %/% step * step)]
#   b[f, on=c("chr", "win"), nomatch=0, allow.cartesian=T]->f
#   f[, .(nfrags=.N, eff_length=sum(length), cov=weighted.mean(cov, length),
# 	gc=(sum(nC)+sum(nG))/(sum(nA)+sum(nC)+sum(nG)+sum(nT))), keyby=.(chr, bin)]->f
#   f[, id :=paste(sep=":", chr, bin)]
#  }
#
#  list(bins=f[], mat=mat[])
# }