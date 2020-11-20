import pandas as pd
import dask.dataframe as dd
import numpy as np


def bin_hic_step(hic, frags, binsize, chrlen: pd.DataFrame, step=None, chrs=None, cores=1):
    if chrs is None:
        chrs = np.unique(chrlen["chr"].values)

    if step is None:
        step = binsize

    bins = chrlen.query("chr in @chrs", local_dict={"chrs": chrs}).groupby("chr")


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
#  chrlen[chr %in% chrs, .(bin=seq(0, length, step)), key=chr][, end := bin + binsize - 1][, start := bin] -> bins
#  bins[, id := paste0(chr, ":", bin)]
#
#  hic[chr1 %in% chrs & chr2 %in% chrs] -> z
#
#  if(step == binsize){
#   z[, bin1 := as.integer(start1 %/% binsize * binsize)]
#   z[, bin2 := as.integer(start2 %/% binsize * binsize)]
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