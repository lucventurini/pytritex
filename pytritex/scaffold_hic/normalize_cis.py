import dask.dataframe as dd
import pandas as pd
import numpy as np


# Intra-chromosomal normalization, taken from HiCNorm https://academic.oup.com/bioinformatics/article/28/23/3131/192582
# normalize_mat_cis<-function(u, v){
#  u_vec<-u[upper.tri(u,diag=F)]
#
# #get cov matrix
#  len_m<-as.matrix(log(v[, eff_length]%o%v[, eff_length]))
#  gcc_m<-as.matrix(log(v[, gc]%o%v[, gc]))
#  map_m<-as.matrix(log(v[, cov]%o%v[, cov]))
#
# #centralize cov matrix of enz, gcc
#  len_m<-(len_m-mean(c(len_m)))/sd(c(len_m))
#  gcc_m<-(gcc_m-mean(c(gcc_m)))/sd(c(gcc_m))
#
# #change matrix into vector
#  len_vec<-len_m[upper.tri(len_m,diag=F)]
#  gcc_vec<-gcc_m[upper.tri(gcc_m,diag=F)]
#  map_vec<-map_m[upper.tri(map_m,diag=F)]
#
# #fit Poisson regression: u~len+gcc+offset(map)
#  fit<-glm(u_vec~len_vec+gcc_vec+offset(map_vec),family="poisson")
#
# #summary(fit)
#  coeff <- round(fit$coeff,4)
#  res <- round(u/exp(coeff[1]+coeff[2]*len_m+coeff[3]*gcc_m+map_m), 4)
#  data.table(id1=colnames(res), res)->res
#  melt(res, id.var="id1", value.name="nlinks_norm", variable.name="id2")->res
#  res[nlinks_norm > 0]
# }

def normalize_mat_cis()



def normalize_cis(binhic: dict, chrs=None, percentile=0, omit_smallest=0):
    if chrs is None:
        chrs = binhic["bin"]["chr"].unique().compute()

    # Fragments
    mf = binhic["bin"].query("chr in @chrs", local_dict={"chrs": chrs})
    # Binned HiCs
    ab = binhic["mat"]

    if percentile > 0:
        threshold = np.percentile(mf["eff_lenght"].compute().values, percentile + 1)
        mf = mf.query("eff_length >= @threshold", local_dict={"threshold": threshold})

    if omit_smallest > 0:
        eff_length = mf["eff_length"].compute().values
        eff_length.sort()
        threshold = eff_length[min(omit_smallest, eff_length.shape[0] - 1)]
        mf = mf.query("eff_length >= @threshold", local_dict={"threshold": threshold})

    #   y[, id1:=stri_c(sep=":", chr1, bin1)]
    #   y[, id2:=stri_c(sep=":", chr2, bin2)]
    left = ab.query("chr1 == chr2")
    merged = dd.merge(left, mf.rename(columns={"chr": "chr1", "bin": "bin1"}), on=["chr1", "bin1"]).merge(
        mf.rename(columns={"chr": "chr1", "bin": "bin1"}), on=["chr1", "bin1"])
    merged = merged.query("~((chr1 == chr2) & (bin1 == bin2))")
    matrix = dd.pivot_table(merged, index=["chr1", "chr2", "bin1"], columns=["bin2"], values=["nlinks"])




# # format intrachromosomal Hi-C matrix for HiCNorm
#  rbindlist(mclapply(mc.cores=ncores, chrs, function(i) {
#   ab[chr1 == chr2 & chr1 == i & id1 %in% mf$id & id2 %in% mf$id]->abf
#   dcast.data.table(abf, id1 ~ id2, value.var="nlinks", fill=as.integer(0))->mat
#
#   u<-as.matrix(mat[, setdiff(names(mat), "id1"), with=F])
#   setkey(mf, id)[colnames(u)]->v
#   normalize_mat_cis(u,v)
#  }))->nhic
#  mf[, data.table(key="id1", id1=id, chr1=chr, bin1=bin)][setkey(nhic, id1)]->nhic
#  mf[, data.table(key="id2", id2=id, chr2=chr, bin2=bin)][setkey(nhic, id2)]->nhic
#  nhic[, dist := abs(bin2 - bin1)][]
# }