# Intra-chromosomal normalization, taken from HiCNorm https://academic.oup.com/bioinformatics/article/28/23/3131/192582
normalize_mat_cis<-function(u, v){
 u_vec<-u[upper.tri(u,diag=F)]

#get cov matrix
 len_m<-as.matrix(log(v[, eff_length]%o%v[, eff_length]))
 gcc_m<-as.matrix(log(v[, gc]%o%v[, gc]))
 map_m<-as.matrix(log(v[, cov]%o%v[, cov]))

#centralize cov matrix of enz, gcc
 len_m<-(len_m-mean(c(len_m)))/sd(c(len_m))
 gcc_m<-(gcc_m-mean(c(gcc_m)))/sd(c(gcc_m))

#change matrix into vector
 len_vec<-len_m[upper.tri(len_m,diag=F)]
 gcc_vec<-gcc_m[upper.tri(gcc_m,diag=F)]
 map_vec<-map_m[upper.tri(map_m,diag=F)]

#fit Poisson regression: u~len+gcc+offset(map)
 fit<-glm(u_vec~len_vec+gcc_vec+offset(map_vec),family="poisson")

#summary(fit)
 coeff <- round(fit$coeff,4)
 res <- round(u/exp(coeff[1]+coeff[2]*len_m+coeff[3]*gcc_m+map_m), 4)
 data.table(id1=colnames(res), res)->res
 melt(res, id.var="id1", value.name="nlinks_norm", variable.name="id2")->res
 res[nlinks_norm > 0]
}

# format intrachromosomal Hi-C matrix for HiCNorm 
normalize_cis<-function(binhic, ncores=1, chrs=NULL, percentile=0, omit_smallest=0){
 if(is.null(chrs)){
  chrs <- unique(binhic$bins$chr)
 }

 mf <- binhic$bins[chr %in% chrs]
 ab <- binhic$mat

 if(percentile > 0){
  mf[eff_length >= quantile(eff_length, 0:100/100)[percentile + 1]]->mf
 }

 if(omit_smallest > 0){
  setorder(mf, eff_length)
  mf[, omit_idx := 1:.N, by=chr]
  mf[omit_idx > omit_smallest][, omit_idx := NULL] -> mf
 }

 rbindlist(mclapply(mc.cores=ncores, chrs, function(i) {
  ab[chr1 == chr2 & chr1 == i & id1 %in% mf$id & id2 %in% mf$id]->abf
  dcast.data.table(abf, id1 ~ id2, value.var="nlinks", fill=as.integer(0))->mat

  u<-as.matrix(mat[, setdiff(names(mat), "id1"), with=F])
  setkey(mf, id)[colnames(u)]->v
  normalize_mat_cis(u,v)
 }))->nhic
 mf[, data.table(key="id1", id1=id, chr1=chr, bin1=bin)][setkey(nhic, id1)]->nhic
 mf[, data.table(key="id2", id2=id, chr2=chr, bin2=bin)][setkey(nhic, id2)]->nhic
 nhic[, dist := abs(bin2 - bin1)][]
}

# format interchromosomal Hi-C matrix for HiCNorm 
normalize_trans<-function(binhic, ncores=1, chrs=NULL, percentile=0){
 if(is.null(chrs)){
  chrs <- unique(binhic$bins$chr)
 }

 mf <- binhic$bins[chr %in% chrs]
 ab <- binhic$mat

 if(percentile > 0){
  mf[eff_length >= quantile(eff_length, 0:100/100)[percentile + 1]]->mf
 }

 rbindlist(mclapply(mc.cores=ncores, chrs, function(i) rbindlist(lapply(setdiff(1:21, i), function(j) {
  ab[chr1 == i & chr2 == j & id1 %in% mf$id & id2 %in% mf$id]->abf
  dcast.data.table(abf, id1 ~ id2, value.var="nlinks", fill=as.integer(0))->mat
  setkey(mf, id)[mat$id1]->v_a
  u<-as.matrix(mat[, setdiff(names(mat), "id1"), with=F])
  setkey(mf, id)[colnames(u)]->v_b
  normalize_mat_trans(u, v_a, v_b)
 }))))->nhic
 mf[, data.table(key="id1", id1=id, chr1=chr, bin1=bin)][setkey(nhic, id1)]->nhic
 mf[, data.table(key="id2", id2=id, chr2=chr, bin2=bin)][setkey(nhic, id2)]->nhic
}