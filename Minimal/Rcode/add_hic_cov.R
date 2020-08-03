# Calculate physical coverage with Hi-C links in sliding windows along the scaffolds
add_hic_cov<-function(assembly, scaffolds=NULL, binsize=1e3, binsize2=1e5, minNbin=50, innerDist=1e5, cores=1){

 info<-assembly$info

 if("mr" %in% colnames(info) | "mri" %in% colnames(info)){
  stop("assembly$info already has mr and/or mri columns; aborting.")
 }

 fpairs<-assembly$fpairs

 if(is.null(scaffolds)){
  scaffolds <- info$scaffold
  null=T
 } else {
  info[scaffold %in% scaffolds]->info
  fpairs[scaffold1 %in% scaffolds]->fpairs
  null=F
 }

 fpairs[scaffold1 == scaffold2 & pos1 < pos2][, .(scaffold = scaffold1, bin1 = pos1 %/% binsize * binsize, bin2 =pos2 %/% binsize * binsize)]->f
 f[bin2 - bin1 > 2*binsize]->f
 f[, i := 1:.N]
 f[, b := paste0(scaffold, ":", bin1 %/% binsize2)]
 setkey(f, b)

 rbindlist(mclapply(mc.cores=cores, unique(f$b), function(j){
  f[j][, .(scaffold=scaffold, bin=seq(bin1+binsize, bin2-binsize, binsize)), key=i][, .(n=.N), key=.(scaffold, bin)]
 }))->ff

 if(nrow(ff) > 0){
  ff[, .(n=sum(n)), key=.(scaffold, bin)]->ff
  info[, .(scaffold, length)][ff, on="scaffold"]->ff
  ff[, d := pmin(bin, (length-bin) %/% binsize * binsize)]
  ff[, nbin := .N, key="scaffold"]
  ff[, mn := mean(n), key=d]
  ff[, r := log2(n/mn)]
  ff[nbin >= minNbin, .(mr=suppressWarnings(min(r))), key=scaffold][order(mr)]->z
  ff[nbin > minNbin & d >= innerDist, .(mri=suppressWarnings(min(r))), key=scaffold]->zi
  z[ff, on="scaffold"]->ff
  zi[ff, on="scaffold"]->ff
  z[info, on="scaffold"]->info_mr
  zi[info_mr, on="scaffold"]->info_mr
 } else {
  copy(info) -> info_mr
  info_mr[, c("mri", "mr") := list(NA, NA)]
 }

 if(null){
  assembly$info=info_mr
  assembly$cov=ff

  assembly$binsize <- binsize
  assembly$minNbin <- minNbin
  assembly$innerDist <- innerDist

  assembly
 } else {
  list(info=info_mr, cov=ff)
 }
}