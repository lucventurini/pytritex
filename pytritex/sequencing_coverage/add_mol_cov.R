library(data.table)

# Calculate the physical coverage with 10X molecules in sliding windows along the scaffolds
add_molecule_cov<-function(assembly, scaffolds=NULL, binsize=200, cores=1){

 info <- assembly$info

 if(is.null(assembly$molecules) || nrow(assembly$molecules) == 0){
  stop("The assembly object does not have a molecule table; aborting.")
 }

 if("mr_10x" %in% colnames(info)){
  stop("assembly$info already has mr_10x column; aborting.")
 }

 if(is.null(scaffolds)){
  scaffolds <- info$scaffold
  copy(assembly$molecules) -> mol
  null=T
 } else {
  info[scaffolds, on="scaffold"] -> info
  assembly$molecules[scaffolds, on="scaffold"] -> mol
  null=F
 }

 mol -> f
 f[, bin1 := start %/% binsize * binsize]
 f[, bin2 := end %/% binsize * binsize]
 f[bin2 - bin1 > 2 * binsize] -> f
 setkey(f, scaffold)
 f[, i := 1:.N]

 rbindlist(mclapply(mc.cores=cores, unique(f$scaffold), function(j){
  f[j][, .(scaffold, bin=seq(bin1+binsize, bin2-binsize, binsize)), key=i][, .(n=.N), key=.(scaffold, bin)]
 })) -> ff

 if(nrow(ff) > 0){
  info[, .(scaffold, length)][ff, on="scaffold"]->ff
  ff[, d := pmin(bin, (length-bin) %/% binsize * binsize)]
  ff[, nbin := .N, key="scaffold"]
  ff[, mn := mean(n), key=d]
  ff[, r := log2(n/mn)]

  ff[, .(mr_10x = min(r)), key=scaffold][info, on="scaffold"] -> info_mr
 } else {
  copy(info) -> info_mr
  info_mr[, mr_10x := NA]
 }

 if(null){
  assembly$info <- info_mr
  assembly$molecule_cov <- ff
  assembly$mol_binsize <- binsize

  assembly
 } else {
  list(info=info_mr, molecule_cov=ff)
 }
}