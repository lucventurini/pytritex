# Find breakpoints in chimeric scaffolds based on drop in 10X molecule coverage
find_10x_breaks<-function(assembly, scaffolds=NULL, interval = 5e4, minNbin = 20, dist = 5e3, ratio = -3){
 cov <- copy(assembly$molecule_cov)
 if(!is.null(scaffolds)){
  cov[scaffolds, on="scaffold"]->cov
 } 
 cov[, b := bin %/% interval * interval]
 cov[nbin >= minNbin & pmin(bin, length - bin) >= dist & r <= ratio]->e
 if(nrow(e) == 0){
  return(NULL)
 }
 e[order(r)][, idx := 1:.N, by=.(scaffold, b)][idx == 1]->e
 setnames(e, "bin", "br")[]
 e[]
}