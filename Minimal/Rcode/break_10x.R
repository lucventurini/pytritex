# Iteratively break scaffolds using 10X physical coverage. Proceed until no more breakpoints are found.
break_10x<-function(assembly, prefix="scaffold_corrected", species="wheat", ratio=-3, interval=5e4, minNbin=20, dist=2e3, slop=1e3, 
		    intermediate=F, ncores=1, maxcycle=Inf){

 if(dist <= 2 * slop){
  dist <- 2 * slop + 1
  cat(paste0("Setting dist to ", 2 * slop + 1, "\n"))
 }
 find_10x_breaks(assembly=assembly, interval=interval, minNbin=minNbin, dist=dist, ratio=ratio) -> breaks
 i <- 1
 lbreaks <- list()
 lbreaks[[1]] <- breaks
 if(intermediate){
  assemblies <- list()
  assemblies[[1]] <- assembly
 }
 while(nrow(breaks) > 0 & i <= maxcycle){
  cat(paste0("Cycle ", i, ": ", breaks[, .N], " break points detected\n"))
  i <- i + 1
  break_scaffolds(breaks=breaks, assembly=assembly, species=species, prefix=paste0(prefix, "_"), slop=slop, cores=ncores) -> assembly
  if(intermediate){
   assemblies[[i]] <- assembly
  }
  find_10x_breaks(assembly=assembly, interval=interval, minNbin=minNbin, dist=dist, ratio=ratio) -> breaks
  if(is.null(breaks)){
   break
  }
  lbreaks[[i]] <- breaks
 }
 if(!intermediate){
  list(assembly=assembly, breaks=lbreaks)
 } else {
  list(assemblies=assemblies, breaks=lbreaks)
 }
}