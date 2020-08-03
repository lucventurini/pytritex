library(data.table)
library(parallel)
library(igraph)
library(stringi)
library(zoo)

source("Rcode/make_super_path.R")

# Set up Hi-C graph structure, determine ends of chromosome from genetic map, run Hi-C ordering for each chromosome
make_super<-function(hl, cluster_info, prefix="super", cores=1, paths=T, path_max=0, known_ends=F, 
		     maxiter=100, verbose=T){

 hl[cluster1 %in% cluster_info[excluded == F]$cluster & cluster2 %in% cluster_info[excluded == F]$cluster]->hl
 hl[cluster1 < cluster2]->e
 graph.edgelist(as.matrix(e[, .(cluster1, cluster2)]), directed=F)->g
 E(g)$weight<-e$weight

 data.table(cluster=V(g)$name, super=paste(prefix, sep="_", clusters(g)$membership))->mem
 cluster_info[mem, on="cluster"]->mem

 mem[, .(super_size=.N, length=.N, chr=unique(na.omit(chr))[1], cM=mean(na.omit(cM))), keyby=super]->info
 mem[, .(cluster1=cluster, super)][hl, on="cluster1"]->e

 list(super_info=info, membership=mem, graph=g, edges=e)->s

 if(paths){
  if(path_max > 0){
   idx<-head(s$super_info[order(-length)], n=path_max)$super
  } else {
   idx<-s$super_info$super
  }
  rbindlist(mclapply(mc.cores=cores, idx, function(i) {
   start <- end <- NULL
   # Take terminal nodes from genetic map
   if(known_ends){
    s$mem[super == i & !is.na(cM)][order(cM)]$cluster->x
    start=x[1]
    end=tail(x,1)
   }
   make_super_path(s, idx=i, start=start, end=end, maxiter=maxiter, verbose=verbose)->x
   cat(paste0("Chromosome ", head(s$mem[super == i]$chr, n=1), " finished.\n"))
   if(verbose){
    cat(paste0("Chromosome ", head(s$mem[super == i]$chr, n=1), " finished.\n"))
   }
   x
  }))[s$membership, on="cluster"]->s$membership
 }

 s
}