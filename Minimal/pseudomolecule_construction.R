library(data.table)
library(parallel)
library(igraph)
library(stringi)
library(zoo)

cat("Disabling multi-threading and optimizations of recent data.table versions\n")
cat("setDTthreads(1)\n")
cat("options(data.table.optimize=1)\n")
setDTthreads(1)
options(data.table.optimize=1)
options(error=function() { traceback(2); if(!interactive()) quit("no", status = 1, runLast = FALSE) })

files.sources = list.files(path="Rcode/utils", pattern="*.R$", full.names=TRUE, recursive=FALSE)
sapply(files.sources, source)
files.sources = c()
files.sources = list.files(path="Rcode", pattern="*.R$", full.names=TRUE, recursive=FALSE)
sapply(files.sources, source)

# Like scaffold_10x, but using a Hi-C map of the scaffolds, instead of a genetic map to prune erroneous edges
scaffold_10x_hic <- function(assembly, prefix="super", min_npairs=5, max_dist=1e5, min_nmol=6, min_nsample=2, popseq_dist=5, max_dist_orientation=5, ncores=1, verbose=T, raw=F, unanchored=T, hic_map=NULL){

 if(verbose){
  cat("Finding links.\n")
 }

 copy(assembly$info) -> info
 if(!is.null(hic_map)){
  info[, c("popseq_chr", "popseq_cM") := list(NULL, NULL)]
  hic_map$agp[, .(popseq_chr = chr, popseq_cM = hic_bin, scaffold)][info, on="scaffold"]->info
 } 

 assembly$molecules[npairs >= min_npairs] -> z
 info[, .(scaffold, scaffold_length=length)][z, on="scaffold"]->z
 z[end <= max_dist | scaffold_length - start <= max_dist]->z
 z[, nsc := length(unique(scaffold)), key=.(barcode, sample)]
 z[nsc >= 2]->z
 z[, .(scaffold1=scaffold, npairs1=npairs, pos1=as.integer((start+end)/2), sample, barcode)]->x 
 z[, .(scaffold2=scaffold, npairs2=npairs, pos2=as.integer((start+end)/2), sample, barcode)]->y
 y[x, on=.(sample, barcode), allow.cartesian=T][scaffold1 != scaffold2]->xy
 xy -> link_pos
 xy[, .(nmol=.N), key=.(scaffold1, scaffold2, sample)]->w
 w[nmol >= min_nmol]->ww
 ww[, .(nsample = length(unique(sample))), key=.(scaffold1, scaffold2)][nsample >= min_nsample]->ww2
 info[, .(scaffold1=scaffold, popseq_chr1=popseq_chr, length1=length, popseq_pchr1=popseq_pchr, popseq_cM1=popseq_cM)][ww2, on="scaffold1"]->ww2
 info[, .(scaffold2=scaffold, popseq_chr2=popseq_chr, length2=length, popseq_pchr2=popseq_pchr, popseq_cM2=popseq_cM)][ww2, on="scaffold2"]->ww2
 ww2[, same_chr := popseq_chr2 == popseq_chr1]
 ww2[, weight := -1 * log10((length1 + length2) / 1e9)]

 if(popseq_dist > 0){
  ww2[(popseq_chr2 == popseq_chr1 & abs(popseq_cM1 - popseq_cM2) <= popseq_dist)] -> links
 } else {
  ww2 -> links
 }

 ex <- c()

 run <- T

 if(verbose){
  cat("Finding initial super-scaffolds (ln 2311).\n")
 }
 # remove nodes until graph is free of branches of length > 1
 if(!raw){
  cat("Non-raw procedure.\n")
  while(run){
   make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
   copy(out$membership) -> m
   copy(out$info) -> res
   m[unique(m[rank > 1][, .(super, bin)]), on=c("super", "bin")]->a
   a[rank == 0] -> add
   if(nrow(add) == 0){
    run <- F
   } else{
    c(ex, add$scaffold) -> ex
   }
  }
 } else {
  cat("Raw procedure.\n") 
  make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
  out$m -> m
  out$info -> res
 }

 if(!raw){
  if(verbose){
   cat("Tip removal.\n")
  }
  # remove short tips of rank 1 
  links[!scaffold2 %in% ex][, .(degree=.N), key=.(scaffold=scaffold1)]->b
  a <- b[m[m[rank == 1][, .(super=c(super, super, super), bin=c(bin, bin-1, bin+1))], on=c("super", "bin")], on="scaffold"]
  a[degree == 1 & length <= 1e4]$scaffold->add
  if(length(add) > 0){
   ex <- c(ex, add)
   make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
   out$m -> m
  }

  # remove short tips/bulges of rank 1 
  m[rank == 1 & length <= 1e4]$scaffold -> add
  ex <- c(ex, add)
  if(length(add) > 0){
   ex <- c(ex, add)
   make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
   out$m -> m
  }

  # resolve length-one-bifurcations at the ends of paths
  m[rank > 0][bin == 2 | super_nbin - 1 == bin ][, .(super, super_nbin, type = bin == 2, scaffold, length, bin0=bin)]->x
  unique(rbind(
  m[x[type == T, .(super, bin0, bin=1)], on=c("super", "bin")],
  m[x[type == F, .(super, bin0, bin=super_nbin)], on=c("super", "bin")]
  ))->a
  a[, .(super, bin0, scaffold2=scaffold, length2=length)][x, on=c("super", "bin0")][, ex := ifelse(length >= length2, scaffold2, scaffold)]$ex -> add

  if(length(add) > 0){
   ex <- c(ex, add)
   make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
   out$membership -> m
  }

  # remove short tips/bulges of rank 1 
  m[rank == 1 & length <= 1e4]$scaffold -> add
  if(length(add) > 0){
   ex <- c(ex, add)
   make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
   out$m -> m
  }
  
  # remove tips of rank 1 
  links[!scaffold2 %in% ex][, .(degree=.N), key=.(scaffold=scaffold1)]->b
  b[m[rank == 1], on="scaffold"][degree == 1]$scaffold -> add
  if(length(add) > 0){
   ex <- c(ex, add)
   make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
   out$m -> m
  }

  # remove remaining nodes of rank > 0
  m[rank > 0]$scaffold -> add
  if(length(add) > 0){
   ex <- c(ex, add)
   make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
   out$m -> m
  }

  if(popseq_dist > 0 & unanchored == T){
   if(verbose){
    cat("Including unanchored scaffolds.\n")
   }
   # use unanchored scaffolds to link super-scaffolds
   ww2[is.na(popseq_chr1), .(scaffold_link=scaffold1, link_length=length1, scaffold1=scaffold2)]->x
   ww2[is.na(popseq_chr1), .(scaffold_link=scaffold1, scaffold2=scaffold2)]->y
   x[y, on="scaffold_link", allow.cartesian=T][scaffold1 != scaffold2]->xy

   m[, .(scaffold1=scaffold, super1=super, chr1=chr, cM1=cM, size1=super_nbin, d1 = pmin(bin - 1, super_nbin - bin))][xy, on="scaffold1"]->xy
   m[, .(scaffold2=scaffold, super2=super, chr2=chr, cM2=cM, size2=super_nbin, d2 = pmin(bin - 1, super_nbin - bin))][xy, on="scaffold2"]->xy
   xy[super2 != super1 & d1 == 0 & d2 == 0 & size1 > 1 & size2 > 1 & chr1 == chr2]->xy
   xy[scaffold1 < scaffold2, .(nscl=.N), scaffold_link][xy, on="scaffold_link"]->xy
   xy[nscl == 1] -> xy
   xy[super1 < super2][, c("n", "g"):=list(.N, .GRP), by=.(super1, super2)][order(-link_length)][!duplicated(g)]->zz

   sel <- zz[, .(scaffold1=c(scaffold_link, scaffold_link, scaffold1, scaffold2),
	  scaffold2=c(scaffold1, scaffold2, scaffold_link, scaffold_link))]
   rbind(links, ww2[sel, on=c("scaffold1", "scaffold2")])->links2

   make_super_scaffolds(links=links2, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
   out$m -> m

   #resolve branches
   m[rank > 0][bin == 2 | super_nbin - 1 == bin ][, .(super, super_nbin, type = bin == 2, scaffold, length, bin0=bin)]->x
   unique(rbind(
   m[x[type == T, .(super, bin0, bin=1)], on=c("super", "bin")],
   m[x[type == F, .(super, bin0, bin=super_nbin)], on=c("super", "bin")]
   ))->a
   a[, .(super, bin0, scaffold2=scaffold, length2=length)][x, on=c("super", "bin0")][, ex := ifelse(length >= length2, scaffold2, scaffold)]$ex -> add

   if(length(add) > 0){
    ex <- c(ex, add)
    make_super_scaffolds(links=links2, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
    out$membership -> m
   }

   # remove remaining nodes of rank > 0
   m[rank > 0]$scaffold -> add
   if(length(add) > 0){
    ex <- c(ex, add)
    make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
   }
  }

  out$m -> m
  out$info -> res
 }

 if(verbose){
  cat("Orienting scaffolds.\n")
 }
 #orient scaffolds
 m[super_nbin > 1, .(scaffold1=scaffold, bin1=bin, super1=super)][link_pos, on="scaffold1", nomatch=0]->a
 m[super_nbin > 1, .(scaffold2=scaffold, bin2=bin, super2=super)][a, on="scaffold2", nomatch=0]->a
 a[super1 == super2 & bin1 != bin2]->a
 a[, d := abs(bin2 - bin1)]->a
 a[d <= max_dist_orientation]->a
 a[, .(nxt=mean(pos1[bin2 > bin1]), prv=mean(pos1[bin2 < bin1])), key=.(scaffold=scaffold1)]->aa
 info[, .(scaffold, length)][aa, on="scaffold"]->aa
 aa[!is.nan(prv) & !is.nan(nxt), orientation := ifelse(prv <= nxt, 1, -1)]
 aa[is.nan(prv), orientation := ifelse(length - nxt <= nxt, 1, -1)]
 aa[is.nan(nxt), orientation := ifelse(length - prv <= prv, -1, 1)]

 aa[, .(orientation, scaffold)][m, on="scaffold"]->m
 m[, oriented := T]
 m[is.na(orientation), oriented := F]
 m[is.na(orientation), orientation := 1]
 setorder(m, super, bin, rank)

 m[, super_pos := 1 + cumsum(c(0, length[-.N])), by=super]

 if(verbose){
  cat("Anchoring super-scaffolds.\n")
 }
 # assign super scaffolds to genetic positions
 m[!is.na(chr), .(nchr=.N), key=.(chr, super)][, pchr := nchr/sum(nchr), by=super][order(-nchr)][!duplicated(super)]->y
 m[!is.na(cM)][y, on=c("chr", "super")]->yy
 y[res, on="super"]->res
 yy[, .(cM=mean(cM), min_cM=min(cM), max_cM=max(cM)), key=super][res, on='super']->res
 setorder(res, -length)

 make_agp_member(m, gap=100)->a

 list(membership=m, info=res, 
      agp=a$agp, agp_bed=a$agp_bed)
}
