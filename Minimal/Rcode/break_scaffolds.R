library(data.table)
library(parallel)
library(igraph)
library(stringi)
library(zoo)


# Break scaffolds at specified points and lift positional information to updated assembly
break_scaffolds<-function(breaks, assembly, prefix, slop, cores=1, species, 
			  regex1="(^.*[^-0-9])(([0-9]+)(-[0-9]+)?$)", 
			  regex2="(^.*[^-0-9])([0-9]+(-[0-9]+)?$)"){
 info <- assembly$info
 cov <- assembly$cov
 fpairs <- assembly$fpairs
 cssaln <- assembly$cssaln
 molecules <- assembly$molecules
 
 br<-copy(breaks)
 info[, .(scaffold, orig_scaffold, old_scaffold=scaffold, orig_start, orig_end, length)]->fai

 cat("Split scaffolds\n")
 j <- 0
 while(nrow(br) > 0){
  j <- j + 1
  o <- nrow(fai)
  fai[br, on="scaffold"] -> br
  br[, orig_br := orig_start + br - 1]
  br[order(scaffold, br)] -> br
  br[duplicated(scaffold)] -> nbr
  br[!duplicated(scaffold)] -> br

  max(as.integer(sub(regex1, "\\3", fai$scaffold)))->maxidx
  br[, idx := 3*1:.N-2]
  br[, scaffold1 := paste0(prefix, maxidx+idx)]
  br[, start1 := 1]
  br[, end1 := pmax(0, br - slop - 1)]
  br[, scaffold2 := paste0(prefix, maxidx+idx+1)]
  br[, start2 := pmax(1, br - slop)]
  br[, end2 := pmin(br + slop - 1, length)]
  br[, scaffold3 := paste0(prefix, maxidx+idx+2)]
  br[, start3 := pmin(length + 1, br + slop)]
  br[, end3 := length]
  br[, length1 := 1 + end1 - start1]
  br[, length2 := 1 + end2 - start2]
  br[, length3 := 1 + end3 - start3]
  rbind(
   br[, .(scaffold=scaffold1, length=length1, orig_scaffold, orig_start=orig_start+start1-1, orig_end=orig_start+end1-1, old_scaffold)],
   br[, .(scaffold=scaffold2, length=length2, orig_scaffold, orig_start=orig_start+start2-1, orig_end=orig_start+end2-1, old_scaffold)],
   br[, .(scaffold=scaffold3, length=length3, orig_scaffold, orig_start=orig_start+start3-1, orig_end=orig_start+end3-1, old_scaffold)],
   fai[!scaffold %in% br$scaffold, .(orig_scaffold, length,  orig_start, orig_end, old_scaffold,
			      scaffold=paste0(prefix, sub(regex2, "\\2", scaffold)))]

  ) -> fai

  cat(paste0("Iteration ", j, " finished. "))
  fai[length > 0]->fai
  cat(paste0("The number of scaffolds increased from ", o, " to ", nrow(fai), ".\n"))
  fai[, .(scaffold, orig_scaffold, orig_start, orig_br=orig_start)][nbr[, .(orig_scaffold, orig_br)], on=c("orig_scaffold", "orig_br"), roll=T]->nbr
  nbr[, br := orig_br - orig_start + 1]
  nbr[, .(scaffold, br)] -> br
 }

 fai[, split := F]
 fai[old_scaffold %in% breaks$scaffold, split := T]
 fai[, old_scaffold := NULL]

 assembly_new<-list(info=fai)

 cat("Transpose cssaln\n")
 copy(cssaln) -> z
 z[, scaffold_length := NULL]
 z[, scaffold := NULL]
 fai[, .(scaffold, scaffold_length=length, orig_scaffold, orig_start, orig_pos=orig_start)][z, on=c("orig_scaffold", "orig_pos"), roll=T]->z
 z[, pos := orig_pos - orig_start + 1]
 z[, orig_start := NULL]
 assembly_new$cssaln <- z
 
 if("fpairs" %in% names(assembly) && nrow(fpairs) > 0){
  cat("Transpose fpairs\n")
  assembly$fpairs[, .(orig_scaffold1, orig_scaffold2, orig_pos1, orig_pos2)]->z
  fai[, .(scaffold1=scaffold, orig_scaffold1=orig_scaffold, orig_start1=orig_start, orig_pos1=orig_start)][z, on=c("orig_scaffold1", "orig_pos1"), roll=T]->z
  fai[, .(scaffold2=scaffold, orig_scaffold2=orig_scaffold, orig_start2=orig_start, orig_pos2=orig_start)][z, on=c("orig_scaffold2", "orig_pos2"), roll=T]->z
  z[, pos1 := orig_pos1 - orig_start1 + 1]
  z[, pos2 := orig_pos2 - orig_start2 + 1]
  z[, orig_start1 := NULL]
  z[, orig_start2 := NULL]
  assembly_new$fpairs <- z
 } else {
  assembly_new$fpairs <- data.table()
 }

 if("molecules" %in% names(assembly) && nrow(molecules) > 0){
  cat("Transpose molecules\n")
  copy(molecules) -> z
  z[, scaffold := NULL]
  fai[, .(scaffold, orig_scaffold, orig_start, s_length=orig_end - orig_start + 1, orig_pos=orig_start)][z, on=c("orig_scaffold", "orig_start"), roll=T]->z
  z[, start := orig_start - orig_pos + 1]
  z[, end := orig_end - orig_pos + 1]
  z[end <= s_length]->z
  z[, orig_pos := NULL]
  z[, s_length  := NULL]
  assembly_new$molecules <- z
 } else {
  assembly_new$molecules <- data.table()
 }

 assembly_new$breaks <- breaks

 cat("Anchor scaffolds\n")
 anchor_scaffolds(assembly_new, popseq=assembly$popseq, species=species)->assembly_new
 
 if("mr_10x" %in% names(assembly_new$info)){
  assembly_new$info[, mr_10x := NULL]
 }

 if("mr" %in% names(assembly_new$info)){
  assembly_new$info[, mr := NULL]
  assembly_new$info[, mri := NULL]
 }

 if("cov" %in% names(assembly) & nrow(fpairs) > 0){
  cat("Hi-C coverage\n")
  add_hic_cov(assembly_new, scaffolds=fai[split == T]$scaffold, binsize=assembly$binsize, minNbin=assembly$minNbin, innerDist=assembly$innerDist, cores=cores)->cov
  assembly$cov[!scaffold %in% breaks$scaffold]->x
  x[, scaffold:=paste0(prefix, sub(regex2, "\\2", scaffold))]
  if(nrow(cov$cov) > 0){
   rbind(x, cov$cov)->assembly_new$cov
  } else {
   x -> assembly_new$cov
  }
  info[!scaffold %in% breaks$scaffold]->x
  x[, scaffold:=paste0(prefix, sub(regex2, "\\2", scaffold))]
  x[, split := F]
  rbind(x[, names(cov$info), with=F], cov$info)->assembly_new$info
 } else {
  assembly_new$cov <- data.table()
 }

 if("molecule_cov" %in% names(assembly) & nrow(molecules) > 0){
  cat("10X molecule coverage\n")
  add_molecule_cov(assembly_new, scaffolds=fai[split == T]$scaffold, binsize=assembly$mol_binsize, cores=cores)->cov
  info[!breaks$scaffold, on="scaffold"]->x
  x[, scaffold := paste0(prefix, sub(regex2, "\\2", scaffold))]
  x[, split := F]
  rbind(x[, names(cov$info), with=F], cov$info)->assembly_new$info
  assembly_new$mol_binsize <- assembly$mol_binsize

  assembly$molecule_cov[!breaks$scaffold, on="scaffold"]->x
  x[, scaffold := paste0(prefix, sub(regex2, "\\2", scaffold))]
  if(nrow(cov$molecule_cov) > 0){
   rbind(x, cov$molecule_cov)->assembly_new$molecule_cov
  } else {
   x -> assembly_new$molecule_cov 
  }
 } else {
  assembly_new$molecule_cov <- data.table()
 }

 assembly_new$binsize <- assembly$binsize
 assembly_new$innerDist <- assembly$innerDist
 assembly_new$minNbin <- assembly$minNbin

 assembly_new
}