# Read alignment of 10X super-scaffolds to another scaffolds. First used for comparison to NRGene assemblies, hence the name.
read_nrgene_aln <- function(assembly, super, paf, fai0=NULL, minq=20){
 paf["P", on="type"]->z
 z[mapq >= minq]->z
 if(!is.null(fai0)){
  fai0[, .(query=scaffold0, scaffold)][z, on="query"]->z
 } else {
  z[, scaffold := query]
 }
 z[, .(orig_scaffold=scaffold, orig_scaffold_start=query_start, orig_scaffold_end=query_end,
       nrgene=reference, nrgene_length=reference_length, orientation,
       nrgene_start=reference_start,  nrgene_end=reference_end, mapq, alnlen)]->zz
 assembly$info[, .(scaffold, scaffold_length=length, orig_scaffold, orig_scaffold_start=orig_start, orig_start)]->b
 b[zz, on=c("orig_scaffold", "orig_scaffold_start"), roll=T]->zz
 zz[, scaffold_start :=  orig_scaffold_start - orig_start + 1]
 zz[, scaffold_end := orig_scaffold_end - orig_start + 1]
 super$membership[, .(scaffold, super, super_pos, super_orientation=orientation)][zz, on="scaffold"]->zz
 zz[super_orientation == 1, super_start := super_pos - 1 + scaffold_start]
 zz[super_orientation == -1, super_end := super_pos + (scaffold_length - scaffold_start)]
 zz[super_orientation == 1, super_end := super_pos - 1 + scaffold_end]
 zz[super_orientation == -1, super_start := super_pos + (scaffold_length - scaffold_end)]
 zz[super_orientation == -1, orientation := -1L * orientation]
 super$info[, .(super, super_length=length, super_chr=chr)][zz, on="super"]->zz
 zz[]
}

# Plot alignment of 10X super-scaffolds to another assemblies. First used for comparison to NRGene assemblies, hence the name.
plot_nrgene_aln<-function(data, nrgene_aln, mem=NULL, file, ncores=1, minlen=2e4){
 copy(data)[, plot_idx := 1:.N] -> data

 plotfu<-function(data, aln, mem=NULL, minlen){
  cat(paste0(data$plot_idx, " ", data$super, "\n"))
  aln[data[, .(super, nrgene)], on=c("nrgene", "super")][alnlen >= minlen] -> xx
  if(nrow(xx) > 0){
   xlim <- c(0,0)
   ylim <- c(0,0)
   quantile(xx$super_start, 0:50/50)[2]/1e6 -> xlim[1]
   quantile(xx$super_end, 0:50/50)[50]/1e6 -> xlim[2]
   quantile(xx$nrgene_start, 0:50/50)[2]/1e6 -> ylim[1]
   quantile(xx$nrgene_end, 0:50/50)[50]/1e6 -> ylim[2]
   xx[, plot(0, xlim=xlim, type='n', ylim=ylim, las=1, xlab="10X super-scaffold position (Mb)",
	     ylab="NRGene scaffold position (Mb)", bty='l')]
   data[, title(main=paste0(sub("super", "super_scaffold", super), ", ", round(super_length/1e6, 1), " Mb\n",
		            sub("scaffold", "NRGene scaffold", nrgene), ", ", round(nrgene_length/1e6, 1), " Mb"))]
   if(!is.null(mem)){
    abline(v=mem[data$super, on="super"]$super_pos/1e6, col="gray")
   }
   xx[, idx := 1:.N]
   xx[orientation == 1, lines(c(super_start/1e6, super_end/1e6), c(nrgene_start/1e6, nrgene_end/1e6)), by=idx]
   xx[orientation == -1, lines(c(super_start/1e6, super_end/1e6), c(nrgene_end/1e6, nrgene_start/1e6)), by=idx]
  } else {
   plot(0, type='n', axes=F, xlab="", ylab="")
   data[, title(main=paste0(sub("super", "super_scaffold", super), ", ", round(super_length/1e6, 1), " Mb\n",
		            sub("scaffold", "NRGene scaffold", nrgene), ", ", round(nrgene_length/1e6, 1), " Mb"))]
  }
 }
 parallel_plot(data=data, group="idx", cores=ncores, aln=nrgene_aln, minlen=minlen, mem=mem,
		file=file, height=700, width=700, res=150, plot_function=plotfu)
}