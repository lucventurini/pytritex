

read_morex_aln <- function(assembly, super, paf, agp, fai0=NULL, minq=20){
 paf["P", on="type"]->z
 z[mapq >= minq]->z
 z[, seq := sub("-[0-9]+$", "", query)]
 if(!is.null(fai0)){
  fai0[, .(reference=scaffold0, scaffold)][z, on="reference"]->z
 } else {
  z[, scaffold := reference]
 }
 agp[z, on="seq"]->z
 z[, .(orig_scaffold=scaffold, orig_scaffold_start=reference_start, orig_scaffold_end=reference_end,
       orientation, agp_chr=chr, 
       agp_start = agp_start - 1 + query_start,
       agp_end = agp_start - 1 + query_end, seq, bac, hic_bin, cluster, mapq, alnlen)]->zz
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
 super$info[, .(super, super_chr=chr)][zz, on="super"]->zz
 zz[]
}


plot_morex_aln<-function(data, morex_aln, file, ncores=1, minlen=5e3){
 copy(data)[, idx := 1:.N] -> data

 plotfu<-function(data, aln, mem=NULL, minlen){
  cat(paste0(data$idx, " ", data$super, "\n"))
  aln[data$super, on="super"][agp_chr == super_chr][alnlen >= minlen] -> xx
  if(nrow(xx) > 0){
   quantile(xx$agp_start, 0:50/50)[c(2, 50)]/1e6 -> ylim
   xlim <- c(0, data$length)/1e6
   xx[, plot(0, xlim=xlim, type='n', ylim=ylim, las=1, xlab="10X super-scaffold position (Mb)",
	     ylab="Morex AGP position (Mb)", bty='l')]
   data[, title(main=paste0(sub("super", "super_scaffold", super), ", ", chr, "H, ", round(length/1e6, 1), " Mb"))]
   if(!is.null(mem)){
    abline(v=mem[data$super, on="super"]$super_pos/1e6, col="gray")
   }
   xx[, idx := 1:.N]
   xx[orientation == 1, lines(c(super_start/1e6, super_end/1e6), c(agp_start/1e6, agp_end/1e6)), by=idx]
   xx[orientation == -1, lines(c(super_start/1e6, super_end/1e6), c(agp_end/1e6, agp_start/1e6)), by=idx]
  } else {
   plot(0, type='n', axes=F, xlab="", ylab="")
   data[, title(main=paste0(sub("super", "super_scaffold", super), ", ", chr, "H, ", round(length/1e6, 1), " Mb",
			   "\nNo alignments to Morex pseudomolecules."))]
  }
 }
 parallel_plot(data=data, group="super", cores=ncores, aln=morex_aln, minlen=minlen,
		file=file, height=700, width=700, res=150, plot_function=plotfu)
}


# Lift alignment coordinate: query to AGP; reference (=IBSC2017 assembly) BACs to pseudonmolecule
lift_morex_aln_agp<-function(paf, morex_agp, scaffold_to_agp){
 paf[, .(seq=sub("-[0-9]+$", "", query), seq_start=query_start, seq_end=query_end,
	       orig_scaffold=reference, orig_scaffold_start = reference_start, orig_scaffold_end = reference_end,
	       aln_orientation = orientation, mapq, alnlen)]->p
 morex_agp[, .(seq, morex_chr=agp_chr, morex_agp_start0=agp_start)][p, on="seq"]->p
 p[, morex_agp_start := morex_agp_start0 - 1 + seq_start]
 p[, morex_agp_end := morex_agp_start0 - 1 + seq_end]
 scaffold_to_agp[, .(orig_scaffold, orig_pos = orig_scaffold_start, orig_scaffold_start, agp_chr, agp_start0=agp_start, agp_end0=agp_end, orientation)][p, on=c("orig_scaffold", "orig_scaffold_start"), roll=T]->p
 p[, aln_orientation := orientation * aln_orientation]
 p[orientation == 1, agp_start := agp_start0 + (orig_scaffold_start - orig_pos)]
 p[orientation == 1, agp_end := agp_start0 + (orig_scaffold_end - orig_pos)]
 p[orientation == -1, agp_start := agp_end0 - (orig_scaffold_end - orig_pos)]
 p[orientation == -1, agp_end := agp_end0 - (orig_scaffold_start - orig_pos)]
 p[, agp_start0 := NULL]
 p[, agp_end0 := NULL]
 p[, orig_pos := NULL]
 p[]
}