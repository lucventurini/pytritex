# helper function to supply scaffolds order on the short of 5AS (a single bin in the Chapman et al. POPSEQ map of bread wheat)
correct_5A<-function(hic_info, bam, map, assembly){
 fread(cmd=paste("samtools view -q 30 -F260", bam, "| cut -f 1,3,4"))->i90sam
 setnames(i90sam, c("i90_marker", "scaffold", "pos"))

 i90map<-fread(map, select=c(1,3,4))
 setnames(i90map, c("i90_marker", "i90_alphachr", "i90_cM"))
 i90map[, n := .N, by=i90_marker]
 i90map[n == 1][, n := NULL][]->i90map
 setnames(chrNames(species="wheat"), c("i90_alphachr", "i90_chr"))[i90map, on="i90_alphachr"]->i90map
 i90map[i90sam, on="i90_marker"][i90_chr == 5, .(orig_scaffold=scaffold, orig_pos=pos, cM=i90_cM)]->z
 copy(z)->i90k_5A

 assembly$info[, .(scaffold, orig_scaffold, orig_start, orig_pos=orig_start)][z, on=c("orig_scaffold", "orig_pos"), roll=T]->z
 z[, .(scaffold, pos=orig_pos-orig_start+1, cM)]->i90k_5A

 z[, .(chr=5, icM=median(cM)), key="scaffold"][hic_info, on=c("scaffold", "chr")]->hic_info
 mx <- hic_info[chr == 5, max(na.omit(icM))]
 hic_info[chr == 5, cM := mx - icM]->hic_info
 hic_info[, icM := NULL][]
}

# Plot the biases computed by find_inversions() along the genome
plot_inversions<-function(inversions, hic_map, bad=c(), species, file, height=5, width=20){
 yy<-inversions$summary
 w<-inversions$ratio
 yy[scaffold %in% bad]->yy2

 pdf(file, height=height, width=width)
 lapply(sort(unique(w$chr)), function(i){
  w[chr == i, plot(bin/1e6, r, pch=20, las=1, bty='l', ylab="r", xlab="genomic position (Mb)", main=alphachr[i])]
  hic_map$agp[gap == T & agp_chr == chrNames(agp=T, species=species)[chr == i]$agp_chr, abline(lwd=1, col='gray', v=(agp_start+agp_end)/2e6)]
  if(nrow(yy2) > 0 & i %in% yy2$chr){
   yy2[chr == i, rect(col="#FF000011", agp_start/1e6, -1000, agp_end/1e6, 1000), by=scaffold]
  }
 })
 dev.off()
}