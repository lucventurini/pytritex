# Export Hi-C map information for use in the R Shiny app
inspector_export<-function(hic_map, assembly, inversions, species, file, trafo=NULL){
 hic_map$hic_1Mb$norm -> x
 if(!is.null(trafo)){
  x[, nlinks_norm := trafo(nlinks_norm)]
 }
 wheatchr(agp=T, species=species)[, .(chr1=chr, agp_chr)][x, on="chr1"]->x
 x[, c("chr1", "chr2", "id2", "id1", "dist") := list(NULL, NULL, NULL, NULL, NULL)]
 
 inversions$ratio[, .(agp_chr, bin, r)] -> y

 copy(hic_map$agp)[, chr := NULL][]->a
 assembly$info[, .(scaffold, popseq_chr=popseq_alphachr)][a, on="scaffold"]->a
 saveRDS(list(agp=a, binsize=1e6, links=x, ratio=y), file=file)
}

# Plot heatmap of interchromosomal Hi-C contact matrix
interchromosomal_matrix_plot<-function(hic_map, binsize=1e6, file, species, resolution=72,
				       height=3000, width=3000, col=NULL, ncol=100, cex=2,
				       oma=c(5,7,5,1)){
 colorRampPalette(c("white", "red"))(ncol)->whitered
 copy(hic_map$hic_1Mb$mat)->z
 z[, nl := log10(nlinks)]
 if(is.null(col)){
  z[, col := whitered[cut(nl, ncol, labels=F)]]
 } else{
  cc <- col
  z[, col := cc] 
 }

 nrow(chrNames(species=species)) -> n

 z[, chr1 := as.character(chr1)]
 z[, chr2 := as.character(chr2)]

 copy(hic_map$chrlen)[, chr := paste(chr)][paste(1:n), on="chr"]$length -> ll

 png(file, res=resolution, height=height, width=width)
 par(mar=c(0.1, 0.1, 0.1, 0.1))
 layout(matrix(ncol=n, nrow=n, 1:(n*n), byrow=T), widths=ll, heights=ll)
 par(oma=oma)
 lapply(1:n, function(i){
  lapply(1:n, function(j){
   z[paste(i), on="chr1"][paste(j), on="chr2"] -> zz
   zz[, plot(0, las=1, xaxt='n', yaxt='n', type='n', xlim=range(bin1/1e6), ylim=range(bin2/1e6), xlab="", ylab="", col=0, main="")]
   zz[, rect((bin1-binsize)/1e6, (bin2-binsize)/1e6, bin1/1e6, bin2/1e6, col=col, border=NA)]
   if(i == 1){
    title(xpd=NA, wheatchr(species=species)[chr == j]$alphachr, line=1, cex.main=cex, font=2)
   }
   if(j == 1){
    mtext(wheatchr(species=species)[chr == i]$alphachr, side=2, line=1.5, cex=cex*0.6, font=2, las=1, xpd=NA)
   }
  })
 })
 dev.off()
}