#function to plot contact matrix, Hi-C physical coverage and directionality bias, one page for each chromosome
big_hic_plot<-function(hic_map, cov, inv, species, file, chrs=NULL, cores=1, breaks=NULL){
 hmap <- hic_map

 ncol=100
 colorRampPalette(c("white", "red"))(ncol)->whitered
 links <- hmap$hic_1Mb$norm
 if(is.null(chrs)){
  chrs <- unique(links$chr1)
 }
 links[chr1 %in% chrs, .(chr=chr1, bin1, bin2, l=log10(nlinks_norm))]->z
 chrNames(agp=T, species=species)[z, on="chr"]->z
 binsize <- min(links[dist > 0]$dist)
 z[, col := whitered[cut(l, ncol, labels=F)]]
 nchr <- length(chrs)

 tempdir() -> d
 tmp <- paste0(d, "/plot")

 mclapply(mc.cores=cores, chrs, function(i){
  f <- paste0(tmp, "_", chrNames(agp=T, species=species)[chr == i, agp_chr], ".png")
  png(f, res=500, height=5/3*2500, width=5/3*1500)
  layout(matrix(1:(3*1), ncol=1), heights=rep(c(3,1,1), 1))
   
   z[chr == i, plot(0, las=1, type='n', bty='l', xlim=range(bin1/1e6), ylim=range(bin2/1e6), xlab="position (Mb)", ylab="position (Mb)", col=0, main=chrNames(species=species)[chr == i, alphachr])]
    hmap$agp[gap == T & agp_chr == chrNames(agp=T, species=species)[chr == i, agp_chr], abline(lwd=1, col='gray', v=(agp_start+agp_end)/2e6)]
    hmap$agp[gap == T & agp_chr == chrNames(agp=T, species=species)[chr == i, agp_chr], abline(lwd=1, col='gray', h=(agp_start+agp_end)/2e6)]
   z[chr == i, rect((bin1-binsize)/1e6, (bin2-binsize)/1e6, bin1/1e6, bin2/1e6, col=col, border=NA)]

   cov[chr == i][, plot(bin/1e6, type='n', las=1, xlab="genomic position (Mb)",
		       ylab="r", bty='l', r)]
   hmap$agp[chr == i, abline(col="gray", v=c(0,agp_end/1e6))]
   if(!is.null(breaks)){
    abline(col="blue", lty=2, v=breaks[chr == i, bin/1e6])
   }
   abline(h=-2, lty=2, col='red')
   cov[chr == i][, points(bin/1e6, r, pch=".")]
   title(main="physical Hi-C coverage")

   inv$ratio[chr == i][, plot(bin/1e6, type='n', las=1, 
			      xlab="genomic position (Mb)",
			      ylab="r", bty='l', r)]
   hmap$agp[chr == i, abline(col="gray", v=c(0,agp_end/1e6))]
   inv$ratio[chr == i][, points(bin/1e6, pch=".", r)]
   title(main="directionality bias")
  dev.off()
  system(paste("convert", f, sub(".png$", ".pdf", f)))
 })

 system(paste0("gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=", file, " ", 
	       paste(paste0(tmp, "_", 
			    setdiff(chrNames(agp=T, species=species)[chr %in% chrs]$agp_chr, "chrUn"), ".pdf"), collapse=' ')))
 unlink(tmp, force=T)
}
