library(data.table)
library(parallel)
library(igraph)
library(stringi)
library(zoo)

source("Rcode/utils/wheatChrom.R")

plot_popseq_carma_tcc<-function(scaffold, breaks=NULL, page, info, cssaln, tcc_pos, span, molcov, species){
  chrNames(species=species)->wheatchr

  i <- scaffold

  if(is.null(tcc_pos) || nrow(tcc_pos) == 0){
   hic <- F
  } else {
   hic <- T
  }

  if(is.null(molcov) || nrow(molcov) == 0){
   tenex <- F
  } else {
   tenex <- T
  }

  nrow <- 1
  if(hic){
   nrow <- nrow + 1
  }
  if(tenex){
   nrow <- nrow + 1
  }

  par(mfrow=c(nrow,3))

  par(oma=c(1,0,3,0))
  if(hic){
   par(mar=c(1,4,4,1))
  } else {
   par(mar=c(4,4,4,1))
  }
  par(cex=1)

  if(!is.null(breaks)){
   br <- breaks[i, on="scaffold", br]/1e6
  } else if (hic & autobreaks){
   span[i, on="scaffold"][order(r)][1, bin]/1e6->br
  } else {
   br <- NULL
  }

  if(hic){
   info[scaffold == i, paste0(", bad Hi-C: ", bad_hic)] -> badhic
   xlab<-""
  } else {
   badhic <- ""
   xlab<-"position in scaffold (Mb)"
  }

  cssaln[i, on="scaffold"]->z
  l=info[i, on="scaffold"][, length]/1e6

  if(species %in% c("rye", "lolium", "barley")){
   ymin <- 8
  } else {
   ymin <- 24
  }

  if(mbscale == F){
   ylab = "POPSEQ chromosome"
  } else {
   ylab = "Morex V2 chromosome"
  }

  plot(xlim=c(0,l), 0, ylim=c(ymin,1), bty='l', col=0, yaxt='n', pch=20, xlab=xlab, ylab=ylab)
  if(!is.null(br)){
   abline(v=br, lwd=3, col="blue")
  }
  if(mbscale == F){
   title("POPSEQ", line=1)
  } else {
   title("Morex V2", line=1)
  }
  z[, points(pos/1e6, popseq_chr, pch=20)]
  axis(2, las=1, wheatchr$chr, wheatchr$alphachr)

  if(mbscale == F){
   info[scaffold == i, title(outer=T, line=0, paste("page: ", page, ", ", scaffold, sep="", " (", round(l,1), " Mb),\n",  
			   "bad POPSEQ: ", bad_popseq, ", bad CARMA: ", bad_sorted, badhic))]
  } else {
   info[scaffold == i, title(outer=T, line=0, paste("page: ", page, "\n", scaffold, sep="", " (", round(l,1), " Mb)"))]
  }

  w<-info[i, on="scaffold"]$popseq_chr
  if(is.na(w)){
   plot(xlab='', ylab='', type='n', 0, axes=F)
  } else {
   if(mbscale == F){
    ylab = "genetic position (cM)"
   } else {
    ylab = "pseudomolecule position (Mb)"
   }
   z[popseq_chr == w][, plot(pch=20, xlab=xlab, ylab=ylab, xlim=c(0,l), bty='l', las=1, pos/1e6, popseq_cM)]
   if(mbscale == F){
    title(main="Genetic positions of CSS contigs\nfrom major chromosome")
   } else {
    title(main="Physical positions of 1 kb single-copy tags\nfrom major chromosome")
   }
   if(hic){
    abline(v=br, lwd=3, col="blue")
   }
  } 

  if(species != 'lolium'){
   if(mbscale == F){
    ylab="flow-sorting chromosome"
    ttt="flow-sorting CARMA"
   } else {
    ylab="chromsome arm"
    ttt="Hi-C based chromosome arm assignment"
   }
   plot(xlim=c(0,l), 0, bty='l', ylim=c(ymin,1), yaxt='n', pch=20, col=0, xlab=xlab, ylab=ylab)
   title(ttt, line=1)
   if(!is.null(br)){
    abline(v=br, lwd=3, col="blue")
   }
   if(species == "wheat"){
    z[, points(pos/1e6, sorted_chr, col=ifelse(sorted_alphachr == "3B" | sorted_arm == "S", 1, 2), pch=20)]
    legend(horiz=T, "bottomleft", pch=19, col=1:2, bg="white", legend=c("short arm / 3B", "long arm"))
   }
   if(species == "rye"){
    z[, points(pos/1e6, sorted_chr, col=1, pch=20)]
   }
   if(species == "barley"){
    z[, points(pos/1e6, sorted_chr, col=ifelse(sorted_alphachr == "1H" | sorted_arm == "S", 1, 2), pch=20)]
    if(mbscale == F){
     llx =  "short arm / 1H"
    } else {
     llx = "short arm"
    }
    legend(horiz=T, "bottomleft", pch=19, col=1:2, bg="white", legend=c(llx, "long arm"))
   }
   axis(2, las=1, wheatchr$chr, wheatchr$alphachr)
  } else {
   plot(0, type='n', xlab="", ylab="", axes=F)
  }

  if(hic){
   par(mar=c(4,4,4,1))
   plot(0, xlim=c(0,l), col=0, bty='l', ylab='Hi-C chromosome', yaxt='n', ylim=c(ymin,1), xlab="position in scaffold (Mb)")
   title("Interchromosomal Hi-C links", line=1)
   tcc_pos[i, on="scaffold1", points(pos1/1e6, col="#00000003", chr2, pch=20)]
   abline(v=br, lwd=3, col="blue")
   axis(2, las=1, wheatchr$chr, wheatchr$alphachr)

   span[i, on="scaffold"]->w
   if(nrow(w[!is.na(n)]) == 0){
    plot(xlab='', ylab='', type='n', 0, axes=F)
   } else {
    w[order(bin), plot(bty='l', bin/1e6, n, log='y', xlim=c(0,l), col=0, ylab='coverage', las=1, xlab="position in scaffold (Mb)")]
    title("Intrascaffold physical Hi-C coverage", line=1)
    abline(v=br, lwd=3, col="blue")
    w[order(bin), points(bin/1e6, n, xlim=c(0,l), type='l', lwd=3, col=1)]
   }

   if(nrow(w[!is.na(n)]) == 0){
    plot(xlab='', ylab='', type='n', 0, axes=F)
   } else {
    w[order(bin), plot(bty='l', bin/1e6, r, xlim=c(0,l), col=0, ylab='log2(observed/expected ratio)', las=1, xlab="position in scaffold (Mb)")]
    title("Hi-C expected vs. observed coverage", line=1)
    abline(v=br, lwd=3, col="blue")
    w[order(bin), points(bin/1e6, r, xlim=c(0,l), type='l', lwd=3, col=1)]
   }
  }

  if(tenex){
   molcov[i, on="scaffold"]->w
   if(nrow(w) == 0){
    plot(xlab='', ylab='', type='n', 0, axes=F)
   } else {
    w[order(bin), plot(bty='l', bin/1e6, n, log='y', xlim=c(0,l), col=0, ylab='coverage', las=1, xlab="position in scaffold (Mb)")]
    title("10X molecule coverage", line=1)
    if(!is.null(br)){
     abline(v=br, lwd=3, col="blue")
    }
    w[order(bin), points(bin/1e6, n, xlim=c(0,l), type='l', lwd=3, col=1)]
   }

   if(nrow(w) == 0){
    plot(xlab='', ylab='', type='n', 0, axes=F)
   } else {
    w[order(bin), plot(bty='l', bin/1e6, r, xlim=c(0,l), col=0, ylab='log2(observed/expected ratio)', las=1, xlab="position in scaffold (Mb)")]
    title("10X expected vs. observed coverage", line=1)
    if(!is.null(br)){
     abline(v=br, lwd=3, col="blue")
    }
    w[order(bin), points(bin/1e6, r, xlim=c(0,l), type='l', lwd=3, col=1)]
   }
  }
 }