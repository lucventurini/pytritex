library(data.table)
library(parallel)
library(igraph)
library(stringi)
library(zoo)


source("Rcode/plot_ccarma.R")

# Create diagnostics plot for putative chimeras
plot_chimeras<-function(scaffolds, assembly, breaks=NULL, file, mindist=0, cores=20, species="wheat", mbscale=F, autobreaks=T){

 info<-assembly$info
 ff<-assembly$cov
 tcc_pos <- assembly$fpairs
 cssaln <- assembly$cssaln
 molcov <- assembly$molecule_cov

 width <- 2500
 height <- 1000
 res <- 150

 if(!is.null(tcc_pos) && nrow(tcc_pos) > 0){
  height <- height + 1000
 }
 if(!is.null(molcov) && nrow(molcov) > 0){
  height <- height + 1000
 }

 scaffolds <- scaffolds[!duplicated(scaffold), .(scaffold)]

 bad <- copy(scaffolds)
 out <- file
 bad[, f:=tempfile(fileext=".png"), by=scaffold]
 bad[, p:=tempfile(fileext=".pdf"), by=scaffold]
 mclapply(mc.cores=cores, 1:nrow(bad), function(i){
  file <- bad[i, f]
  pdf <- bad[i, p]
  s<-bad[i,scaffold]
  cat(paste0(i, " ", bad[i, scaffold], "\n"))
  png(file, height=height, res=res, width=width)
  plot_popseq_carma_tcc(s, breaks=breaks, page=i, info, cssaln, tcc_pos, ff[d >= mindist], molcov, species)
  dev.off()
  system(paste("convert", file, pdf))
  unlink(file)
 })
 system(paste0("gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=", out, " ", paste(bad$p, collapse=' ')))
 unlink(bad$p)
}