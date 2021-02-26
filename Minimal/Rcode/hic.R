# Wrapper function for Hi-C mapping
make_hic_info<-function(cluster_info, super_global, chrs){
 s<-super_global$super_info
 s[!duplicated(s$chr),]->s
 s[chr %in% chrs]->s

 super_global$membership[, .(cluster, super, bin, rank, backbone)]->tmp
 tmp[super %in% s$super]->tmp
 tmp[, super := NULL]
 setnames(tmp, c("cluster", "hic_bin", "hic_rank", "hic_backbone"))
 tmp[cluster_info, on="cluster"]->cluster_info
 cluster_info[order(chr, hic_bin, hic_rank, cluster)]
}


# Wrapper function for Hi-C mapping
make_hic_map<-function(hic_info, links, ncores=1, maxiter=100, known_ends=T){

 copy(links)->hl
 copy(hic_info)->info

 setnames(hl, c("scaffold1", "scaffold2"), c("cluster1", "cluster2"))
 setnames(info, "scaffold", "cluster")
 setkey(info, "cluster")
 chrs <- info[!is.na(chr), unique(chr)]

 make_hic_info(info, 
  super_global<-make_super(hl, cluster_info=info, cores=ncores, maxiter=maxiter,
			   known_ends=known_ends, path_max=length(chrs)), chrs=chrs)->res
 res[order(chr, hic_bin)][, .(scaffold=cluster, chr, cM, hic_bin, hic_backbone, hic_rank)][!is.na(hic_bin)]
}


# Create a heatmap plot of Hi-C contact matrix
contact_matrix<-function(hic_map, links, file, species, chrs=NULL, boundaries=T, grid=NULL, ncol=100, trafo=NULL, v=NULL){
 colorRampPalette(c("white", "red"))(ncol)->whitered

 if(is.null(chrs)){
  chrs <- unique(links$chr1) 
 }

 binsize <- min(links[dist > 0]$dist)

 links[chr1 %in% chrs, .(chr=chr1, bin1, bin2, l=log10(nlinks_norm))]->z
 if(is.null(trafo)){
  z[, col := whitered[cut(l, ncol, labels=F)]]
 } else {
  z[, col := whitered[cut(trafo(l), ncol, labels=F)]]
 }

 chrNames(agp=T, species=species)[z, on="chr"]->z

 pdf(file)
 lapply(chrs, function(i){
  z[chr == i, plot(0, las=1, type='n', bty='l', xlim=range(bin1/1e6), ylim=range(bin2/1e6), xlab="position (Mb)", ylab="position (Mb)", col=0, main=chrNames(species=species)[chr == i, alphachr])]
  if(boundaries){
   hic_map$agp[gap == T & agp_chr == chrNames(agp=T, species=species)[chr == i, agp_chr], abline(lwd=1, col='gray', v=(agp_start+agp_end)/2e6)]
   hic_map$agp[gap == T & agp_chr == chrNames(agp=T, species=species)[chr == i, agp_chr], abline(lwd=1, col='gray', h=(agp_start+agp_end)/2e6)]
  }
  if(!is.null(grid)){
   max(hic_map$agp[gap == T & agp_chr == chrNames(agp=T, species=species)[chr == i, agp_chr]]$agp_end)->end
   abline(v=seq(0, end, grid)/1e6, col="blue", lty=2)
  }
  if(!is.null(v)){
   abline(v=v/1e6, col="blue", lty=2)
  }
  z[chr == i, rect((bin1-binsize)/1e6, (bin2-binsize)/1e6, bin1/1e6, bin2/1e6, col=col, border=NA)]
 })
 dev.off()
}


#calculate Hi-C physical coverage for pseudomolecules
hic_cov_psmol <- function(hic_map, binsize=1e3, binsize2=1e5, maxdist=1e6, 
			  cores=1){
 fpairs <- copy(hic_map$links)
 info <- hic_map$chrlen
 setnames(fpairs, c("start1", "start2"), c("pos1", "pos2"))

 fpairs[chr1 == chr2 & pos1 < pos2][, .(chr = chr1, bin1 = pos1 %/% binsize * binsize, bin2 =pos2 %/% binsize * binsize)]->f
 f[bin2 - bin1 > 2*binsize & bin2 - bin1 <= maxdist]->f
 f[, i := 1:.N]
 f[, b := paste0(chr, ":", bin1 %/% binsize2)]
 setkey(f, b)

 rbindlist(mclapply(mc.cores=cores, unique(f$b), function(j){
  f[j][, .(chr, bin=seq(bin1+binsize, bin2-binsize, binsize)), key=i][, .(n=.N), key=.(chr, bin)]
 }))->ff

 ff[, .(n=sum(n)), key=.(chr, bin)]->ff
 info[, .(chr, length)][ff, on="chr"]->ff
 ff[, d := pmin(bin, (length-bin) %/% binsize * binsize)]
 ff[, nbin := .N, key="chr"]
 ff[, mn := mean(n), key=d]
 ff[, r := log2(n/mn)][]
}
