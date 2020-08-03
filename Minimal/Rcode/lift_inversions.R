# Find inverted scaffolds based on directionality biases
find_inversions<-function(hic_map, links, species, chrs=NULL, cores=1, winsize=15, maxdist=1e8, threshold=40, factor=100){
 
 if(is.null(chrs)){
  unique(links$chr1) -> chrs
 }

 links[chr1 %in% chrs]->zz
 zz[, l := log(factor*nlinks_norm/sum(nlinks_norm)*1e6)]->zz
 hic_map$chrlen[, .(chr1=chr, length)][zz, on="chr1"]->zz
 zz[dist <= bin1 & dist <= bin2 & (length - bin1) >= dist & (length - bin2) >= dist]->zz
 zz[l >= 0 & dist <= maxdist, .(r = sum(sign(bin1 - bin2) * l)), key=.(chr=chr1, bin=bin1)]->w
 chrNames(agp=T, species=species)[w, on="chr"]->w

 data.table(w, rbindlist(mclapply(mc.cores=cores, chrs, function(i){
  ww<-w[chr == i][order(bin)]
  data.table(
   cc=ww[, rollapply(r, width=winsize, FUN=function(x) cor(x, 1:length(x)), align="left", fill=NA)],
   lm=ww[, rollapply(r, width=winsize, FUN=function(x) lm(data=.(a=x, b=1:winsize/winsize), a~b)$coefficient[2], align="left",fill=NA)]
  )
 })))->w

 hic_map$agp[gap == F, .(agp_chr, bin=agp_start, scaffold)]->x
 copy(w)->y
 y[, bin := bin + 1]
 x[y, on=c("agp_chr", "bin"), roll=T]->y
 y[agp_chr != "chrUn"]->y
 y[, .(s=sd(r)), key=scaffold][!is.na(s)][order(-s)]->yy
 hic_map$hic_map[, .(scaffold, chr=consensus_chr, hic_bin)][yy, on="scaffold"]->yy
 hic_map$agp[, .(scaffold, agp_start, agp_end)][yy, on="scaffold"]->yy
 list(ratio=w, summary=yy)
}
 
# Flip the orientation of groups of scaffolds
correct_multi_inversions<-function(hic_map, ranges, species){
 copy(ranges)->y
 chrNames(agp=T, species=species)[, .(consensus_chr=chr, agp_chr)][hic_map$hic_map, on="consensus_chr"]->m
 y[, i:=1:.N]

 y[, .(b=seq(start,end)), by=.(agp_chr, i)][, .N, key=.(agp_chr, b)][N > 1]->dups
 if(nrow(dups) > 0){
  stop(paste("Overlapping ranges specified: bin(s)",
        paste(dups[, paste(sep=":", agp_chr, b)], collapse=", "), "are contained in more than one inversion."))
 }

 y[, .(agp_chr=agp_chr[1], hic_bin = start:end, new_bin = end:start, new=T), by=i][, i := NULL][m, on=c("agp_chr", "hic_bin")][is.na(new), new := F]->m
 m[new == T, hic_bin := new_bin]
 m[new == T, consensus_orientation := consensus_orientation * -1]
 m[new == T & is.na(consensus_orientation), consensus_orientation := -1]
 m[, c("new", "new_bin", "agp_chr") := list(NULL, NULL, NULL)]
 make_agp(m, gap_size=hic_map$gap_size, species=species)->a

 list()->new
 new$hic_map <- m
 new$agp <- a$agp
 new$gap_size <- copy(hic_map$gap_size)
 new$chrlen <- copy(hic_map$chrlen)
 new$binsize <- copy(hic_map$binsize)
 new$hic_map_bin <- copy(hic_map$hic_map_bin)
 new$max_cM_dist <- copy(hic_map$max_cM_dist)
 new$min_nfrag_bin <- copy(hic_map$min_nfrag_bin)
 new$agp_bed <- a$agp_bed
 new$corrected_multi_inversions <- copy(ranges)
 new
}