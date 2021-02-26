
# Call Hi-C function for ordering, orient by ordering parts of scaffolds, create AGP table
hic_map<-function(info, assembly, frags, species, ncores=1, min_nfrag_scaffold=50, max_cM_dist = 20, 
		  binsize=5e5, min_nfrag_bin=30, gap_size=100, maxiter=100, orient=T, agp_only=F,
		  map=NULL, known_ends=T, orient_old=F, min_binsize=1e5, min_nbin=5){
 if(!agp_only){
  copy(info)->hic_info
  hic_info[, excluded := nfrag < min_nfrag_scaffold]
  
  assembly$fpairs[scaffold1 != scaffold2, .(nlinks=.N), key=.(scaffold1, scaffold2)]->hl
  hic_info[, .(scaffold1=scaffold, chr1=chr, cM1=cM)][hl, nomatch=0, on="scaffold1"]->hl
  hic_info[, .(scaffold2=scaffold, chr2=chr, cM2=cM)][hl, nomatch=0, on="scaffold2"]->hl
  hl[chr1 == chr2]->hl
  hl<-hl[abs(cM1-cM2) <= max_cM_dist | is.na(cM1) | is.na(cM2)]
  hl[, weight:=-log10(nlinks)]

  cat("Scaffold map construction started.\n")
  make_hic_map(hic_info=hic_info, links=hl, ncores=ncores, known_ends=known_ends)->hic_map
  cat("Scaffold map construction finished.\n")

  if(orient){
   if(orient_old){
    options(scipen = 1000)
    frags[, .(nfrag=.N), keyby=.(scaffold, pos = start %/% binsize * binsize)]->fragbin
    fragbin[, id := paste(sep=":", scaffold, pos)]
    fragbin<-hic_info[excluded == F, .(scaffold, chr, cM)][fragbin, on="scaffold", nomatch=0]

    assembly$fpairs[, .(nlinks=.N), keyby=.(scaffold1, pos1 = pos1 %/% binsize * binsize, scaffold2, pos2 = pos2 %/% binsize * binsize)]->binl
    binl[, id1 := paste(sep=":", scaffold1, pos1)]
    binl[, id2 := paste(sep=":", scaffold2, pos2)]
    binl[id1 != id2]->binl
    fragbin[, .(id1=id, chr1=chr, cM1=cM)][binl, on="id1"]->binl
    fragbin[, .(id2=id, chr2=chr, cM2=cM)][binl, on="id2"]->binl
    binl[, c("scaffold1", "scaffold2", "pos1", "pos2") := list(NULL, NULL, NULL, NULL)]
    setnames(binl, c("id1", "id2"), c("scaffold1", "scaffold2"))

    cat("Scaffold bin map construction started.\n")
    fragbin[, .(scaffold=id, nfrag, chr, cM)]->hic_info_bin
    hic_info_bin[, excluded:=nfrag < min_nfrag_bin]
    binl[chr1 == chr2 & (abs(cM1-cM2) <= max_cM_dist | is.na(cM1) | is.na(cM2))]->binl
    binl[, weight:=-log10(nlinks)]

    make_hic_map(hic_info=hic_info_bin, links=binl, ncores=ncores, maxiter=maxiter, known_ends=known_ends)->hic_map_bin
    cat("Scaffold bin map construction finished.\n")

    w<-hic_map_bin[!is.na(hic_bin), .(id=scaffold, scaffold=sub(":.*$", "", scaffold), pos=as.integer(sub("^.*:", "", scaffold)), chr, hic_bin)]
    w<-w[, .(gbin=mean(na.omit(hic_bin)),
	     hic_cor=as.numeric(suppressWarnings(cor(method='s', hic_bin, pos, use='p')))), keyby=scaffold][!is.na(hic_cor)]
    hic_map[!is.na(hic_bin) & scaffold %in% w$scaffold][order(chr, hic_bin)]->z0
    z0[,.(scaffold1=scaffold[1:(.N-2)], scaffold2=scaffold[2:(.N-1)], scaffold3=scaffold[3:(.N)]), by=chr]->z
    z0[, data.table(key="scaffold1", scaffold1=scaffold, hic_bin1=hic_bin)][setkey(z, "scaffold1")]->z
    z0[, data.table(key="scaffold2", scaffold2=scaffold, hic_bin2=hic_bin)][setkey(z, "scaffold2")]->z
    z0[, data.table(key="scaffold3", scaffold3=scaffold, hic_bin3=hic_bin)][setkey(z, "scaffold3")]->z
    w[, data.table(key="scaffold1", scaffold1=scaffold, gbin1=gbin)][setkey(z, "scaffold1")]->z
    w[, data.table(key="scaffold2", scaffold2=scaffold, gbin2=gbin)][setkey(z, "scaffold2")]->z
    w[, data.table(key="scaffold3", scaffold3=scaffold, gbin3=gbin)][setkey(z, "scaffold3")]->z
    z[, cc:= apply(z[, .(hic_bin1, hic_bin2, hic_bin3, gbin1, gbin2, gbin3)],1,function(x) {
		    suppressWarnings(cor(x[1:3], x[4:6]))
		    })]
    z[, data.table(key="scaffold", scaffold=scaffold2, cc=ifelse(cc > 0, 1, -1))]->ccor
    ccor[w]->m
    m[, hic_orientation:=ifelse(hic_cor > 0, 1 * cc, -1 * cc)]
    m[, .(scaffold, hic_cor, hic_invert=cc, hic_orientation)][hic_map, on="scaffold"]->hic_map_oriented

    setnames(hic_map_oriented, "chr", "consensus_chr")
    setnames(hic_map_oriented, "cM", "consensus_cM")
    hic_map_oriented[, consensus_orientation := hic_orientation]
   } else {
    options(scipen = 1000)
    assembly$info[, .(scaffold, binsize=pmax(min_binsize, length %/% min_nbin))][frags, on='scaffold']->f
    f[, .(nfrag=.N), keyby=.(scaffold, binsize, pos = start %/% binsize * binsize)]->fragbin
    fragbin[, id := paste(sep=":", scaffold, pos)]
    fragbin<- hic_map[, .(scaffold, chr, cM=hic_bin)][fragbin, on="scaffold", nomatch=0]

    unique(fragbin[, .(scaffold1=scaffold, binsize1=binsize)])[assembly$fpairs, on='scaffold1']->fp
    unique(fragbin[, .(scaffold2=scaffold, binsize2=binsize)])[fp, on='scaffold2']->fp
    fp[, .(nlinks=.N), keyby=.(scaffold1, pos1 = pos1 %/% binsize1 * binsize1, scaffold2, pos2 = pos2 %/% binsize2 * binsize2)]->binl
    binl[, id1 := paste(sep=":", scaffold1, pos1)]
    binl[, id2 := paste(sep=":", scaffold2, pos2)]
    binl[id1 != id2]->binl

    fragbin[, .(id1=id, chr1=chr, cM1=cM)][binl, on="id1"]->binl
    fragbin[, .(id2=id, chr2=chr, cM2=cM)][binl, on="id2"]->binl
    binl[, c("scaffold1", "scaffold2", "pos1", "pos2") := list(NULL, NULL, NULL, NULL)]
    setnames(binl, c("id1", "id2"), c("scaffold1", "scaffold2"))

    fragbin[, .(scaffold=id, nfrag, chr, cM)]->hic_info_bin
    hic_info_bin[, excluded:=nfrag < min_nfrag_bin]
    binl[chr1 == chr2 & (abs(cM1-cM2) <= 2 | is.na(cM1) | is.na(cM2))]->binl
    binl[, weight:=-log10(nlinks)]

    make_hic_map(hic_info=hic_info_bin, links=binl, ncores=ncores, maxiter=maxiter, known_ends=known_ends)->hic_map_bin

    w<-hic_map_bin[!is.na(hic_bin), .(id=scaffold, scaffold=sub(":.*$", "", scaffold), pos=as.integer(sub("^.*:", "", scaffold)), chr, hic_bin)]
    w<-w[, .(gbin=mean(na.omit(hic_bin)),
	     hic_cor=as.numeric(suppressWarnings(cor(method='s', hic_bin, pos, use='p')))), keyby=scaffold][!is.na(hic_cor)]
    hic_map[!is.na(hic_bin) & scaffold %in% w$scaffold][order(chr, hic_bin)]->z0
    z0[,.(scaffold1=scaffold[1:(.N-2)], scaffold2=scaffold[2:(.N-1)], scaffold3=scaffold[3:(.N)]), by=chr]->z
    z0[, data.table(key="scaffold1", scaffold1=scaffold, hic_bin1=hic_bin)][setkey(z, "scaffold1")]->z
    z0[, data.table(key="scaffold2", scaffold2=scaffold, hic_bin2=hic_bin)][setkey(z, "scaffold2")]->z
    z0[, data.table(key="scaffold3", scaffold3=scaffold, hic_bin3=hic_bin)][setkey(z, "scaffold3")]->z
    w[, data.table(key="scaffold1", scaffold1=scaffold, gbin1=gbin)][setkey(z, "scaffold1")]->z
    w[, data.table(key="scaffold2", scaffold2=scaffold, gbin2=gbin)][setkey(z, "scaffold2")]->z
    w[, data.table(key="scaffold3", scaffold3=scaffold, gbin3=gbin)][setkey(z, "scaffold3")]->z
    z[, cc:= apply(z[, .(hic_bin1, hic_bin2, hic_bin3, gbin1, gbin2, gbin3)],1,function(x) {
		    suppressWarnings(cor(x[1:3], x[4:6]))
		    })]
    z[, data.table(key="scaffold", scaffold=scaffold2, cc=ifelse(cc > 0, 1, -1))]->ccor
    ccor[w]->m
    m[, hic_orientation:=ifelse(hic_cor > 0, 1 * cc, -1 * cc)]
    m[, .(scaffold, hic_cor, hic_invert=cc, hic_orientation)][hic_map, on="scaffold"]->hic_map_oriented
    hic_map_oriented[is.na(hic_orientation), hic_orientation := ifelse(hic_cor > 0, 1, -1)]
    setnames(hic_map_oriented, "chr", "consensus_chr")
    setnames(hic_map_oriented, "cM", "consensus_cM")
    hic_map_oriented[, consensus_orientation := hic_orientation]
   }
  } else {
   hic_map_oriented<-copy(hic_map)
   setnames(hic_map_oriented, "chr", "consensus_chr")
   setnames(hic_map_oriented, "cM", "consensus_cM")
   hic_map_oriented[, consensus_orientation := as.numeric(NA)]
   hic_map_oriented[, hic_cor := as.numeric(NA)]
   hic_map_oriented[, hic_invert := as.numeric(NA)]
   hic_map_oriented[, hic_orientation := as.numeric(NA)]
   hic_map_bin <- NA
  }

  if("orientation" %in% names(hic_info)){
   hic_info[, .(scaffold, old_orientation=orientation)][hic_map_oriented, on="scaffold"]->hic_map_oriented
   hic_map_oriented[!is.na(old_orientation), consensus_orientation := old_orientation]
   hic_map_oriented[, old_orientation := NULL]
  }

  hic_map_oriented[assembly$info, on="scaffold"]->hic_map_oriented
 } else {
  hic_map_oriented <- map$hic_map
  hic_map_bin <- map$hic_map_bin
  min_nfrag_scaffold <- map$min_nfrag_scaffold
  binsize <- map$binsize
  max_cM_dist <- map$max_cM_dist
  min_nfrag_bin <- map$min_nfrag_bin
  gap_size <- map$gap_size
 }

 make_agp(hic_map_oriented, gap_size=gap_size, species=species)->a

 a$agp[, .(length=sum(scaffold_length)), key=agp_chr]->chrlen
 chrlen[, alphachr := sub("chr", "", agp_chr)]
 chrNames(species=species)[chrlen, on="alphachr"]->chrlen
 chrlen[, truechr := !grepl("Un", alphachr)]
 chrlen[order(!truechr, chr)]->chrlen
 chrlen[, offset := cumsum(c(0, length[1:(.N-1)]))]
 chrlen[, plot_offset := cumsum(c(0, length[1:(.N-1)]+1e8))]

 list(agp=a$agp, agp_bed=a$agp_bed, chrlen=chrlen, hic_map=hic_map_oriented, hic_map_bin=hic_map_bin)->res
 invisible(lapply(sort(c("min_nfrag_scaffold", "max_cM_dist", "binsize", "min_nfrag_bin", "gap_size")), function(i){
  res[[i]] <<- get(i)
 }))
 res
}


# Add Hi-C link information to a hic_map object
add_psmol_fpairs<-function(assembly, hic_map, nucfile, map_10x=NULL, assembly_10x=NULL, cov=NULL){
 if(is.null(map_10x)){
  assembly$fpairs[, .(scaffold1, scaffold2, pos1, pos2)] -> z
 } else {
  assembly_10x$fpairs[, .(scaffold1, scaffold2, pos1, pos2)] -> z
 }
 hic_map$agp[agp_chr != "chrUn", .(chr, scaffold, orientation=orientation, agp_start, agp_end)]->a
 a[is.na(orientation), orientation := 1]
 setnames(copy(a), paste0(names(a), 1))[z, on="scaffold1"]->z
 setnames(copy(a), paste0(names(a), 2))[z, on="scaffold2"]->z
 z[orientation1 == 1, pos1 := agp_start1 - 1 + pos1]
 z[orientation1 == -1, pos1 := agp_end1 + 1 - pos1]
 z[orientation2 == 1, pos2 := agp_start2 - 1 + pos2]
 z[orientation2 == -1, pos2 := agp_end2 + 1 - pos2]
 z[!is.na(chr1) & !is.na(chr2), .(chr1, chr2, start1=pos1, start2=pos2)]->links

 n<-c("orig_scaffold", "orig_start", "orig_end", "frag_id", "nA", "nC", "nG", "nT", "nN", "length")
 nuc<-fread(nucfile, select=c(1:4,7:11,13), head=T, col.names=n)
 assembly$info[, .(scaffold, orig_scaffold, orig_start, off=orig_start)][nuc, on=c("orig_scaffold", "orig_start"), roll=T]->z
 z[, start := orig_start - off + 1]
 z[, end := orig_end - off + 1]
 if(!is.null(map_10x)){
  map_10x$agp[gap == F, .(super, orientation, super_start, super_end, scaffold)][z, on="scaffold"]->z
  z[orientation == 1, start := super_start - 1 + start]
  z[orientation == 1, end := super_start - 1 + end]
  z[orientation == -1, start := super_end - end + 1]
  z[orientation == -1, end := super_end - start + 1]
  z[, c("orientation", "super_start", "super_end", "scaffold") := list(NULL, NULL, NULL, NULL)]
  setnames(z, "super", "orig_scaffold")
  assembly_10x$info[, .(scaffold, start=orig_start, orig_start, orig_scaffold)][z, on=c("orig_scaffold", "start"), roll=T]->z
  z[, start := start - orig_start + 1]
  z[, end := end - orig_start + 1]
  z[, orig_start := NULL]
  z[, orig_scaffold := NULL]
 }
 a[z, on="scaffold", nomatch=0]->z
 z[orientation == 1, start := agp_start - 1 + start]
 z[orientation == -1, start := agp_end + 1 - start]
 z[orientation == 1, end := agp_start - 1 + end]
 z[orientation == -1, end := agp_end + 1 - end]
 z[, c("chr", "start", "end", n[4:length(n)]), with=F]->z
 z[, cov := 1]
 z->frags
 
 copy(hic_map) -> res
 res$links <- links
 res$frags <- frags
 res
}