# Convert original scaffolds positions to AGP positions
orig_scaffold_to_agp<-function(assembly, map_10x=NULL, hic_map){
 if(!is.null(map_10x)){
  super <- map_10x
  assembly$info[, .(orig_scaffold, orig_scaffold_start=orig_start, 
			      orig_scaffold_end=orig_end, scaffold)]->zz
  super$agp[, .(scaffold, scaffold_length=length, super, super_start, super_end,
		       super_orientation=orientation)][zz, on="scaffold"]->zz
  hic_map$agp[scaffold != "gap", .(agp_chr, agp_start0 = agp_start, agp_end0 = agp_end, super=scaffold, agp_orientation = ifelse(is.na(orientation), 1, orientation))]->a
  a[zz, on="super"]->zz
  zz[, orientation := agp_orientation * super_orientation]
  zz[agp_orientation == 1, agp_start := agp_start0 - 1 + super_start] 
  zz[agp_orientation == 1, agp_end := agp_start0 - 1 + super_end] 
  zz[agp_orientation == -1, agp_start := agp_end0 + 1 - super_end] 
  zz[agp_orientation == -1, agp_end := agp_end0 + 1 - super_start] 
  zz[, agp_start0 := NULL]
  zz[, agp_end0 := NULL]
  zz[]
 } else {
  assembly$info[, .(orig_scaffold, orig_scaffold_start=orig_start, 
			       orig_scaffold_end=orig_end, scaffold)]->zz
  hic_map$agp[scaffold != "gap", .(agp_chr, agp_start, agp_end, scaffold, agp_orientation = ifelse(is.na(orientation), 1, orientation))]->a
  a[zz, on="scaffold"]->zz
  zz[, orientation := agp_orientation]
  zz[]
 }
}

# Flip AGP orientation of specified scaffolds
correct_inversions<-function(hic_map, scaffolds, species){
 copy(hic_map$hic_map)->z
 z[scaffold %in% scaffolds, consensus_orientation := ifelse(is.na(consensus_orientation), -1, -1 * consensus_orientation)]
 make_agp(z, gap_size=hic_map$gap_size, species=species)->a

 list()->new
 new$hic_map <- z
 new$agp <- a$agp
 new$gap_size <- copy(hic_map$gap_size)
 new$chrlen <- copy(hic_map$chrlen)
 new$binsize <- copy(hic_map$binsize)
 new$hic_map_bin <- copy(hic_map$hic_map_bin)
 new$max_cM_dist <- copy(hic_map$max_cM_dist)
 new$min_nfrag_bin <- copy(hic_map$min_nfrag_bin)
 new$agp_bed <- a$agp_bed
 new$corrected_inversions <- scaffolds
 new
}


# Convert Hi-C map table into an AGP
make_agp<-function(hic_map_oriented, gap_size=100, species){

 hic_map_oriented[, .(scaffold, chr = consensus_chr,
	  popseq_cM=ifelse(consensus_chr == popseq_chr | is.na(consensus_chr), popseq_cM, NA),
	  scaffold_length = length, hic_bin, orientation=consensus_orientation)]->z
 chrNames(agp=T, species=species)[z, on="chr"]->z

 z[, agp_chr := "chrUn"]
 z[!is.na(hic_bin), agp_chr := sub("NA", "Un", paste0("chr", alphachr))]
 z[, alphachr := NULL]
 z[order(agp_chr, hic_bin, chr, popseq_cM, -scaffold_length)]->z
 z[, index := 2*1:.N-1]
 z[, gap := F]
 rbind(z, data.table(scaffold="gap", gap=T, chr=NA, popseq_cM=NA, scaffold_length = gap_size, hic_bin = NA, orientation = NA, agp_chr=z$agp_chr, index=z$index+1))->z
 z[order(index)][, head(.SD, .N-1), by=agp_chr]->z
 z[, agp_start := cumsum(c(0, scaffold_length[1:(.N-1)]))+1, by = agp_chr]
 z[, agp_end := cumsum(scaffold_length), by = agp_chr]

 z[, .(scaffold=scaffold, bed_start=0, bed_end=scaffold_length, name=scaffold, score=1, strand=ifelse(is.na(orientation) | orientation == 1, "+", "-"), agp_chr=agp_chr)]->agp_bed

 list(agp=z, agp_bed=agp_bed)
}


make_agp_membership<-function(membership, gap_size=100){

  membership[, .(scaffold, length, super, bin, orientation)]->z

  setorder(z, super, bin, -length)
  z[, index := 2*1:.N-1]
  z[, gap := F]
  rbind(z, data.table(scaffold="gap", gap=T, super=z$super, bin=NA, length = gap_size,orientation = NA, index=z$index+1))->z
  z[order(index)][, head(.SD, .N-1), by=super]->z
  z[, n := .N, key=super]
  z[n > 1, super_start := cumsum(c(0, length[1:(.N-1)])) + 1, by = super]
  z[n == 1, super_start := 1]
  z[, super_end := cumsum(length), by = super]
  z[, n := NULL]

  z[, .(scaffold=scaffold, bed_start=0, bed_end=length, name=scaffold, score=1, strand=ifelse(is.na(orientation) | orientation == 1, "+", "-"), super=super)]->agp_bed

  list(agp=z, agp_bed=agp_bed)
 }