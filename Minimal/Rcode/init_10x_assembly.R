# Initialize an assembly object for super-scaffolds created by scaffold_10x(). Lift positional information from scaffolds to super-scaffolds.
init_10x_assembly<-function(assembly, map_10x, molecules=F){
 super <- map_10x

 copy(assembly$cssaln)->z
 z[, orig_scaffold := NULL]
 z[, orig_scaffold_length := NULL]
 z[, orig_pos := NULL]
 z[, scaffold_length := NULL]
 super$agp[, .(scaffold, super, super_start, super_end, orientation)][z, on="scaffold"]->z
 z[orientation == 1, pos := super_start - 1 + pos]
 z[orientation == -1, pos := super_end - pos + 1]
 z[, scaffold := NULL]
 setnames(z, "super", "scaffold")
 z[, c("super_start", "super_end", "orientation") := list(NULL, NULL, NULL)]
 super$info[, .(scaffold=super, scaffold_length=length)][z, on="scaffold"]->z
 z->s_cssaln

 if(molecules){
  copy(assembly$molecules)->z
  z[, c("orig_scaffold", "orig_start") := list(NULL, NULL)]
  super$agp[, .(scaffold, super, super_start, super_end, orientation)][z, on="scaffold"]->z
  z[orientation == 1, start := super_start - 1 + start]
  z[orientation == 1, end := super_start - 1 + end]
  z[orientation == -1, start := super_end - end + 1]
  z[orientation == -1, end := super_end - start + 1]
  z[, scaffold := NULL]
  z[, c("super_start", "super_end", "orientation") := list(NULL, NULL, NULL)]
  setnames(z, "super", "scaffold")
  z -> s_molecules
 } else {
  s_molecules <- NULL
 }

 if(!is.null(assembly$fpairs) && nrow(assembly$fpairs) > 0){
  copy(assembly$fpairs)->z
  z[, c("orig_scaffold1", "orig_pos1") := list(NULL, NULL)]
  z[, c("orig_scaffold2", "orig_pos2") := list(NULL, NULL)]
  z[, c("chr1", "chr2") := list(NULL, NULL)]
  super$agp[, .(scaffold1=scaffold, super, super_start, super_end, orientation)][z, on="scaffold1"]->z
  z[orientation == 1, pos1 := super_start - 1 + pos1]
  z[orientation == -1, pos1 := super_end - pos1 + 1]
  z[, scaffold1 := NULL]
  setnames(z, "super", "scaffold1")
  z[, c("super_start", "super_end", "orientation") := list(NULL, NULL, NULL)]
  super$agp[, .(scaffold2=scaffold, super, super_start, super_end, orientation)][z, on="scaffold2"]->z
  z[orientation == 1, pos2 := super_start - 1 + pos2]
  z[orientation == -1, pos2 := super_end - pos2 + 1]
  z[, scaffold2 := NULL]
  setnames(z, "super", "scaffold2")
  z[, c("super_start", "super_end", "orientation") := list(NULL, NULL, NULL)]
 } else {
  z <- NULL
 }

 init_assembly(fai=super$agp_bed[, .(length=sum(bed_end - bed_start)), key=.(scaffold=super)], cssaln=s_cssaln, molecules=s_molecules, fpairs=z) 
}