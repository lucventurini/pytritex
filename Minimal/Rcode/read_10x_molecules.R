# Read the files as produced by run_10x_mapping.zsh
read_10x_molecules<-function(files, ncores=1){
 f <- files
 rbindlist(mclapply(mc.cores=ncores, names(f), function(i){
  fread(cmd=paste('gzip -cd', f[i]), head=F, col.names=c("scaffold", "start", "end", "barcode", "npairs"))->z
  z[, sample := i]
 })) -> mol
 mol[, start := mol$start + 1]
 mol[, length := mol$end - mol$start + 1]
 mol[]
}

# Initialized an assembly object with genetic map/flowsorting information, 10X and Hi-C links
init_assembly<-function(fai, cssaln, fpairs=NULL, molecules=NULL, rename=NULL){

 copy(fai)->info
 info[, orig_start := 1]
 info[, orig_end := length]
 if(is.null(rename)) {
  info[, orig_scaffold := scaffold]
 } else {
  setnames(info, "scaffold", "orig_scaffold")
  rename[info, on="orig_scaffold"] -> info
 }

 copy(cssaln)->z
 if(is.null(rename)) {
  z[, orig_scaffold := scaffold] 
 } else {
  setnames(z, "scaffold", "orig_scaffold")
  rename[z, on="orig_scaffold"] -> z
 }
 z[, orig_pos := pos]
 z[, orig_scaffold_length := scaffold_length]

 if(!is.null(molecules)){
  copy(molecules)->y
  y[, orig_scaffold := scaffold]
  y[, orig_start := start]
  y[, orig_end := end]
  y[, scaffold := NULL]
  info[, .(orig_scaffold, scaffold)][y, on="orig_scaffold"]->y
 } else {
  y <- data.table()
 }

 if(!is.null(fpairs)){
  copy(fpairs)->tcc
  tcc[, orig_scaffold1 := scaffold1]
  tcc[, orig_pos1 := pos1]
  tcc[, orig_scaffold2 := scaffold2]
  tcc[, orig_pos2 := pos2]
  tcc[, scaffold1 := NULL]
  tcc[, scaffold2 := NULL]
  info[, .(orig_scaffold1=orig_scaffold, scaffold1=scaffold)][tcc, on="orig_scaffold1"]->tcc
  info[, .(orig_scaffold2=orig_scaffold, scaffold2=scaffold)][tcc, on="orig_scaffold2"]->tcc
 } else {
  tcc <- data.table()
 }

 list(info=info, cssaln=z, fpairs=tcc, molecules=y)
}