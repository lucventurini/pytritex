lift_nucfile<-function(assembly, map_10x, nucfile, outfile){
 n <- c("orig_scaffold", "orig_start", "orig_end", "frag_id", "nA", "nC", "nG", "nT", "nN", "length")
 nuc <- fread(nucfile, select=c(1:4,7:11,13), head=T, col.names=n)
 assembly$info[, .(scaffold, orig_scaffold, orig_start, off=orig_start)][nuc, on=c("orig_scaffold", "orig_start"), roll=T]->z
 z[, start := orig_start - off + 1]
 z[, end := orig_end - off + 1]
 map_10x$agp[gap == F, .(super, orientation, super_start, super_end, scaffold)][z, on="scaffold"]->z
 z[orientation == 1, start := super_start - 1 + start]
 z[orientation == 1, end := super_start - 1 + end]
 z[orientation == -1, start := super_end - end + 1]
 z[orientation == -1, end := super_end - start + 1]
 z[, .(super, start - 1, end, frag_id, ".", ".", nA, nC, nG, nT, nN, ".", length)]->zz
 fwrite(zz, file=outfile, col.names=T, row.names=F, sep="\t", quote=F)
}