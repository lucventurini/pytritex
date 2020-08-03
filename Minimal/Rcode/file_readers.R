library(data.table)
library(parallel)
library(igraph)
library(stringi)
library(zoo)

# Read BAM files with alignment of genetic markers (originally Wheat CSS contigs, the "marker sequences" of the POPSEQ map). Merge with POPSEQ data. Extract flow-sorting information from contig names/
# Specific to bread wheat and CSS POPSEQ map.
read_cssaln<-function(bam, popseq, fai, minqual=30, minlen=1000){
 cssaln <- fread(cmd=paste("samtools view -q ", minqual, "-F260", bam, "| cut -f 1,3,4"), col.names=c("css_contig", "scaffold", "pos")) 

 cssaln[, sorted_lib := sub("_[0-9]+$", "", sub("ta_contig_", "", css_contig))]
 cssaln[, sorted_alphachr := sub("L|S$", "", sorted_lib)]
 cssaln[, sorted_subgenome := sub("[1-7]", "", sorted_alphachr)]
 cssaln[, sorted_arm := sub("[1-7][A-D]", "", sorted_lib)]

 popseq[cssaln, on="css_contig"]->cssaln

 chrNames(species="wheat")->wheatchr
 setnames(copy(wheatchr), c("popseq_alphachr", "popseq_chr"))[cssaln, on="popseq_alphachr"] -> cssaln
 setnames(copy(wheatchr), c("sorted_alphachr", "sorted_chr"))[cssaln, on="sorted_alphachr"] -> cssaln
 fai[, .(scaffold, scaffold_length = length)][cssaln, on="scaffold"]->cssaln
 cssaln[css_contig_length >= minlen]
}

# Genetic function for reading minimap2 alignment of genetic maps and merging them genetic positional infromation. First used for barley. cv Morex, hence the name.
read_morexaln_minimap<-function(paf, popseq, minqual=30, minlen=500, prefix=T){
 fread(cmd=paste0("zgrep tp:A:P ", paf, " | awk -v l=", minlen," -v q=", minqual, " '$2 >= l && $12 >= q' | cut -f 1,6-8"),
              head=F, col.names=c("css_contig", "scaffold", "scaffold_length", "pos"))->z
 if(prefix){
  z[, css_contig := paste0("morex_", css_contig)]
 }
 popseq[z, on="css_contig"]
}

# Generic function to read PAF files
read_paf<-function(file, primary_only=T, save=F){
 fread(head=F, cmd=paste("zcat", file, "| cut -f -13"), 
       col.names=c("query", "query_length", "query_start", "query_end",
		   "orientation", "reference", "reference_length",
		   "reference_start", "reference_end", "matches", "alnlen", "mapq", "type"))->z
 z[, c("query_start", "query_end", "reference_start", "reference_end") := list(query_start + 1, query_end + 1, reference_start + 1, reference_end + 1)]
 z[, orientation := as.integer(ifelse(orientation == "+", 1, -1))]
 z[, type := sub("tp:A:", "", type)]
 if(primary_only){
  z["P", on="type"]->z
 }
 if(save){
  saveRDS(z, file=sub("paf.gz$", "Rds", file))
 }
 z[]
}

# Aggregate alignments for each (reference, query) pair
summarize_paf<-function(paf){
 setnames(dcast(paf[, .(l=sum(alnlen)), key=.(query, reference, orientation)], query + reference ~ orientation, value.var="l", fill=0), c("-1", "1"), c("lenrev", "lenfw"))->fr
 paf[, .(query_start=min(query_start), query_end=max(query_end),
       reference_start=min(reference_start), reference_end=max(reference_end),
       matches=sum(matches), alnlen=sum(alnlen), naln=.N), 
    key=.(query, query_length, reference, reference_length)]->paf_summary
 fr[paf_summary, on=c("query", "reference")]->paf_summary
 setorder(paf_summary, -alnlen)
 paf_summary[, idx := paste(1:.N)][]
}


# Read BED files with positions of restriction fragments on the input assembly, lift coordinates to updated assemblies
read_fragdata<-function(info, map_10x=NULL, assembly_10x=NULL, file){
 fragbed<-fread(file, head=F, col.names=c("orig_scaffold", "start", "end"))
 fragbed[, length := end - start]
 fragbed[, start := start + 1]
 info[, .(scaffold, start=orig_start, orig_start, orig_scaffold)][fragbed, on=c("orig_scaffold", "start"), roll=T]->fragbed
 fragbed[, start := start - orig_start + 1]
 fragbed[, end := end - orig_start + 1]
 fragbed[, orig_start := NULL]
 fragbed[, orig_scaffold := NULL]
 if(!is.null(assembly_10x)){
  map_10x$agp[gap == F, .(super, orientation, super_start, super_end, scaffold)][fragbed, on="scaffold"]->fragbed
  fragbed[orientation == 1, start := super_start - 1 + start]
  fragbed[orientation == 1, end := super_start - 1 + end]
  fragbed[orientation == -1, start := super_end - end + 1]
  fragbed[orientation == -1, end := super_end - start + 1]
  fragbed[, c("orientation", "super_start", "super_end", "scaffold") := list(NULL, NULL, NULL, NULL)]
  setnames(fragbed, "super", "orig_scaffold")

  assembly_10x$info[, .(scaffold, start=orig_start, orig_start, orig_scaffold)][fragbed, on=c("orig_scaffold", "start"), roll=T]->fragbed
  fragbed[, start := start - orig_start + 1]
  fragbed[, end := end - orig_start + 1]
  fragbed[, orig_start := NULL]
  fragbed[, orig_scaffold := NULL]

  fragbed[, .(nfrag = .N), keyby=scaffold][assembly_10x$info, on="scaffold"][is.na(nfrag), nfrag := 0]->z
 } else {
  fragbed[, .(nfrag = .N), keyby=scaffold][info, on="scaffold"][is.na(nfrag), nfrag := 0]->z
 }
 list(bed=fragbed[], info=z[])
}