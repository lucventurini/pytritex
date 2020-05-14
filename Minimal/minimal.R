#source TRITEX R code 
source('../pseudomolecule_construction.R')
options(error=function() { traceback(2); if(!interactive()) quit("no", status = 1, runLast = FALSE) })

#wheat POPSEQ map
readRDS(file="../180130_wheat_CSS_POPSEQ.Rds") -> popseq

#read contig lengths
f <- 'minimal.fa.fai'
fread(f, head=F, select=1:2, col.names=c("scaffold", "length")) -> fai

#read alignments of Chinese Spring CSS cpontigs
f <- 'popseq_eiref.paf.gz'
read_morexaln_minimap(paf=f, popseq=popseq, minqual=30, minlen=500, prefix=F)->cssaln

#read 10X alignments
dirtenex <- '.'
# sort(fread(paste("find", dirtenex, "-type f | sort | grep 'molecules.tsv.gz$'"), head=F)$V1)->f
names <- c("LIB35275_190111", "LIB35275_190306", 
	      "LIB35276_190111", "LIB35276_190306",
	      "LIB35277_190111", "LIB35277_190306",
	      "LIB35278_190111", "LIB35278_190306");
paster <- function(name){paste(
"10X", paste("Sample_PRO1373_S20", name, sep="_"), paste("Sample_PRO1373_S20", name, "molecules.tsv.gz", sep="_"), sep="/")}
lapply(names, paster) -> f;
names(f) <- names;
read_10x_molecules(files=f, ncores=8) -> molecules

#read Hi-C alignments
dirhic <- '.' 
fread(paste('/bin/ls HiC/HiC_fragment_pairs.tsv.gz | xargs zcat'),
  header=F, col.names=c("scaffold1", "pos1", "scaffold2", "pos2"))->fpairs

#initialize assembly object
init_assembly(fai=fai, cssaln=cssaln, molecules=molecules, fpairs=fpairs) -> assembly
saveRDS(assembly, file="assembly_init.Rds")
anchor_scaffolds(assembly = assembly, popseq=popseq, species="wheat") -> assembly
add_molecule_cov(assembly = assembly, binsize=200, cores=10) -> assembly
add_hic_cov(assembly, binsize=5e3, binsize2=5e4, minNbin=50, innerDist=3e5, cores=10)->assembly
saveRDS(assembly, file="alchemy_v0.Rds")

# #find potential breaks and plot chimeras
find_10x_breaks(assembly=assembly, interval=5e4, minNbin=20, dist=1e4, ratio=-3) -> breaks
saveRDS(breaks, file="alchemy_v0_breaks.Rds")
#breaks[d >= 1e4][order(-d)][1:100]->b
#plot_chimeras(assembly=assembly, scaffolds=b, breaks=b, species="wheat",
# 		   file="20200414_assembly_v0_chimeras.pdf", cores=1)

# #break chimeras
break_10x(assembly,
          prefix="scaffold_corrected",
          ratio=-3,
          interval=5e4,
          minNbin=20,
          dist=2e3,
          slop=2e2,
          species="wheat", intermediate=F, ncores=8)->a
cat("Broken chimeras\n")

# #save corrected assembly
a$assembly->assembly_v1
saveRDS(assembly_v1, file="20200414_assembly_v1.Rds")
saveRDS(a$breaks, file="20200414_assembly_breaks.Rds")

quit()
# #run 10X scaffolding for a grid of  parameters
# #npairs : minimum # of aligned read pairs per molecules
# #nmol : minimum # molecules to accept a scaffold join
# #nsample : no. of independent 10X sample supporting a link
# #dist : maximum distance of molecules from scaffold end, smaller values are more stringent
setnames(data.table(expand.grid(2:3, 2:3, 1:2, 6:9*1e4)), c("npairs", "nmol", "nsample", "dist"))->grid
cat("Starting grid evaluation\n")
scaffolder<-function(i, grid, assembly_v1){
grid[i, npairs]->n
 grid[i, nmol]->m
 grid[i, nsample]->s
 grid[i, dist]->d
 scaffold_10x(assembly=assembly_v1, prefix="scaffold_10x", 
	      min_npairs=n, max_dist=d, min_nmol=m, popseq_dist=5, 
	      max_dist_orientation=5, min_nsample=s, unanchored=F, ncores=1)->z
 cat(paste(n, " ", m, " ", s, " ", n50(z$info$length), "\n"))
 z
}

cat(paste("Assembly v1:", nrow(assembly_v1), ncol(assembly_v1), colnames(assembly_v1), "\n"))

lapply(1:nrow(grid), scaffolder, grid=grid, assembly_v1=assembly_v1)->res

# #exclude failed trials, sort results by n50
which(unlist(lapply(res, length)) ==  4) -> idx

res[idx]->res2
grid[idx]->grid2
n50_list <- lapply(res2, function(i) n50(i$info$length))
n50_u <- as.numeric(n50_list);

data.table(grid2, index=1:length(res2), n50=n50_u)[order(-n50)][1:10]
cat("Excluded failed trials\n")

# #pick best value and save
cat(paste(nrow(res2), ncol(res2), "\n"))
res2[[9]] -> assembly_v1_10x 
saveRDS(assembly_v1_10x, file="20200414_assembly_v1_10x.Rds")

# #init super scaffold assembly object
# init_10x_assembly(assembly=assembly_v1, map_10x=assembly_v1_10x, molecules=F)->assembly_v2
# anchor_scaffolds(assembly = assembly_v2, popseq=popseq, species="wheat") -> assembly_v2
# add_hic_cov(assembly_v2, binsize=1e4, binsize2=1e6, minNbin=100, innerDist=3e5, cores=16)->assembly_v2
# saveRDS(assembly_v2, file="20200414_assembly_v2.Rds")

# # read DpnII fragment data
# f <- 'minimal_DpnII_fragments_30bp.bed'
# read_fragdata(info=assembly_v1$info, file=f, map_10x=assembly_v1_10x, assembly_10x=assembly_v2)->frag_data

# # exclude scaffolds <= 300 kb from Hi-C map construction
# frag_data$info[!is.na(hic_chr) & length >= 3e5, .(scaffold, nfrag, chr=hic_chr, cM=popseq_cM)]->hic_info

# # make Hi-C map
# hic_map(info=hic_info, assembly=assembly_v2, frags=frag_data$bed, species="wheat", ncores=8,
# 	min_nfrag_scaffold=50, max_cM_dist = 50,
# 	binsize=1e5, min_nfrag_bin=20, gap_size=100)->hic_map_v1

# # add links data and compute contact matrices
# f <- 'minimal_DpnII_fragments_30bp_split.nuc.txt'
# add_psmol_fpairs(assembly=assembly_v1, hic_map=hic_map_v1, map_10x=assembly_v1_10x, 
# 		 assembly_10x=assembly_v2, nucfile=f)->hic_map_v1

# bin_hic_step(hic=hic_map_v1$links, frags=hic_map_v1$frags, binsize=1e6,
# 	    chrlen=hic_map_v1$chrlen, chrs=1:21, cores=21)->hic_map_v1$hic_1Mb
# normalize_cis(hic_map_v1$hic_1Mb, ncores=8, percentile=0, omit_smallest=1)->hic_map_v1$hic_1Mb$norm

# find_inversions(hic_map=hic_map_v1, links=hic_map_v1$hic_1Mb$norm, species="wheat", cores=30)->inv

# hic_cov_psmol(hic_map=hic_map_v1, binsize=1e3, binsize2=1e6, maxdist=1e6, cores=40) -> cov

# # find potential breaks and make a summary plot showing contact matrix, inversion statistic and directionality bias
# setorder(cov, r)[r < -2][, .(r=r[1]), key=.(chr, bin %/% 1e5 * 1e5)]->br0
# hic_map_v1$agp[, .(chr, scaffold, scaffold_length=agp_end-agp_start+1, agp_start, agp_end, orientation, bin=agp_start)][br0, on=c("chr", "bin"), roll=T]->br 
# br[orientation == 1, br := bin - agp_start]
# br[orientation == -1, br := agp_end - bin]
# hic_map_v1$chrlen[,.(chr, length)][br, on="chr"]->br
# br[, d := pmin(bin, length - bin)]
# br[, ds := pmin(br, scaffold_length - br)]
# copy(br) -> br1
# br1[ds >= 1e3 & d >= 1e5]->br

# big_hic_plot(hic_map=hic_map_v1, cov=cov, inv=inv, species="wheat", cores=21,
#   file="20200414_hic_map_v1_plot.pdf", breaks=br)
