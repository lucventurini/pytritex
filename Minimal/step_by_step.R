# Load the R functions for pseudomolecule construction.
source("/software/testing/tritex/src/tritexassembly.bitbucket.io/R/pseudomolecule_construction.R")
# source("./pseudomolecule_construction.R")

# Import a genetic map of sufficient density. In barley, we use the POPSEQ map.
readRDS('../../180130_wheat_CSS_POPSEQ.Rds') -> popseq
setnames(popseq, sub("morex", "css", names(popseq)))

# Import the table of scaffold lengths.
f <- 'minimal.fa.fai'
fread(f, head=F, select=1:2, col.names=c("scaffold", "length")) -> fai

# Alignment of genetic markers used for the genetic map. In this example, the Morex WGS assembly by IBSC (2012).
f <- 'minimal.fa.gz_csCSS.paf.gz'
read_morexaln_minimap(paf=f, popseq=popseq, minqual=30, minlen=500, prefix=T)->morexaln

# Read the list of Hi-C links.
dirhic <- './HiC'
fread(cmd=paste('find', dirhic, '| grep "_fragment_pairs.tsv.gz$" | xargs zcat'),
  header=F, col.names=c("scaffold1", "pos1", "scaffold2", "pos2"))->fpairs


# Read the list of 10X molecules.
dirtenex <- './10X/'
fread(cmd=paste("find", dirtenex, "-type f | grep 'molecules.tsv.gz$'"), head=F)$V1->f
# names(f) <- c("S1", "S2")
read_10x_molecules(files=f, ncores=length(f)) -> molecules

init_assembly(fai=fai, cssaln=morexaln, molecules=molecules, fpairs=fpairs) -> assembly
anchor_scaffolds(assembly = assembly, popseq=popseq, species="wheat") -> assembly
add_molecule_cov(assembly = assembly, cores=30) -> assembly
add_hic_cov(assembly, cores=30)->assembly
saveRDS(assembly, file="assembly_v1.Rds") 