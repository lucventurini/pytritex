import pandas as pd
import numpy as np
from .utils.chrnames import chrNames


def make_super_scaffolds(links, info: pd.DataFrame, excluded=(), ncores=1, prefix=None):
    info2 = info[["scaffold", "popseq_chr", "popseq_cM", "length"]][:].rename(columns={
        "popseq_chr": "chr", "popseq_cM": "cM"})
    excluded_scaffolds = excluded
    input_df = info2[info2.scaffold.isin(excluded_scaffolds)].rename(
        columns={"scaffold": "cluster"}
    )
    hl = links.copy()


make_super_scaffolds<-function(links, info, excluded=c(), ncores, prefix){
  info[, .(scaffold, chr=popseq_chr, cM=popseq_cM, length)]->info2
  excluded -> excluded_scaffolds
  info2[, excluded := scaffold %in% excluded_scaffolds] -> input
  setnames(input, "scaffold", "cluster")
  setnames(copy(links), c("scaffold1", "scaffold2"), c("cluster1", "cluster2"))->hl
  make_super(hl=hl, cluster_info=input, verbose=F, prefix=prefix, cores=ncores, paths=T, path_max=0, known_ends=F, maxiter=100)->s
  copy(s$membership) -> m
  setnames(m, "cluster", "scaffold")

  max(as.integer(sub(paste0(prefix, "_"), "", s$super_info$super))) -> maxidx
  rbind(m, info[!m$scaffold, on="scaffold"][, .(scaffold, bin=1, rank=0, backbone=T, chr=popseq_chr, cM=popseq_cM,
					length=length, excluded=scaffold %in% excluded_scaffolds, super=paste0(prefix, "_", maxidx + 1:.N))])->m

  m[, .(n=.N, nbin=max(bin), max_rank=max(rank), length=sum(length)), key=super]->res
  res[, .(super, super_size=n, super_nbin=nbin)][m, on="super"]->m
  list(membership=m, info=res)
 }




def anchor_scaffolds(assembly: dict,
                     popseq,
                     species=None,
                     sorted_percentile=95,
                     popseq_percentile=90,
                     hic_percentile=98):
    if species is None:
        raise KeyError(
            "Parameter 'species' is NULL. Please set 'species' to one of"
"\"wheat\", \"barley\", \"oats\", \"lolium\", \"sharonensis\" or \"rye\".")
    elif species not in ("wheat", "barley", "rye", "oats", "sharonensis", "lolium"):
        raise KeyError(
            "Parameter 'species' is not valid. Please set 'species' to one of"
"\"wheat\", \"barley\", \"oats\", \"lolium\", \"sharonensis\" or \"rye\"."
        )
    wheatchr = chrNames(species=species)
    wheatchr.columns = ("popseq_alphachr", "popseq_chr")
    fai = assembly["info"]
    cssaln = assembly["cssaln"]
    if "fpairs" not in assembly:
        fpairs = None
        hic = False
    else:
        fpairs = assembly["fpairs"]
        hic = (assembly["fpairs"].shape[0] > 0)
    z = cssaln[(cssaln.]
# Use POSPEQ, Hi-C and flow-sorting infromation to assign scaffolds to approximate chromosomal locations
anchor_scaffolds < -function(assembly, popseq, species=NULL,
                             sorted_percentile=95,
                             popseq_percentile=90,
                             hic_percentile=98)
{

assembly$fpairs -> fpairs
}

# Assignment of CARMA chromosome
cssaln[! is.na(sorted_alphachr),.N, keyby =.(scaffold, sorted_alphachr)]->z
z[order(-N)][,.(Ncss=sum(N), sorted_alphachr = sorted_alphachr[1], sorted_Ncss1 = N[1],
                                                                                  sorted_alphachr2 = sorted_alphachr[
                                                                                                         2], sorted_Ncss2 =
N[2]), keyby = scaffold]->z
z[, sorted_pchr := sorted_Ncss1 / Ncss]
z[, sorted_p12 := sorted_Ncss2 / sorted_Ncss1]

# Assignment of CARMA chromosome arm
cssaln[sorted_arm == "L",.(NL=.N), keyby =.(scaffold, sorted_alphachr)]->al
cssaln[sorted_arm == "S",.(NS=.N), keyby =.(scaffold, sorted_alphachr)]->as
al[z, on = c("scaffold", "sorted_alphachr")]->z
as[z, on = c("scaffold", "sorted_alphachr")]->z
z[ is.na(NL), NL := 0]
z[ is.na(NS), NS := 0]
z[, sorted_arm := ifelse(NS >= NL, "S", "L")]
z[sorted_alphachr % in % c("1H", "3B"), sorted_arm := NA]
z[sorted_arm == "S", sorted_parm := NS / sorted_Ncss1]
z[sorted_arm == "L", sorted_parm := NL / sorted_Ncss1]
setnames(copy(wheatchr), c("sorted_alphachr", "sorted_chr"))[z, on = "sorted_alphachr"]->z
setnames(copy(wheatchr), c("sorted_alphachr2", "sorted_chr2"))[z, on = "sorted_alphachr2"]->z
z[fai, on = "scaffold"]->info
info[ is.na(Ncss), Ncss := 0]
info[ is.na(NS), NS := 0]
info[ is.na(NL), NL := 0]
info[ is.na(sorted_Ncss1), sorted_Ncss1 := 0]
info[ is.na(sorted_Ncss2), sorted_Ncss2 := 0]

# Assignment of POPSEQ genetic positions
popseq[! is.na(popseq_alphachr),.(css_contig, popseq_alphachr, popseq_cM)][cssaln[,.(css_contig,
                                                                                     scaffold)], on = "css_contig", nomatch = 0]->z
z[,.(.N, popseq_cM=mean(popseq_cM), popseq_cM_sd = ifelse(length(popseq_cM) > 1, sd(popseq_cM), 0), popseq_cM_mad = mad(
    popseq_cM)), keyby =.(scaffold, popseq_alphachr)]->zz
zz[, popseq_Ncss := sum(N), by = scaffold]->zz
zz[order(-N)][,.(popseq_alphachr=popseq_alphachr[1], popseq_Ncss1=N[1],
popseq_alphachr2=popseq_alphachr[2], popseq_Ncss2=N[2]), keyby = scaffold]->x
zz[,.(scaffold, popseq_alphachr, popseq_Ncss, popseq_cM, popseq_cM_sd, popseq_cM_mad)][x, on = c("scaffold",
                                                                                                 "popseq_alphachr")]->zz
zz[, popseq_pchr := popseq_Ncss1 / popseq_Ncss]
zz[, popseq_p12 := popseq_Ncss2 / popseq_Ncss1]
wheatchr[zz, on = "popseq_alphachr"]->zz
setnames(copy(wheatchr), paste0(names(wheatchr), 2))[zz, on = "popseq_alphachr2"]->zz
zz[info, on = "scaffold"]->info
info[ is.na(popseq_Ncss), popseq_Ncss := 0]
info[ is.na(popseq_Ncss1), popseq_Ncss1 := 0]
info[ is.na(popseq_Ncss2), popseq_Ncss2 := 0]

# Chromosome assignment based on Hi-C links
if (hic){
info[!(popseq_chr != sorted_chr)][, .(scaffold, chr=popseq_chr)]->info0
setnames(copy(info0), names(info0), sub("$", "1", names(info0)))[fpairs, on="scaffold1"]->tcc_pos
setnames(copy(info0), names(info0), sub("$", "2", names(info0)))[tcc_pos, on="scaffold2"]->tcc_pos
tcc_pos[! is.na(chr1), .N, key=.(scaffold=scaffold2, hic_chr=chr1)]->z
z[order(-N)][, .(Nhic=sum(N), hic_chr=hic_chr[1], hic_N1=N[1],
hic_chr2=hic_chr[2], hic_N2=N[2]), keyby=scaffold]->zz
zz[, hic_pchr := hic_N1 / Nhic]
zz[, hic_p12 := hic_N2 / hic_N1]
zz[info, on="scaffold"]->info
info[is.na(Nhic), Nhic := 0]
info[is.na(hic_N1), hic_N1 := 0]
info[is.na(hic_N2), hic_N2 := 0]
}

if (hic){
measure < - c("popseq_chr", "hic_chr", "sorted_chr")
} else {
measure < - c("popseq_chr", "sorted_chr")
}
melt(info, id.vars="scaffold", measure.vars=measure, variable.factor=F, variable.name="map", na.rm=T, value.name="chr")->w
w[, .N, key=.(scaffold, chr)]->w
w[order(-N), .(Nchr_ass = sum(N), Nchr_ass_uniq =.N), keyby=scaffold]->w
w[info, on="scaffold"]->info
info[is.na(Nchr_ass), Nchr_ass := 0]
info[is.na(Nchr_ass_uniq), Nchr_ass_uniq := 0]

x < -info[Ncss >= 30, quantile(na.omit(sorted_p12), 0:100 / 100)][sorted_percentile + 1]
info[, bad_sorted := (sorted_p12 >= x & sorted_Ncss2 >= 2)]
x < -info[popseq_Ncss >= 30, quantile(na.omit(popseq_p12), 0:100 / 100)][popseq_percentile + 1]
info[, bad_popseq := (popseq_p12 >= x & popseq_Ncss2 >= 2)]

info[ is.na(bad_sorted), bad_sorted := F]
info[ is.na(bad_popseq), bad_popseq := F]

if (hic){
x < -info[Nhic >= 30, quantile(na.omit(hic_p12), 0:100 / 100)][hic_percentile + 1]
info[Nhic >= 30, bad_hic := hic_p12 >= x & hic_N2 >= 2]
info[ is.na(bad_hic), bad_hic := F]
}

melt(info, id.vars = "scaffold", measure.vars = grep(value=T, "bad_", names(
    info)), variable.factor = F, variable.name = "map", na.rm = T, value.name = "bad")[bad == T]->w
w[,.(Nbad=.N), key = scaffold]->w
w[info, on = "scaffold"]->info
info[ is.na(Nbad), Nbad := 0]

assembly$info < - info
assembly$popseq < - popseq
if (hic)
{
    assembly$fpairs < - tcc_pos
}
assembly
}