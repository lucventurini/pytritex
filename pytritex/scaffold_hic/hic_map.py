import dask.dataframe as dd
import pandas as pd
import numpy as np
import os
from typing import Union
from .hic_map_constructor import make_hic_map
from dask.distributed import Client


def hic_map(assembly: dict, client: Client,
            fragment_data, species, ncores=1, min_nfrag_scaffold=50, max_cM_dist = 20,
            min_length=3e5, binsize=5e5, min_nfrag_bin=30, gap_size=100, maxiter=100, orient=True, agp_only=False,
            map=None, known_ends=True, min_binsize=1e5, min_nbin=5):

    #   copy(info)->hic_info
    #   hic_info[, excluded := nfrag < min_nfrag_scaffold]
    #
    #   assembly$fpairs[scaffold1 != scaffold2, .(nlinks=.N), key=.(scaffold1, scaffold2)]->hl
    #   hic_info[, .(scaffold1=scaffold, chr1=chr, cM1=cM)][hl, nomatch=0, on="scaffold1"]->hl
    #   hic_info[, .(scaffold2=scaffold, chr2=chr, cM2=cM)][hl, nomatch=0, on="scaffold2"]->hl
    #   hl[chr1 == chr2]->hl
    #   hl<-hl[abs(cM1-cM2) <= max_cM_dist | is.na(cM1) | is.na(cM2)]
    #   hl[, weight:=-log10(nlinks)]
    hic_info = dd.read_parquet(fragment_data["info"], infer_divisions=True)
    hic_info = hic_info.query("hic_chr == hic_chr & super_length >= @min_length",
                              local_dict={"min_length": min_length})[["nfrag", "hic_chr", "popseq_cM"]]
    hic_info["excluded"] = hic_info.eval("nfrag < @mnf", local_dict={"mnf": min_nfrag_scaffold})
    hl = assembly["fpairs"].eval("scaffold_index1 != scaffold_index2").persist()
    hl = hl.merge(hl.groupby(["scaffold_index1", "scaffold_index2"]).size().to_frame("nlinks").compute(),
                  on=["scaffold_index1", "scaffold_index2"]).persist()

    left = hic_info[["chr", "cM"]]
    left1 = left.rename(columns={"chr": "chr1", "cM": "cM1"})
    left1.index = left1.index.rename("scaffold_index1")
    left2 = hic_info.rename(columns={"chr": "chr2", "cM": "cM2"})
    left2.index = left1.index.rename("scaffold_index2")
    hl = hl.merge(left1, on="scaffold_index1", how="inner").merge(left2, on="scaffold_index2", how="inner").query(
        "chr1 == chr2")
    hl["cMDist"] = hl.eval("cMDist = cM1 - cM2")
    hl = hl.query("cMDist <= @max_cM_dist | cMDist != cMDist", local_dict={"max_cM_dist": max_cM_dist})
    # TODO: why is the log minused? Should it not be plus?
    hl["weight"] = -np.log10(hl["nlinks"])
    # The following is already implemented
    #   cat("Scaffold map construction started.\n")
    #   make_hic_map(hic_info=hic_info, links=hl, ncores=ncores, known_ends=known_ends)->hic_map_bin
    #   cat("Scaffold map construction finished.\n")
    hic_map_bin = make_hic_map(hic_info=hic_info, links=hl, client=client, ncores=ncores, known_ends=known_ends)

    hic_map_oriented = orient_hic_map()

    # Now determine the orientation?

    #     w<-hic_map_bin[!is.na(hic_bin), .(id=scaffold, scaffold=sub(":.*$", "", scaffold), pos=as.integer(sub("^.*:", "", scaffold)), chr, hic_bin)]
    #     w<-w[, .(gbin=mean(na.omit(hic_bin)),
    # 	     hic_cor=as.numeric(suppressWarnings(cor(method='s', hic_bin, pos, use='p')))), keyby=scaffold][!is.na(hic_cor)]
    #     hic_map[!is.na(hic_bin) & scaffold %in% w$scaffold][order(chr, hic_bin)]->z0

    fai = dd.read_parquet(assembly["fai"], infer_divisions=True)

    left = hic_map_bin.query("hic_bin != hic_bin")[["chr", "hic_bin"]]
    left = left.merge(fai[["orig_scaffold_index", "orig_start"]], on="scaffold_index")
    left["gbin"] =



    #     z0[,.(scaffold1=scaffold[1:(.N-2)], scaffold2=scaffold[2:(.N-1)], scaffold3=scaffold[3:(.N)]), by=chr]->z
    #     z0[, data.table(key="scaffold1", scaffold1=scaffold, hic_bin1=hic_bin)][setkey(z, "scaffold1")]->z
    #     z0[, data.table(key="scaffold2", scaffold2=scaffold, hic_bin2=hic_bin)][setkey(z, "scaffold2")]->z
    #     z0[, data.table(key="scaffold3", scaffold3=scaffold, hic_bin3=hic_bin)][setkey(z, "scaffold3")]->z
    #     w[, data.table(key="scaffold1", scaffold1=scaffold, gbin1=gbin)][setkey(z, "scaffold1")]->z
    #     w[, data.table(key="scaffold2", scaffold2=scaffold, gbin2=gbin)][setkey(z, "scaffold2")]->z
    #     w[, data.table(key="scaffold3", scaffold3=scaffold, gbin3=gbin)][setkey(z, "scaffold3")]->z
    #     z[, cc:= apply(z[, .(hic_bin1, hic_bin2, hic_bin3, gbin1, gbin2, gbin3)],1,function(x) {
    # 		    suppressWarnings(cor(x[1:3], x[4:6]))
    # 		    })]
    #     z[, data.table(key="scaffold", scaffold=scaffold2, cc=ifelse(cc > 0, 1, -1))]->ccor
    #     ccor[w]->m
    #     m[, hic_orientation:=ifelse(hic_cor > 0, 1 * cc, -1 * cc)]
    #     m[, .(scaffold, hic_cor, hic_invert=cc, hic_orientation)][hic_map, on="scaffold"]->hic_map_oriented
    #
    #     setnames(hic_map_oriented, "chr", "consensus_chr")
    #     setnames(hic_map_oriented, "cM", "consensus_cM")
    #     hic_map_oriented[, consensus_orientation := hic_orientation]



    dd.to_parquet()

    #     assembly$info[, .(scaffold, binsize=pmax(min_binsize, length %/% min_nbin))][frags, on='scaffold']->f
    #     f[, .(nfrag=.N), keyby=.(scaffold, binsize, pos = start %/% binsize * binsize)]->fragbin
    #     fragbin[, id := paste(sep=":", scaffold, pos)]
    #     fragbin<- hic_map[, .(scaffold, chr, cM=hic_bin)][fragbin, on="scaffold", nomatch=0]
    fragged = assembly["info"].assign(
        binsize=np.maximum(min_binsize, assembly["info"]["length"] // min_nbin)
    )[["binsize"]].merge(frags, on="scaffold_index")
    fragbin = fragged.assign(
        pos=(fragged["start"] // binsize) * binsize).groupby(["scaffold_index", "binsize", "pos"]).size().toframe("nfrag")
    fragbin = hic_map[["chr", "hic_bin"]].rename(columns={"hic_bin": "cM"}).merge(fragbin, on="scaffold_index",
                                                                                  how="inner")
    #     unique(fragbin[, .(scaffold1=scaffold, binsize1=binsize)])[assembly$fpairs, on='scaffold1']->fp
    #     unique(fragbin[, .(scaffold2=scaffold, binsize2=binsize)])[fp, on='scaffold2']->fp
    #     fp[, .(nlinks=.N), keyby=.(scaffold1, pos1 = pos1 %/% binsize1 * binsize1, scaffold2, pos2 = pos2 %/% binsize2 * binsize2)]->binl
    #     binl[, id1 := paste(sep=":", scaffold1, pos1)]
    #     binl[, id2 := paste(sep=":", scaffold2, pos2)]
    #     binl[id1 != id2]->binl

    left = fragbin[["binsize"]].reset_index(drop=False).rename(columns={"binsize": "binsize1",
                                                                        "scaffold_index": "scaffold_index1"}).drop_duplicates()
    fp = left.merge(assembly["fpairs"], on="scaffold_index1")
    if fp.index.name == "scaffold_index1":
        fp = fp.reset_index(drop=False)
    left = fragbin[["binsize"]].reset_index(drop=False).rename(columns={"binsize": "binsize2",
                                                                        "scaffold_index": "scaffold_index2"}).drop_duplicates()
    fp = left.merge(fp, on="scaffold_index2")
    if fp.index.name == "scaffold_index2":
        fp = fp.reset_index(drop=False)
    # Calculate number of links per bin?
    binl = fp.assign(pos1=(fp["pos1"] / fp["binsize1"]) * fp["binsize1"],
                     pos2=(fp["pos2"] / fp["binsize2"]) * fp["binsize2"]).groupby(
        ["scaffold_index1", "pos1", "scaffold_index2", "pos2"]).size().to_frame("nlinks")
    # Ignore self-links
    binl = binl.query("scaffold_index1 != scaffold_index2 | pos1 != pos2")
    # Merge with fragbin
    binl = fragbin.rename(columns={"scaffold_index": "scaffold_index1", "pos": "pos1", "chr": "chr1", "cM": "cM1"})[
        ["scaffold_index1", "chr1", "cM1"]].merge(binl, on=["scaffold_index1", "pos1"])
    binl = fragbin.rename(columns={"scaffold_index": "scaffold_index2", "pos": "pos2", "chr": "chr2", "cM": "cM2"})[
        ["scaffold_index2", "chr2", "cM2"]].merge(binl, on=["scaffold_index2", "pos2"])
    #     fragbin[, .(id1=id, chr1=chr, cM1=cM)][binl, on="id1"]->binl
    #     fragbin[, .(id2=id, chr2=chr, cM2=cM)][binl, on="id2"]->binl
    #     binl[, c("scaffold1", "scaffold2", "pos1", "pos2") := list(NULL, NULL, NULL, NULL)]
    #     setnames(binl, c("id1", "id2"), c("scaffold1", "scaffold2"))

    # fragbin[, .(scaffold=id, nfrag, chr, cM)]->hic_info_bin
    #     hic_info_bin[, excluded:=nfrag < min_nfrag_bin]
    #     binl[chr1 == chr2 & (abs(cM1-cM2) <= 2 | is.na(cM1) | is.na(cM2))]->binl
    #     binl[, weight:=-log10(nlinks)]
    #
    #     make_hic_map(hic_info=hic_info_bin, links=binl, ncores=ncores, maxiter=maxiter, known_ends=known_ends)->hic_map_bin

    hic_info_bin = fragbin[["scaffold_index", "pos", "nfrag", "chr", "cM"]]
    hic_info_bin["excluded"] = hic_info_bin["nfrag"] < min_nfrag_bin
    # TODO Why 2 centimorgans max?
    binl = binl.query("chr1 == chr2 & (cM1 != cM1 | cM2 != cM2 | (-2 <= cM1 -cM2 <= 2))")[:]
    # TODO why is the weight negatively correlated with the number of links?
    binl["weight"] = -np.log10(binl["nlinks"])
    # TODO why are we recalculating the map?
    hic_map_bin = make_hic_map(hic_info=hic_info_bin, links=binl, client=client, ncores=ncores, known_ends=known_ends)
