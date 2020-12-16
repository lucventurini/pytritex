import dask.dataframe as dd
import pandas as pd
import numpy as np
import os
from typing import Union
from .hic_map_constructor import make_hic_map
from .orient_hic_map import orient_hic_map
from dask.distributed import Client
from .make_agp import make_agp


def _merge_with_hic_info(hl, hic_info):
    if "chr1" in hl.columns:
        hl = hl.drop(["chr1", "chr2"], axis=1)
    left = hic_info[["chr", "cM", "excluded"]]
    left1 = left.rename(columns=dict((_, _ + "1") for _ in left.columns))
    left1.index = left1.index.rename("scaffold_index1")
    left2 = left.rename(columns=dict((_, _ + "2") for _ in left.columns))
    left2.index = left1.index.rename("scaffold_index2")
    idx_name = hl.index.name
    assert idx_name is not None
    hl = hl.reset_index(drop=False).merge(
        left1, on="scaffold_index1", how="inner").merge(
        left2, on="scaffold_index2", how="inner").query("chr1 == chr2").eval("cMDist = cM1 - cM2").set_index(idx_name)
    hl["cMDist"] = hl["cMDist"].abs()
    return hl


def calculate_hic_link_weights(assembly:dict, save_dir: Union[None, str]):
    fpairs = assembly["fpairs"]
    
    if isinstance(fpairs, str):
        fpairs = dd.read_parquet(fpairs, infer_divisions=True)

    columns = ["scaffold_index1", "scaffold_index2", "pos1", "pos2"]
    hl = fpairs.query("scaffold_index1 != scaffold_index2")[columns]
    # Now we have to equalise the counts. Notice how we are reversing the order
    top = hl.query("scaffold_index1 < scaffold_index2")
    bottom = hl.query("scaffold_index1 > scaffold_index2")
    bottom = bottom.rename(columns={"scaffold_index1": "scaffold_index2a", "scaffold_index2": "scaffold_index1a",
                                    "pos1": "pos2a", "pos2": "pos1a"})
    merged = top.merge(bottom, left_on=columns, right_on=[_ + "a" for _ in columns], how="outer") 
    bait = merged["scaffold_index1"].isna()
    for column in ["scaffold_index1", "scaffold_index2", "pos1", "pos2"]:
        merged[column] = merged[column].mask(bait, merged[column + "a"])
    merged = merged.drop([_ + "a" for _ in columns], axis=1)
    nlinks = merged.groupby(["scaffold_index1", "scaffold_index2"]).size().to_frame("nlinks").reset_index(drop=False)
    nlinks["hic_index"] = 1
    nlinks["hic_index"] = nlinks["hic_index"].cumsum()
    nlinks = nlinks.set_index("hic_index", sorted=True)
    nlinks["weight"] = 1 + np.log10(nlinks["nlinks"])
    if save_dir is not None:
        link_folder = os.path.join(save_dir, "weighted_links")
        dd.to_parquet(nlinks, link_folder, compression="gzip")
        assembly["weighted_links"] = link_folder
        link_folder = os.path.join(save_dir, "nr_fpairs")
        dd.to_parquet(merged, link_folder, compression="gzip")
        assembly["nr_fpairs"] = link_folder
    else:
        assembly["weighted_links"] = nlinks
        assembly["nr_fpairs"] = merged
    return assembly
    

def hic_map(assembly: dict, client: Client,
            fragment_data, species, ncores=1, min_nlinks=10, min_nfrag_scaffold=50, max_cM_dist = 20,
            save_dir=None,
            min_length=3e5, binsize=5e5, min_nfrag_bin=30, gap_size=100, maxiter=100, orient=True, agp_only=False,
            map=None, known_ends=True, min_binsize=1e5, min_nbin=5):

    #   copy(info)->hic_info
    #   hic_info[, excluded := nfrag < min_nfrag_scaffold]
    #   assembly$fpairs[scaffold1 != scaffold2, .(nlinks=.N), key=.(scaffold1, scaffold2)]->hl
    #   hic_info[, .(scaffold1=scaffold, chr1=chr, cM1=cM)][hl, nomatch=0, on="scaffold1"]->hl
    #   hic_info[, .(scaffold2=scaffold, chr2=chr, cM2=cM)][hl, nomatch=0, on="scaffold2"]->hl
    #   hl[chr1 == chr2]->hl
    #   hl<-hl[abs(cM1-cM2) <= max_cM_dist | is.na(cM1) | is.na(cM2)]
    #   hl[, weight:=-log10(nlinks)]
    frag_info = dd.read_parquet(fragment_data["info"], infer_divisions=True)
    hic_info = dd.read_parquet(assembly["info"], infer_divisions=True)
    hic_info = hic_info.merge(frag_info[["nfrag"]], on="scaffold_index")
    hic_info["nfrag"] = hic_info["nfrag"].fillna(0)

    hic_info = hic_info.rename(columns={"popseq_chr": "chr", "popseq_cM": "cM"}).query(
        "chr == chr & length >= @min_length", local_dict={"min_length": min_length})[["chr", "cM", "nfrag", "length"]]
    hic_info["excluded"] = hic_info.eval("nfrag < @mnf", local_dict={"mnf": min_nfrag_scaffold})
    hl = assembly["weighted_links"]
    if isinstance(hl, str):
        hl = dd.read_parquet(hl, infer_divisions=True)

    hl = hl.query("nlinks >= @min_nlinks", local_dict={"min_nlinks": min_nlinks})

    hl = _merge_with_hic_info(hl, hic_info)
    hl = hl.query("excluded1 == False & excluded2 == False & (cMDist <= @max_cM_dist | cMDist != cMDist)",
                  local_dict={"max_cM_dist": max_cM_dist})

    # The following is already implemented
    #   cat("Scaffold map construction started.\n")
    #   make_hic_map(hic_info=hic_info, links=hl, ncores=ncores, known_ends=known_ends)->hic_map_bin
    #   cat("Scaffold map construction finished.\n")
    hic_map_bin = make_hic_map(hic_info=hic_info, links=hl, client=client, ncores=ncores, known_ends=known_ends,
                               save_dir=save_dir)

    # Now orient the HiC map
    hic_map_oriented = orient_hic_map(info=assembly["info"], assembly=assembly, hic_map=hic_map_bin,
                                      frags=fragment_data["info"], client=client, min_nfrag_bin=min_nfrag_bin,
                                      cores=ncores,
                                      maxiter=maxiter, orient_old=False, min_nbin=min_nbin, min_binsize=min_binsize)
    if save_dir is not None:
        dd.to_parquet(hic_map_oriented, os.path.join(save_dir, "hic_map_oriented"), compression="gzip")
        return os.path.join(save_dir, "hic_map_oriented")
    return hic_map_oriented


    # make_agp(hic_map_oriented, gap_size=gap_size, species=species)->a
    #
    #  a$agp[, .(length=sum(scaffold_length)), key=agp_chr]->chrlen
    #  chrlen[, alphachr := sub("chr", "", agp_chr)]
    #  chrNames(species=species)[chrlen, on="alphachr"]->chrlen
    #  chrlen[, truechr := !grepl("Un", alphachr)]
    #  chrlen[order(!truechr, chr)]->chrlen
    #  chrlen[, offset := cumsum(c(0, length[1:(.N-1)]))]
    #  chrlen[, plot_offset := cumsum(c(0, length[1:(.N-1)]+1e8))]
    #
    #  list(agp=a$agp, agp_bed=a$agp_bed, chrlen=chrlen, hic_map=hic_map_oriented, hic_map_bin=hic_map_bin)->res
    #  invisible(lapply(sort(c("min_nfrag_scaffold", "max_cM_dist", "binsize", "min_nfrag_bin", "gap_size")), function(i){
    #   res[[i]] <<- get(i)
    #  }))
    #  res

    # agp = make_agp(hic_map_oriented, gap_size=gap_size)


    # if("orientation" %in% names(hic_info)){
    #    hic_info[, .(scaffold, old_orientation=orientation)][hic_map_oriented, on="scaffold"]->hic_map_oriented
    #    hic_map_oriented[!is.na(old_orientation), consensus_orientation := old_orientation]
    #    hic_map_oriented[, old_orientation := NULL]
    #   }
    #
    #   hic_map_oriented[assembly$info, on="scaffold"]->hic_map_oriented

    #     assembly$info[, .(scaffold, binsize=pmax(min_binsize, length %/% min_nbin))][frags, on='scaffold']->f
    #     f[, .(nfrag=.N), keyby=.(scaffold, binsize, pos = start %/% binsize * binsize)]->fragbin
    #     fragbin[, id := paste(sep=":", scaffold, pos)]
    #     fragbin<- hic_map[, .(scaffold, chr, cM=hic_bin)][fragbin, on="scaffold", nomatch=0]
    # fragged = assembly["info"].assign(
    #     binsize=np.maximum(min_binsize, assembly["info"]["length"] // min_nbin)
    # )[["binsize"]].merge(frags, on="scaffold_index")
    # fragbin = fragged.assign(
    #     pos=(fragged["start"] // binsize) * binsize).groupby(["scaffold_index", "binsize", "pos"]).size().toframe("nfrag")
    # fragbin = hic_map[["chr", "hic_bin"]].rename(columns={"hic_bin": "cM"}).merge(fragbin, on="scaffold_index",
    #                                                                               how="inner")
    # #     unique(fragbin[, .(scaffold1=scaffold, binsize1=binsize)])[assembly$fpairs, on='scaffold1']->fp
    # #     unique(fragbin[, .(scaffold2=scaffold, binsize2=binsize)])[fp, on='scaffold2']->fp
    # #     fp[, .(nlinks=.N), keyby=.(scaffold1, pos1 = pos1 %/% binsize1 * binsize1, scaffold2, pos2 = pos2 %/% binsize2 * binsize2)]->binl
    # #     binl[, id1 := paste(sep=":", scaffold1, pos1)]
    # #     binl[, id2 := paste(sep=":", scaffold2, pos2)]
    # #     binl[id1 != id2]->binl

    # left = fragbin[["binsize"]].reset_index(drop=False).rename(columns={"binsize": "binsize1",
    #                                                                     "scaffold_index": "scaffold_index1"}).drop_duplicates()
    # fp = left.merge(assembly["fpairs"], on="scaffold_index1")
    # if fp.index.name == "scaffold_index1":
    #     fp = fp.reset_index(drop=False)
    # left = fragbin[["binsize"]].reset_index(drop=False).rename(columns={"binsize": "binsize2",
    #                                                                     "scaffold_index": "scaffold_index2"}).drop_duplicates()
    # fp = left.merge(fp, on="scaffold_index2")
    # if fp.index.name == "scaffold_index2":
    #     fp = fp.reset_index(drop=False)
    # # Calculate number of links per bin?
    # binl = fp.assign(pos1=(fp["pos1"] / fp["binsize1"]) * fp["binsize1"],
    #                  pos2=(fp["pos2"] / fp["binsize2"]) * fp["binsize2"]).groupby(
    #     ["scaffold_index1", "pos1", "scaffold_index2", "pos2"]).size().to_frame("nlinks")
    # # Ignore self-links
    # binl = binl.query("scaffold_index1 != scaffold_index2 | pos1 != pos2")
    # # Merge with fragbin
    # binl = fragbin.rename(columns={"scaffold_index": "scaffold_index1", "pos": "pos1", "chr": "chr1", "cM": "cM1"})[
    #     ["scaffold_index1", "chr1", "cM1"]].merge(binl, on=["scaffold_index1", "pos1"])
    # binl = fragbin.rename(columns={"scaffold_index": "scaffold_index2", "pos": "pos2", "chr": "chr2", "cM": "cM2"})[
    #     ["scaffold_index2", "chr2", "cM2"]].merge(binl, on=["scaffold_index2", "pos2"])
    # #     fragbin[, .(id1=id, chr1=chr, cM1=cM)][binl, on="id1"]->binl
    # #     fragbin[, .(id2=id, chr2=chr, cM2=cM)][binl, on="id2"]->binl
    # #     binl[, c("scaffold1", "scaffold2", "pos1", "pos2") := list(NULL, NULL, NULL, NULL)]
    # #     setnames(binl, c("id1", "id2"), c("scaffold1", "scaffold2"))

    # # fragbin[, .(scaffold=id, nfrag, chr, cM)]->hic_info_bin
    # #     hic_info_bin[, excluded:=nfrag < min_nfrag_bin]
    # #     binl[chr1 == chr2 & (abs(cM1-cM2) <= 2 | is.na(cM1) | is.na(cM2))]->binl
    # #     binl[, weight:=-log10(nlinks)]
    # #
    # #     make_hic_map(hic_info=hic_info_bin, links=binl, ncores=ncores, maxiter=maxiter, known_ends=known_ends)->hic_map_bin

    # hic_info_bin = fragbin[["scaffold_index", "pos", "nfrag", "chr", "cM"]]
    # hic_info_bin["excluded"] = hic_info_bin["nfrag"] < min_nfrag_bin
    # # TODO Why 2 centimorgans max?
    # binl = binl.query("chr1 == chr2 & (cM1 != cM1 | cM2 != cM2 | (-2 <= cM1 -cM2 <= 2))")[:]
    # # TODO why is the weight negatively correlated with the number of links?
    # binl["weight"] = 1 + np.log10(binl["nlinks"])
    # # TODO why are we recalculating the map?
    # hic_map_bin = make_hic_map(hic_info=hic_info_bin, links=binl, client=client, ncores=ncores, known_ends=known_ends)
