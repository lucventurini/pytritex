import pandas as pd
import numpy as np
import dask.dataframe as dd
import os
import scipy.stats as sps
from .hic_map_constructor import make_hic_map
from dask.distributed import Client


def orient_hic_map(assembly: dict, hic_map: dd.DataFrame, frags: dd.DataFrame, client: Client,
                   min_nbin=5, min_binsize=1e5, binsize=5e5, min_nfrag_bin=30,):

    # assembly$info[, .(scaffold, binsize=pmax(min_binsize, length %/% min_nbin))][frags, on='scaffold']->f
    if isinstance(assembly["fai"], str):
        fai = dd.read_parquet(assembly["fai"], infer_divisions=True)
    else:
        fai = assembly["fai"]
    fai["binsize"] = fai.eval("length / @min_nbin", local_dict={"min_nbin": min_nbin})
    fai["min_nbinsize"] = min_binsize
    fai["binsize"] = fai[["binsize", "min_nbinsize"]].max(axis=1)
    fai = fai.drop("min_nbinsize")
    f = fai.merge(frags, on="scaffold_index")
    f["pos"] = f["start"] // f["binsize"] * f["binsize"]
    #     f[, .(nfrag=.N), keyby=.(scaffold, binsize, pos = start %/% binsize * binsize)]->fragbin
    #     fragbin[, id := paste(sep=":", scaffold, pos)]
    #     fragbin<- hic_map[, .(scaffold, chr, cM=hic_bin)][fragbin, on="scaffold", nomatch=0]
    fragbin = f.groupby(["scaffold_index", "binsize", "pos"]).size().to_frame(
        "nfrag").reset_index(drop=False).set_index("scaffold_index")
    fragbin = hic_map[["chr", "cM"]].rename(columns={"cM": "hic_bin"}).merge(fragbin, on="scaffold_index", how="inner")
    # Now merge with the HiC links
    left = fragbin[["binsize"]].rename(columns={"binsize": "binsize1"})
    left.index = left.index.rename("scaffold_index1")
    left = left.drop_duplicates()

    if isinstance(assembly["fpairs"], str):
        fpairs = dd.read_parquet(assembly["fpairs"], infer_divisions=True)
    else:
        fpairs = assembly["fpairs"]

    #     unique(fragbin[, .(scaffold1=scaffold, binsize1=binsize)])[assembly$fpairs, on='scaffold1']->fp
    #     unique(fragbin[, .(scaffold2=scaffold, binsize2=binsize)])[fp, on='scaffold2']->fp
    fp = left.merge(fpairs.set_index("scaffold_index1"), on="scaffold_index1", how="right").reset_index(drop=False)
    left = left.rename(columns={"binsize1": "binsize2"})
    left.index = left.index.rename("scaffold_index2")
    fp = left.merge(fp.set_index("scaffold_index2"), on="scaffold_index2", how="right").reset_index(drop=False)

    # Calculate binl?
    # TODO: do we actually have pos1 and pos2?
    #     fp[, .(nlinks=.N), keyby=.(scaffold1, pos1 = pos1 %/% binsize1 * binsize1, scaffold2, pos2 = pos2 %/% binsize2 * binsize2)]->binl
    #     binl[, id1 := paste(sep=":", scaffold1, pos1)]
    #     binl[, id2 := paste(sep=":", scaffold2, pos2)]
    #     binl[id1 != id2]->binl
    fp["pos1"] = fp["pos1"] // fp["binsize1"] * fp["binsize1"]
    fp["pos1"] = fp["pos2"] // fp["binsize2"] * fp["binsize2"]
    binl = fp.groupby(["scaffold_index1", "pos1", "scaffold_index2", "pos2"]).size().to_frame("nlinks").reset_index(
        drop=False)

    binl = binl.query("scaffold_index1 != scaffold_index2 | pos1 != pos2")

    # Now merge the fragbin?
    #
    #     fragbin[, .(id1=id, chr1=chr, cM1=cM)][binl, on="id1"]->binl
    #     fragbin[, .(id2=id, chr2=chr, cM2=cM)][binl, on="id2"]->binl
    left = fragbin[["chr", "cM"]].rename(columns={"chr": "chr1", "cM": "cM1"})
    left.index = left.index.rename("scaffold_index1")
    binl = left.merge(binl.set_index("scaffold_index1"), on="scaffold_index1").reset_index(drop=False)

    left = fragbin[["chr", "cM"]].rename(columns={"chr": "chr2", "cM": "cM2"})
    left.index = left.index.rename("scaffold_index2")
    binl = left.merge(binl.set_index("scaffold_index2"), on="scaffold_index2").reset_index(drop=False)
    # TODO: Completely IGNORE the madness of creating a new "scaffold" index based on scaffold + ":" + pos
    #     fragbin[, .(scaffold=id, nfrag, chr, cM)]->hic_info_bin
    #     hic_info_bin[, excluded:=nfrag < min_nfrag_bin]
    hic_info_bin = fragbin[["nfrag", "chr", "cM"]]
    fragbin = fragbin.eval("excluded = nfrag < @min_nfrag_bin", local_dict={"min_nfrag_bin": min_nfrag_bin})
    # binl[chr1 == chr2 & (abs(cM1-cM2) <= 2 | is.na(cM1) | is.na(cM2))]->binl
    binl["cMdist"] = (binl["cM1"] - binl["cM2"]).abs()
    binl = binl.query("chr1 == chr2 | cMDist != cMDist | cMDist <= 2")
    # binl[, weight:=-log10(nlinks)]
    binl["weight"] = -binl["nlinks"].apply(np.log10, meta=float)
    # make_hic_map(hic_info=hic_info_bin, links=binl,
    #  ncores=ncores, maxiter=maxiter, known_ends=known_ends)->hic_map_bin
    hic_map_bin = make_hic_map(hic_info=hic_info_bin, links=binl, client=client,
                               save_dir=None, ncores=1, maxiter=100, known_ends=True)
    # Yet again another table merging in preparation here ...

    #  w<-hic_map_bin[!is.na(hic_bin),
    #                .(id=scaffold, scaffold=sub(":.*$", "", scaffold), pos=as.integer(sub("^.*:", "", scaffold)), chr, hic_bin)]
    # w<-w[, .(gbin=mean(na.omit(hic_bin)),
    # 	     hic_cor=as.numeric(suppressWarnings(cor(method='s', hic_bin, pos, use='p')))), keyby=scaffold][!is.na(hic_cor)]

    w = hic_map_bin.query("hic_bin == hic_bin")[["pos", "chr", "hic_bin"]]
    grouped = w.groupby("scaffold_index")
    cors = grouped[["hic_bin", "pos"]].apply(lambda group: sps.spearmanr(group["hic_bin"].values, group["pos"].values)[0],
                                             meta=float).to_frame("hic_cor")
    grouped = grouped["hic_bin"].mean().to_frame("gbin").merge(cors)
    w = w[["chr", "hic_bin"]].merge(grouped, left_on="scaffold_index", right_on="scaffold_index")
    if w.index.name != "scaffold_index":
        assert w.index.name is None
        w = w.set_index("scaffold_index")

    # Another dimension now
    # hic_map[!is.na(hic_bin) & scaffold %in% w$scaffold][order(chr, hic_bin)]->z0
    #     z0[,.(scaffold1=scaffold[1:(.N-2)], scaffold2=scaffold[2:(.N-1)], scaffold3=scaffold[3:(.N)]), by=chr]->z
    #     z0[, data.table(key="scaffold1", scaffold1=scaffold, hic_bin1=hic_bin)][setkey(z, "scaffold1")]->z
    #     z0[, data.table(key="scaffold2", scaffold2=scaffold, hic_bin2=hic_bin)][setkey(z, "scaffold2")]->z
    #     z0[, data.table(key="scaffold3", scaffold3=scaffold, hic_bin3=hic_bin)][setkey(z, "scaffold3")]->z
    z0 = hic_map.query("hic_bin == hic_bin").loc[w.index].compute().sort_values(["chr", "hic_bin"])

    #     z0[,.(scaffold1=scaffold[1:(.N-2)],
    #     scaffold2=scaffold[2:(.N-1)], scaffold3=scaffold[3:(.N)]), by=chr]->z
    z = z0.reset_index(drop=False)
    zgrouped = z.groupby("chr")
    z["scaffold_index1"] = zgrouped["scaffold_index"].shift(0)
    z["scaffold_index2"] = zgrouped["scaffold_index"].shift(-1)
    z["scaffold_index3"] = zgrouped["scaffold_index"].shift(-2)
    z = z.dropna(subset=["scaffold_index1", "scaffold_index2"])
    # Now merge with w
    #     w[, data.table(key="scaffold1", scaffold1=scaffold, gbin1=gbin)][setkey(z, "scaffold1")]->z
    left = w[["gbin"]].rename(columns={"gbin": "gbin1"})
    left.index = left.index.rename("scaffold_index1")
    z = left.merge(z.set_index("scaffold_index1"), on="scaffold_index1").reset_index(drop=False)
    #     w[, data.table(key="scaffold2", scaffold2=scaffold, gbin2=gbin)][setkey(z, "scaffold2")]->z
    left = w[["gbin"]].rename(columns={"gbin": "gbin2"})
    left.index = left.index.rename("scaffold_index2")
    z = left.merge(z.set_index("scaffold_index2"), on="scaffold_index2").reset_index(drop=False)
    #     w[, data.table(key="scaffold3", scaffold3=scaffold, gbin3=gbin)][setkey(z, "scaffold3")]->z
    left = w[["gbin"]].rename(columns={"gbin": "gbin3"})
    left.index = left.index.rename("scaffold_index3")
    z = left.merge(z.set_index("scaffold_index3"), on="scaffold_index3").reset_index(drop=False)
    #     z[, cc:= apply(z[, .(hic_bin1, hic_bin2, hic_bin3, gbin1, gbin2, gbin3)],1,function(x) {
    # 		    suppressWarnings(cor(x[1:3], x[4:6]))
    # 		    })]
    

    #     z[, data.table(key="scaffold", scaffold=scaffold2, cc=ifelse(cc > 0, 1, -1))]->ccor
    #     ccor[w]->m
    #     m[, hic_orientation:=ifelse(hic_cor > 0, 1 * cc, -1 * cc)]
    #     m[, .(scaffold, hic_cor, hic_invert=cc, hic_orientation)][hic_map, on="scaffold"]->hic_map_oriented
    #     hic_map_oriented[is.na(hic_orientation), hic_orientation := ifelse(hic_cor > 0, 1, -1)]
    #     setnames(hic_map_oriented, "chr", "consensus_chr")
    #     setnames(hic_map_oriented, "cM", "consensus_cM")
    #     hic_map_oriented[, consensus_orientation := hic_orientation]

    pass