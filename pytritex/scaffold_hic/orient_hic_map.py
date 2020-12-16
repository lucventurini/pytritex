import pandas as pd
import numpy as np
import dask.dataframe as dd
import os
from .hic_map_constructor import make_hic_map
from dask.distributed import Client
from scipy.stats import spearmanr, pearsonr


def orient_hic_map(info, assembly: dict, hic_map: dd.DataFrame, frags: dd.DataFrame, client: Client,
                   min_nfrag_bin=30, cores=1, save_dir=None,
                   maxiter=100, orient_old=False, min_nbin=5, min_binsize=1e5):

    # assembly$info[, .(scaffold, binsize=pmax(min_binsize, length %/% min_nbin))][frags, on='scaffold']->f
    if isinstance(assembly["fai"], str):
        fai = dd.read_parquet(assembly["fai"], infer_divisions=True)
    else:
        fai = assembly["fai"]
    fai["binsize"] = fai.eval("length / @min_nbin", local_dict={"min_nbin": min_nbin})
    fai["min_nbinsize"] = min_binsize
    fai["binsize"] = fai[["binsize", "min_nbinsize"]].max(axis=1).drop("min_nbinsize", axis=1)
    f = fai.merge(frags, on="scaffold_index")
    f["pos"] = f["frag_start"] // f["binsize"] * f["binsize"]
    if orient_old is False:
        hic_map_oriented = _new_orientation(info, assembly, hic_map, frags, client, min_nfrag_bin=min_nfrag_bin,
                                            cores=cores, maxiter=maxiter, min_nbin=min_nbin, min_binsize=min_binsize)
    else:
        hic_map_oriented = _use_old_orient(info, assembly, hic_map, frags, client, min_nfrag_bin=min_nfrag_bin,
                                            cores=cores, maxiter=maxiter, min_nbin=min_nbin, min_binsize=min_binsize)

    # hic_map_oriented[assembly$info, on="scaffold"]->hic_map_oriented
    hic_map_oriented = hic_map_oriented.merge(assembly["info"], on="scaffold_index")
    if save_dir is not None:
        dd.to_parquet(hic_map_oriented, os.path.join(save_dir, "hic_map_oriented"))
        return os.path.join(save_dir, "hic_map_oriented")

    return hic_map_oriented


def _use_old_orient(*args, **kwargs):
    raise NotImplementedError()


def _new_orientation(info, assembly: dict, hic_map: dd.DataFrame, frags: dd.DataFrame, client: Client,
                     cores=1, binsize=5e5, min_nfrag_bin=30,
                     known_ends=True, maxiter=100, min_nbin=5, min_binsize=1e5):

    # info = assembly["info"]
    if isinstance(info, str):
        info = dd.read_parquet(info, infer_divisions=True)
    else:
        assert isinstance(info, (pd.DataFrame, dd.DataFrame))

    if info.index is None:
        info = info.set_index("scaffold_index")
    else:
        assert info.index.name == "scaffold_index"

    #     assembly$info[, .(scaffold, binsize=pmax(min_binsize, length %/% min_nbin))][frags, on='scaffold']->f
    temp_info = np.maximum(info[["length"]] // min_nbin, min_binsize).rename(columns={"length": "binsize"})
    fragbin = dd.merge(temp_info, frags, how="right", on="scaffold_index")
    #     f[, .(nfrag=.N), keyby=.(scaffold, binsize, pos = start %/% binsize * binsize)]->fragbin
    #     fragbin[, id := paste(sep=":", scaffold, pos)]
    #     fragbin<- hic_map[, .(scaffold, chr, cM=hic_bin)][fragbin, on="scaffold", nomatch=0]
    fragbin["pos"] = fragbin["frag_start"] // binsize * binsize
    fragbin = fragbin.groupby(["scaffold_index", "binsize", "pos"]).size().to_frame("nfrag").reset_index(drop=False)
    fragbin = dd.merge(hic_map[["chr", "hic_bin"]].rename(colums={"hic_bin": "cM"}), fragbin, on="scaffold_index",
                       how="inner")
    #     unique(fragbin[, .(scaffold1=scaffold, binsize1=binsize)])[assembly$fpairs, on='scaffold1']->fp
    left = fragbin[["binsize"]].reset_index(drop=False).unique().set_index("scaffold_index")
    left = left.rename(columns={"binsize": "binsize1"})
    left.index = left.index.rename("scaffold_index1")
    fragged_fpairs = dd.merge(left, assembly["fpairs"], on="scaffold_index1")
    #     unique(fragbin[, .(scaffold2=scaffold, binsize2=binsize)])[fp, on='scaffold2']->fp
    left = left.rename(columns={"binsize": "binsize2"})
    left.index = left.index.rename("scaffold_index2")
    fragged_fpairs = dd.merge(left, fragged_fpairs, on="scaffold_index2")
    # fp[, .(nlinks=.N), keyby=.(scaffold1, pos1 = pos1 %/% binsize1 * binsize1, scaffold2, pos2 = pos2 %/% binsize2 * binsize2)]->binl
    fragged_fpairs = fragged_fpairs.assign(
        pos1=fragged_fpairs["pos1"] // fragged_fpairs["binsize1"] * fragged_fpairs["binsize1"],
        pos2=fragged_fpairs["pos2"] // fragged_fpairs["binsize2"] * fragged_fpairs["binsize2"])
    #     binl[, id1 := paste(sep=":", scaffold1, pos1)]
    #     binl[, id2 := paste(sep=":", scaffold2, pos2)]
    #     binl[id1 != id2]->binl
    binned_fpairs = fragged_fpairs.groupby(
        ["scaffold_index1", "scaffold_index2", "pos1", "pos2"]).size().to_frame("nlinks")
    binned_fpairs = binned_fpairs.query("chr1 != chr2 | pos1 != pos2")  # remove self-links
    #     fragbin[, .(id1=id, chr1=chr, cM1=cM)][binl, on="id1"]->binl
    #     fragbin[, .(id2=id, chr2=chr, cM2=cM)][binl, on="id2"]->binl
    left = fragbin[["chr", "pos", "cM"]].rename(columns={"chr": "chr1", "pos": "pos1", "cM": "cM1"})
    left.index = left.index.rename("scaffold_index1")
    binned_fpairs = dd.merge(left, binned_fpairs, on=["scaffold_index1", "pos1"], how="right")
    left = fragbin[["chr", "pos", "cM"]].rename(columns={"chr": "chr2", "pos": "pos2", "cM": "cM2"})
    left.index = left.index.rename("scaffold_index2")
    binned_fpairs = dd.merge(left, binned_fpairs, on=["scaffold_index2", "pos2"], how="right")
    #     fragbin[, .(scaffold=id, nfrag, chr, cM)]->hic_info_bin
    #     hic_info_bin[, excluded:=nfrag < min_nfrag_bin]
    #     binl[chr1 == chr2 & (abs(cM1-cM2) <= 2 | is.na(cM1) | is.na(cM2))]->binl
    #     binl[, weight:=-log10(nlinks)]
    hic_info_bin = fragbin[["pos", "nfrag", "chr", "cM"]].eval("excluded = nfrag < @min_nfrag_bin",
                                                               local_dict={"min_nfrag_bin": min_nfrag_bin})
    binned_fpairs = binned_fpairs.query(
        "chr1 == chr2 & (cM1 != cM1 | cM2 != cM2 | (cM1 - cM2 >= -2 & cM1 - cM2 <= 2)").eval(
        "weight = - log(nlinks) / log(10)")
    #  make_hic_map(hic_info=hic_info_bin, links=binl, ncores=ncores, maxiter=maxiter, known_ends=known_ends)->hic_map_bin
    hic_map_bin = make_hic_map(hic_info=hic_info_bin, links=binned_fpairs, ncores=cores, maxiter=maxiter,
                               save_dir=None,
                               known_ends=known_ends, client=client)

    # TODO what is this w?
    #     w<-hic_map_bin[!is.na(hic_bin), .(id=scaffold, scaffold=sub(":.*$", "", scaffold), pos=as.integer(sub("^.*:", "", scaffold)), chr, hic_bin)]
    #     w<-w[, .(gbin=mean(na.omit(hic_bin)),
    # 	     hic_cor=as.numeric(suppressWarnings(cor(method='s', hic_bin, pos, use='p')))), keyby=scaffold][!is.na(hic_cor)]
    w = hic_map_bin.query("hic_bin == hic_bin")[["pos", "chr", "hic_bin"]]
    grouped = w.groupby("scaffold_index")
    gbin = grouped["hic_bin"].mean().to_frame("gbin")
    hic_cor = grouped[["hic_bin", "pos"]].apply(lambda group: spearmanr(group["hic_bin"].values, group["pos"].values)[0],
                                                meta=float)
    w = dd.merge(gbin, hic_cor, on="scaffold_index").query("hic_cor == hic_cor")
    w_scaffolds = np.unique(w.index.values)
    #     hic_map[!is.na(hic_bin) & scaffold %in% w$scaffold][order(chr, hic_bin)]->z0
    z0 = hic_map.loc[w_scaffolds].query("hic_bin == hic_bin").compute().sort_values(["chr", "hic_bin"])
    z0 = z0.reset_index(drop=False)
    # z0[,.(scaffold1=scaffold[1:(.N-2)], scaffold2=scaffold[2:(.N-1)],
    # scaffold3=scaffold[3:(.N)]), by=chr]->z
    scaf_dfs = dict()
    for num in (1, 2, 3):
        shifted = 3 - num
        scaf_dfs[num] = z0.groupby("chr")[["scaffold_index"]].shift(shifted).rename(
            columns={"scaffold_index": "scaffold_index{num}".format(num=num)})

    # Bring back the information regarding the hic_bin
    # z0[, data.table(key="scaffold1", scaffold1=scaffold, hic_bin1=hic_bin)][setkey(z, "scaffold1")]->z
    # z0[, data.table(key="scaffold2", scaffold2=scaffold, hic_bin2=hic_bin)][setkey(z, "scaffold2")]->z
    # z0[, data.table(key="scaffold3", scaffold3=scaffold, hic_bin3=hic_bin)][setkey(z, "scaffold3")]->z
    z = pd.concat(scaf_dfs.values(), axis=1).dropna().astype(int)
    for num in (1, 2, 3):
        left = z0[["hic_bin"]].rename(columns={"hic_bin": "hic_bin{num}".format(num=num)})
        left.index = left.index.rename("scaffold_index{num}".format(num=num))
        z = left.merge(z.set_index("scaffold_index{num}".format(num=num)),
                       on="scaffold_index{num}".format(num=num)).reset_index(drop=False)

    # Bring back the genomic bin
    #  w[, data.table(key="scaffold1", scaffold1=scaffold, gbin1=gbin)][setkey(z, "scaffold1")]->z
    #  w[, data.table(key="scaffold2", scaffold2=scaffold, gbin2=gbin)][setkey(z, "scaffold2")]->z
    #  w[, data.table(key="scaffold3", scaffold3=scaffold, gbin3=gbin)][setkey(z, "scaffold3")]->z
    for num in (1, 2, 3):
        left = w[["gbin"]].rename(columns={"gbin": "gbin{num}".format(num=num)})
        left.index = left.index.rename("scaffold_index{num}".format(num=num))
        z = left.merge(z.set_index("scaffold_index{num}".format(num=num)),
                       on="scaffold_index{num}".format(num=num)).reset_index(drop=False)

    # Now calculate the correlation between genomic bins and hic_bins

    #     z[, cc:= apply(z[, .(hic_bin1, hic_bin2, hic_bin3, gbin1, gbin2, gbin3)],1,function(x) {
    # 		    suppressWarnings(cor(x[1:3], x[4:6]))
    # 		    })]
    #     z[, data.table(key="scaffold", scaffold=scaffold2, cc=ifelse(cc > 0, 1, -1))]->ccor
    #     ccor[w]->m
    #     m[, hic_orientation:=ifelse(hic_cor > 0, 1 * cc, -1 * cc)]
    left_cols = z[[f"hic_bin{num}" for num in (1, 2, 3)]].values
    right_cols = z[[f"gbin{num}" for num in (1, 2, 3)]].values
    cols = np.hstack([left_cols.compute(), right_cols.compute()])
    z["cc"] = dd.from_array(
        np.array([pearsonr(cols[index, :3], cols[index, 3:])[0] for index in np.arange(cols.shape[0], dtype=int)]))
    # Grab the central scaffold
    ccor = z[["scaffold_index2", "cc"]].rename(columns={"scaffold_index2": "scaffold_index"}
                                               ).set_index("scaffold_index")
    ccor["cc"] = ccor["cc"].mask(ccor["cc"] >= 0, 1)
    ccor["cc"] = ccor["cc"].mask(ccor["cc"] < 0, -1)
    m = ccor.merge(w, on="scaffold_index")
    # Reverse the orientation if the hic_cor was negative
    m["hic_orientation"] = m["cc"].mask(m["hic_cor"] <= 0, -1 * m["cc"])

    # m[, .(scaffold, hic_cor, hic_invert=cc, hic_orientation)][hic_map, on="scaffold"]->hic_map_oriented
    hic_map_oriented = m[["hic_cor", "cc", "hic_orientation"]].rename(columns={"cc": "hic_invert"}).merge(
        hic_map, on="scaffold_index")
    #     hic_map_oriented[is.na(hic_orientation), hic_orientation := ifelse(hic_cor > 0, 1, -1)]
    hic_map_oriented["hic_orientation"] = hic_map_oriented["hic_orientation"].mask(
        hic_map_oriented["hic_orientation"].isna() & hic_map_oriented["hic_cor"] > 0, 1)
    hic_map_oriented["hic_orientation"] = hic_map_oriented["hic_orientation"].mask(
        hic_map_oriented["hic_orientation"].isna() & hic_map_oriented["hic_cor"] <= 0, -1)
    #     setnames(hic_map_oriented, "chr", "consensus_chr")
    #     setnames(hic_map_oriented, "cM", "consensus_cM")
    #     hic_map_oriented[, consensus_orientation := hic_orientation]
    hic_map_oriented = hic_map_oriented.rename(columns={"chr": "consensus_chr",
                                                        "cM": "consensus_cM"}).eval(
        "consensus_orientation = hic_orientation")
    return hic_map_oriented
