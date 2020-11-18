import pandas as pd
import dask.dataframe as dd
from dask.distributed import Client
import numpy as np
from ..utils.rolling_join import rolling_join


def add_psmol_fpairs(assembly: dict, hic_map: dict, nucfile: str,
                     map_10x=None, assembly_10x=None, cov=None):
    #  if(is.null(map_10x)){
    #   assembly$fpairs[, .(scaffold1, scaffold2, pos1, pos2)] -> z
    #  } else {
    #   assembly_10x$fpairs[, .(scaffold1, scaffold2, pos1, pos2)] -> z
    #  }
    #  hic_map$agp[agp_chr != "chrUn", .(chr, scaffold, orientation=orientation, agp_start, agp_end)]->a
    #  a[is.na(orientation), orientation := 1]
    #  setnames(copy(a), paste0(names(a), 1))[z, on="scaffold1"]->z
    #  setnames(copy(a), paste0(names(a), 2))[z, on="scaffold2"]->z
    #  z[orientation1 == 1, pos1 := agp_start1 - 1 + pos1]
    #  z[orientation1 == -1, pos1 := agp_end1 + 1 - pos1]
    #  z[orientation2 == 1, pos2 := agp_start2 - 1 + pos2]
    #  z[orientation2 == -1, pos2 := agp_end2 + 1 - pos2]
    #  z[!is.na(chr1) & !is.na(chr2), .(chr1, chr2, start1=pos1, start2=pos2)]->links

    fpairs_keys = ["scaffold_index1", "scaffold_index2", "pos1", "pos2"]
    if map_10x is None:
        fpairs = assembly["fpairs"][fpairs_keys]
    else:
        fpairs = assembly_10x["fpairs"][fpairs_keys]

    # TODO double check all these keys
    agp = hic_map["agp"].query("agp_chr != 'chrUn' & agp_chr == agp_chr")[
        ["chr", "scaffold_index", "orientation", "agp_start", "agp_end"]]

    assert isinstance(agp, (pd.DataFrame, dd.DataFrame))
    agp["orientation"] = agp["orientation"].fillna(1)
    a1 = agp.rename(columns=[column + "1" for column in agp.columns])
    a1.index = a1.index.rename(a1.index + "1")
    fpairs = a1.merge(fpairs, on="scaffold_index1")
    if fpairs.index == "scaffold_index1":
        fpairs = fpairs.reset_index(drop=False)
    a2 = agp.rename(columns=[column + "2" for column in agp.columns])
    a2.index = a1.index.rename(a1.index + "2")
    fpairs = a2.merge(fpairs, on="scaffold_index2")
    if fpairs.index == "scaffold_index2":
        fpairs = fpairs.reset_index(drop=False)
    bait = fpairs["orientation"] == 1
    fpairs["pos1"] = fpairs["pos1"].mask(bait, fpairs["agp_start1"] - 1 + fpairs["pos1"])
    fpairs["pos2"] = fpairs["pos2"].mask(bait, fpairs["agp_start2"] - 1 + fpairs["pos2"])
    bait = fpairs["orientation"] == -1
    fpairs["pos1"] = fpairs["pos1"].mask(bait, fpairs["agp_start1"] + 1 - fpairs["pos1"])
    fpairs["pos2"] = fpairs["pos2"].mask(bait, fpairs["agp_start2"] + 1 - fpairs["pos2"])
    links = fpairs.query("chr1 == chr1 & chr2 == chr2")[["chr1", "chr2", "pos1", "pos2"]].rename(
        columns={"start1": "pos1", "start2": "pos2"})
    #  n<-c("orig_scaffold", "orig_start", "orig_end", "frag_id", "nA", "nC", "nG", "nT", "nN", "length")
    #  nuc<-fread(nucfile, select=c(1:4,7:11,13), head=T, col.names=n)
    #  assembly$info[, .(scaffold, orig_scaffold, orig_start, off=orig_start)][
    #  nuc, on=c("orig_scaffold", "orig_start"), roll=T]->z
    #  z[, start := orig_start - off + 1]
    #  z[, end := orig_end - off + 1]

    colnames = ["orig_scaffold", "orig_start", "orig_end", "frag_id", "nA", "nC", "nG", "nT", "nN", "length"]
    nuc = dd.read_csv(nucfile, usecols=list(range(3)) + list(range(6, 11)) + [12],
                      header=0, names=colnames)
    left = assembly["info"].reset_index(drop=False).set_index("orig_scaffold")[
        ["scaffold_index", "orig_scaffold_index", "orig_start"]]
    left["off"] = left["orig_start"]
    frags = rolling_join(left, nuc, on="orig_scaffold", by="orig_start")
    frags["start"] = frags["orig_start"] - frags["off"] + 1
    frags["end"] = frags["orig_end"] - frags["off"] + 1

    if map_10x is not None:
        #  if(!is.null(map_10x)){
        #   map_10x$agp[gap == F, .(super, orientation, super_start, super_end, scaffold)][z, on="scaffold"]->z
        #   z[orientation == 1, start := super_start - 1 + start]
        #   z[orientation == 1, end := super_start - 1 + end]
        #   z[orientation == -1, start := super_end - end + 1]
        #   z[orientation == -1, end := super_end - start + 1]
        #   z[, c("orientation", "super_start", "super_end", "scaffold") := list(NULL, NULL, NULL, NULL)]
        left = map_10x["agp"].query("gap == False")[
            ["super", "orientation", "super_start", "super_end", "scaffold_index"]]
        frags = left.merge(frags, on="scaffold_index")
        bait = (frags["orientation"] == 1)
        frags["start"] = frags["start"].mask(bait, frags["super_start"] - 1 + frags["start"])
        frags["end"] = frags["start"].mask(bait, frags["super_start"] - 1 + frags["end"])
        bait = (frags["orientation"] == -1)
        frags["start"] = frags["start"].mask(bait, frags["super_end"] + 1 - frags["end"])
        frags["end"] = frags["start"].mask(bait, frags["super_end"] + 1 - frags["start"])
        frags = frags.drop(["orientation", "super_start", "super_end", "scaffold_index"], axis=1)
        #   setnames(z, "super", "orig_scaffold")
        frags = frags.rename(columns={"super": "orig_scaffold_index"})  # TODO Not sure at all about this
        # #   assembly_10x$info[, .(scaffold, start=orig_start, orig_start, orig_scaffold)][
        # z, on=c("orig_scaffold", "start"), roll=T]->z
        left = assembly_10x["info"][["orig_scaffold_index", "orig_start"]]
        left["start"] = left["orig_start"]
        left = left.reset_index(drop=False).set_index("orig_scaffold_index")
        #   z[, start := start - orig_start + 1]
        #   z[, end := end - orig_start + 1]
        #   z[, orig_start := NULL]
        #   z[, orig_scaffold := NULL]
        frags = rolling_join(left, frags, on="orig_scaffold_index", by="start")
        frags["start"] = frags["start"] - frags["orig_start"] + 1
        frags["end"] = frags["end"] - frags["orig_start"] + 1
        frags = frags.drop(["orig_start"], axis=1)
        frags = frags.set_index("scaffold_index")

    #  a[z, on="scaffold", nomatch=0]->z
    frags = agp.merge(frags, on="scaffold_index", how="right")
    #  z[orientation == 1, start := agp_start - 1 + start]
    #  z[orientation == -1, start := agp_end + 1 - start]
    #  z[orientation == 1, end := agp_start - 1 + end]
    #  z[orientation == -1, end := agp_end + 1 - end]
    bait = frags["orientation"] == 1
    frags["start"] = frags["start"].mask(bait, frags["agp_start"] - 1 + frags["start"])
    frags["end"] = frags["end"].mask(bait, frags["agp_start"] - 1 + frags["end"])
    bait = frags["orientation"] == -1
    frags["start"] = frags["start"].mask(bait, frags["agp_end"] + 1 - frags["start"])
    frags["end"] = frags["end"].mask(bait, frags["agp_end"] - 1 - frags["end"])
    frags = frags[["chr", "start", "end", colnames[3:]]]
    frags["cov"] = 1

    res = hic_map.copy()
    res["links"] = links
    res["frags"] = frags
    return res
