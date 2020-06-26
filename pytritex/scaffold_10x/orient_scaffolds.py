import pandas as pd
import numpy as np
from .make_agp import make_agp
import dask.dataframe as dd
import dask.array as da
from dask.distributed import Client
from dask.delayed import delayed
import os

#  m[super_nbin > 1, .(scaffold1=scaffold, bin1=bin, super1=super)][link_pos, on="scaffold1", nomatch=0]->a
#  m[super_nbin > 1, .(scaffold2=scaffold, bin2=bin, super2=super)][a, on="scaffold2", nomatch=0]->a
#  a[super1 == super2 & bin1 != bin2]->a
#  a[, d := abs(bin2 - bin1)]->a
#  a[d <= max_dist_orientation]->a
#  a[, .(nxt=mean(pos1[bin2 > bin1]), prv=mean(pos1[bin2 < bin1])), key=.(scaffold=scaffold1)]->aa
#  info[, .(scaffold, length)][aa, on="scaffold"]->aa
#  aa[!is.nan(prv) & !is.nan(nxt), orientation := ifelse(prv <= nxt, 1, -1)]
#  aa[is.nan(prv), orientation := ifelse(length - nxt <= nxt, 1, -1)]
#  aa[is.nan(nxt), orientation := ifelse(length - prv <= prv, -1, 1)]
#
#  aa[, .(orientation, scaffold)][m, on="scaffold"]->m
#  m[, oriented := T]
#  m[is.na(orientation), oriented := F]
#  m[is.na(orientation), orientation := 1]
#  setorder(m, super, bin, rank)
#
#  m[, super_pos := 1 + cumsum(c(0, length[-.N])), by=super]
#
#  if(verbose){
#   cat("Anchoring super-scaffolds.\n")
#  }
#  # assign super scaffolds to genetic positions
#  m[!is.na(chr), .(nchr=.N), key=.(chr, super)][, pchr := nchr/sum(nchr), by=super][order(-nchr)][!duplicated(super)]->y
#  m[!is.na(cM)][y, on=c("chr", "super")]->yy
#  y[res, on="super"]->res
#  yy[, .(cM=mean(cM), min_cM=min(cM), max_cM=max(cM)), key=super][res, on='super']->res
#  setorder(res, -length)


def super_position(group):
    orig_index = group.index
    group = group.sort_values(["bin", "rank"])
    # return group[["length"]].values
    vals = (group["length"].shift(1, fill_value=0).cumsum() + 1).astype(np.int)
    vals = vals[orig_index]
    return vals.values


def chrom_percentage(group):
    chr_sum = group["nchr"].sum()
    return (group["nchr"] / chr_sum).values


def orient_scaffolds(info: str, res: str,
                     membership: str,
                     link_pos: str,
                     client: Client,
                     max_dist_orientation: float,
                     save_dir: str):

    membership = dd.read_parquet(membership, infer_divisions=True)
    link_pos = dd.read_parquet(link_pos, infer_divisions=True)
    info = dd.read_parquet(info, infer_divisions=True)
    res = dd.read_parquet(res, infer_divisions=True)

    # #  m[super_nbin > 1, .(scaffold1=scaffold, bin1=bin, super1=super)][link_pos, on="scaffold1", nomatch=0]->a
    m_greater_one = membership.query("super_nbin > 1")
    assert m_greater_one.shape[0].compute() > 0, membership.head()
    left = m_greater_one[["bin", "super"]].rename(columns={"bin": "bin1", "super": "super1"})
    left.index = left.index.rename("scaffold_index1")
    left = client.scatter(left)
    link_pos = client.scatter(link_pos)
    func = delayed(dd.merge)(left, link_pos, on="scaffold_index1", how="inner")
    association = client.scatter(client.compute(func).result())

    left = m_greater_one[["bin", "super"]].rename(columns={"bin": "bin2", "super": "super2"})
    left.index = left.index.rename("scaffold_index2")
    left = client.scatter(left)
    func = delayed(dd.merge)(left, association, on="scaffold_index2", how="inner")
    association = client.compute(func).result()

    association = association.query("(super1 == super2) & (bin1 != bin2)")
    association = association.eval("d = abs(bin2 - bin1)")
    association = association[association["d"] <= max_dist_orientation]
    association = association.rename(columns={"scaffold_index1": "scaffold_index"})
    assert association.shape[0].compute() > 0

    grouped = association.groupby("scaffold_index")
    nxt = grouped.apply(
        lambda group: group.loc[group.eval("bin2 > bin1"), "pos1"].mean(), meta=float).to_frame("nxt")
    prev = grouped.apply(
        lambda group: group.loc[group.eval("bin2 < bin1"), "pos1"].mean(), meta=float).to_frame("prv")
    association = dd.concat([prev, nxt], axis=1)
    # assert aa.shape[0].compute() > 0
    association = info[["length"]].merge(
        association, left_index=True, right_index=True).reset_index(drop=False)
    assert association.shape[0].compute() > 0
    # "x != x" means: return the np.nan indices
    idx1 = association.eval("(prv == prv) & (nxt == nxt)")
    index = association.index.compute()
    orientation = pd.Series([np.nan] * index.shape[0], index=index)
    orientation.loc[idx1] = association.loc[idx1].eval("prv <= nxt").compute()
    idx2 = association.eval("(prv != prv) & (nxt == nxt)")
    orientation.loc[idx2] = association.loc[idx2].eval("length - nxt <= nxt").compute()
    idx3 = association.eval("(prv == prv) & (nxt != nxt)")
    orientation.loc[idx3] = association.loc[idx3].eval("prv < length - prv").compute()
    orientation = orientation.map({True: 1, False: -1})
    chunks = tuple(association.map_partitions(len).compute().values.tolist())
    orientation = da.from_array(orientation, chunks=chunks)
    association["orientation"] = orientation
    association = association.set_index("scaffold_index")
    # assert "scaffold_index" in membership.columns
    assert "scaffold_index" == association.index.name == membership.index.name, (association.index.name,
                                                                        membership.index.name)
    membership = dd.merge(association[["orientation"]],
                          membership, left_index=True, right_index=True)

    # Now assign an "orientation" value to each (1 if + or unknown, -1 for -)
    # and an "oriented" value (strictly boolean)
    chunks = tuple(membership.map_partitions(len).compute().values.tolist())
    orientation = membership["orientation"].compute()
    oriented = pd.Series([True] * orientation.shape[0], index=orientation,
                         dtype=bool)
    idx1 = orientation.isna()
    oriented.loc[idx1] = False
    orientation.loc[idx1] = 1
    orientation = da.from_array(orientation, chunks=chunks)
    oriented = da.from_array(oriented, chunks=chunks)
    membership["orientation"] = orientation
    membership["oriented"] = oriented
    # Next step: assign to each scaffold a base-pair position in the super-scaffold, "super_pos"
    # This position should be determined by bin and rank.

    super_pos = membership.groupby("super").apply(super_position, meta=np.int).compute()
    super_pos = np.concatenate(super_pos.values)
    chunks = tuple(membership.map_partitions(len).compute().values.tolist())
    super_pos = da.from_array(super_pos, chunks=chunks)
    membership["super_pos"] = super_pos

    # Now we have to assign the genetic positions, ie anchor
    anchored = membership[~membership["chr"].isna()]
    nchr = anchored.groupby(["chr", "super"])["length"].size().to_frame("nchr").reset_index(drop=False)
    nchr_sum = nchr.groupby("super")["nchr"].transform("sum", meta=np.int)
    nchr["pchr"] = nchr["nchr"] / nchr_sum.compute()
    nchr["max_nchr"] = nchr.groupby("super")["nchr"].transform("max", meta=np.int).values
    nchr = nchr.query("nchr == max_nchr").drop_duplicates(subset="super").drop("max_nchr", axis=1)
    nchr_with_cms = dd.merge(membership[~membership["cM"].isna()],
                             nchr, on=["chr", "super"])
    grouped_nchr_with_cms = nchr_with_cms.groupby("super")
    nchr_with_cms["min_cM"] = grouped_nchr_with_cms["cM"].transform("min")
    nchr_with_cms["max_cM"] = grouped_nchr_with_cms["cM"].transform("max")
    nchr_with_cms["cM"] = grouped_nchr_with_cms["cM"].transform("mean")
    res = dd.merge(nchr, res, on=["super"])
    assert "super" in res.columns()
    res = res.persist()
    res = dd.merge(nchr_with_cms, res, on="super")
    res_name = os.path.join(save_dir, "res")
    dd.to_parquet(res, res_name, compression="gzip", compute=True, engine="pyarrow")
    mem_name = os.path.join(save_dir, "membership")
    dd.to_parquet(membership, mem_name, compression="gzip", compute=True, engine="pyarrow")
    return mem_name, res_name
