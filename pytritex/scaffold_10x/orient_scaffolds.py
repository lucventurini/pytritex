import pandas as pd
import numpy as np
from .make_agp import make_agp


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


def orient_scaffolds(info: pd.DataFrame, res: pd.DataFrame,
                     membership: pd.DataFrame,
                     link_pos: pd.DataFrame,
                     max_dist_orientation: float, verbose=False):

    if verbose:
        print("Orienting scaffold indices")

    # #  m[super_nbin > 1, .(scaffold1=scaffold, bin1=bin, super1=super)][link_pos, on="scaffold1", nomatch=0]->a
    m_greater_one = membership.query("super_nbin > 1")
    assert m_greater_one.shape[0] > 0, membership.head()
    left = m_greater_one.loc[:, ["scaffold_index", "bin", "super"]].add_suffix("1")
    a = left.merge(link_pos, on="scaffold_index1", how="inner")
    left = m_greater_one.loc[:, ["scaffold_index", "bin", "super"]].add_suffix("2")
    a = left.merge(a, on="scaffold_index2", how="inner")
    a.query("(super1 == super2) & (bin1 != bin2)", inplace=True)
    a.eval("d = abs(bin2 - bin1)", inplace=True)
    a.query("d <= @max_dist_orientation", inplace=True)
    aa = a.rename(columns={"scaffold_index1": "scaffold_index"})
    assert aa.shape[0] > 0
    aa = aa.groupby("scaffold_index").agg(
        nxt=("pos1", lambda s: s.loc[aa.eval("bin2 > bin1")].mean()),
        prv=("pos1", lambda s: s.loc[aa.eval("bin2 < bin1")].mean()),)
    assert aa.shape[0] > 0
    aa = info[["scaffold_index", "length"]].set_index("scaffold_index").merge(
        aa, left_index=True, right_index=True).reset_index(drop=False)
    assert aa.shape[0] > 0
    # "x != x" means: return the np.nan indices
    idx1 = aa.eval("(prv == prv) & (nxt == nxt)")
    aa.loc[idx1, "orientation"] = aa.loc[idx1].eval("prv <= nxt")
    idx2 = aa.eval("(prv != prv) & (nxt == nxt)")
    aa.loc[idx2, "orientation"] = aa.loc[idx2].eval("length - nxt <= nxt")
    idx3 = aa.eval("(prv == prv) & (nxt != nxt)")
    aa.loc[idx3, "orientation"] = aa.loc[idx3].eval("prv < length - prv")
    aa.loc[:, "orientation"] = aa["orientation"].map({True: 1, False: -1})

    assert "scaffold_index" in membership.columns
    assert "scaffold_index" in aa.columns
    membership = aa[["scaffold_index",
                     "orientation"]].merge(membership, right_on="scaffold_index", left_on="scaffold_index")
    assert "scaffold_index" in membership.columns
    assert membership.shape[0] > 0
    membership.loc[:, "oriented"] = True
    membership.loc[membership["orientation"].isna(), "oriented"] = False
    membership.loc[membership["orientation"].isna(), "orientation"] = 1
    membership.sort_values(["super", "bin", "rank"])
    membership.loc[:, "super_pos"] = membership.groupby("super").transform
    membership.sort_values(["super", "bin", "rank"], inplace=True)
    super_pos = membership.groupby("super")["length"].shift(1, fill_value=0).transform("cumsum") + 1
    membership.loc[:, "super_pos"] = super_pos
    if verbose is True:
        print("Anchoring the super-scaffolds")

    # assign super scaffolds to genetic positions
    y = membership.query("chr == chr").groupby(["chr", "super"]).agg("size").to_frame("nchr")
    y.loc[:, "pchr"] = (y.groupby(level=1)["nchr"].shift(0) / y.groupby(level=1)["nchr"].transform("sum"))
    y = y.reset_index(drop=False)
    y = y.sort_values("nchr", ascending=False).drop_duplicates(subset="super")
    yy = membership.query("cM == cM").merge(y, on=["chr", "super"])
    res = y.merge(res, on="super", how="right")
    left = yy.groupby("super").agg(cM=("cM", "mean"), min_cM=("cM", "min"), max_cM=("cM", "max"))
    res = left.merge(res, left_index=True, right_on="super", how="right").sort_values("length", ascending=False)
    return membership, res
