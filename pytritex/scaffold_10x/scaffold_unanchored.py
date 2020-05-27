import pandas as pd
import numpy as np
from pytritex.graph_utils.make_super_scaffolds import make_super_scaffolds


def _scaffold_unanchored(links, excluded, membership, info, link_pos, ncores=1,
                         prefix=None, verbose=False):
    # use unanchored scaffolds to link super-scaffolds
    # ww2 is link_pos
    # ww2[is.na(popseq_chr1), .(scaffold_link=scaffold1, link_length=length1, scaffold1=scaffold2)]->x
    # ww2[is.na(popseq_chr1), .(scaffold_link=scaffold1, scaffold2=scaffold2)]->y
    # x[y, on="scaffold_link", allow.cartesian=T][scaffold1 != scaffold2]->xy
    left_unanchored = link_pos.query("popseq_chr1 != popseq_chr1")[["scaffold_index1", "length1", "scaffold_index2"]]
    left_side = left_unanchored[:]
    left_side.columns = ["scaffold_link", "link_length", "scaffold_index1"]
    left_side.set_index("scaffold_link")
    right_side = left_unanchored.loc[:, "scaffold_index1", "scaffold_index2"]
    right_side.columns = ["scaffold_link", "scaffold_index2"]
    right_side = right_side.set_index("scaffold_link")
    linkages = left_side.merge(right_side, on="scaffold_link", how="outer").query(
        "scaffold1 != scaffold2").reset_index(drop=False)

    # m[, .(scaffold1=scaffold, super1=super, chr1=chr, cM1=cM, size1=super_nbin, d1 = pmin(bin - 1, super_nbin - bin))][xy, on="scaffold1"]->xy
    # m[, .(scaffold2=scaffold, super2=super, chr2=chr, cM2=cM, size2=super_nbin, d2 = pmin(bin - 1, super_nbin - bin))][xy, on="scaffold2"]->xy
    distances = membership.eval("d = min(bin - 1, super_nbin - bin)")
    distances = distances[["scaffold_index", "super", "chr", "cM", "super_nbin", "d"]]
    left_side = distances[:].add_suffix("1")
    linkages = left_side.merge(linkages, on=["scaffold_index1"], how="right")
    left_side = distances[:].add_suffix("2")
    #         xy[super2 != super1 & d1 == 0 & d2 == 0 & size1 > 1 & size2 > 1 & chr1 == chr2]->xy
    #         xy[scaffold1 < scaffold2, .(nscl=.N), scaffold_link][xy, on="scaffold_link"]->xy
    #         xy[nscl == 1] -> xy
    linkages = left_side.merge(linkages, on=["scaffold_index1"], how="right")
    linkages = linkages.query("super2 != super1 & d1 == 0 & d2 == 0 & size1 > 1 & size2 > 1 & chr1 == chr2")
    nscl = linkages.query("scaffold_index1 < scaffold_index2").groupby("scaffold_link").size().to_frame("nscl")
    linkages = nscl.merge(linkages, on="scaffold_link").query("nscl == 1")
    # xy[super1 < super2][, c("n", "g"):=list(.N, .GRP), by=.(super1, super2)][order(-link_length)][!duplicated(g)]->zz
    keys = ["super1", "super2"]
    grouped = linkages.query("super1 < super2").groupby(keys)
    group_number = grouped.ngroup().to_frame("g")
    link_count = linkages.eval("g = @group_number").merge(grouped.size().to_frame("n"), on=keys).sort_values(
        "link_length", ascending=False).drop_duplicates("g")
    if "super1" not in link_count.columns:
        # Reset the index
        link_count = link_count.reset_index(drop=False)
    _scaffold_link = link_count["scaffold_link"].values
    _scaffold1 = link_count["scaffold_index1"].values
    _scaffold2 = link_count["scaffold_index2"].values
    sel = pd.DataFrame().assign(
        scaffold_index1=np.vstack([_scaffold_link, _scaffold_link, _scaffold1, _scaffold2]),
        scaffold_index2=np.vstack([_scaffold1, _scaffold2, _scaffold_link, _scaffold_link])
    )
    links2 = pd.concat(
        [links,
         link_pos.merge(sel, how="right", on=["scaffold_index1", "scaffold_index2"])])

    out = make_super_scaffolds(links=links2, prefix=prefix, info=info, excluded=excluded,
                               ncores=ncores)
    membership = out["membership"]
    res = out["info"]
    return membership, res
