import pandas as pd
import numpy as np
from pytritex.graph_utils.make_super_scaffolds import make_super_scaffolds
import time


def _scaffold_unanchored(links, excluded, membership, info, sample_count, ncores=1,
                         prefix=None, verbose=False):

    if verbose:
        print(time.ctime(), "Starting to unanchor unassigned scaffolds")

    bait1 = sample_count["popseq_chr1"].isna()
    bait2 = sample_count["popseq_chr2"].isna()

    unanchored = sample_count.loc[(bait1 & ~bait2) | (bait2 & ~bait1)][
        ["scaffold_index1", "length1", "scaffold_index2"]].copy()

    left_side = unanchored.copy()
    left_side.columns = ["scaffold_link", "link_length", "scaffold_index1"]
    left_side.set_index("scaffold_link")
    right_side = unanchored.copy()[["scaffold_index1", "scaffold_index2"]]
    right_side.columns = ["scaffold_link", "scaffold_index2"]
    right_side = right_side.set_index("scaffold_link")
    linkages = left_side.merge(right_side, on="scaffold_link", how="outer").query(
        "scaffold_index1 != scaffold_index2").reset_index(drop=False)
    if linkages.shape[0] == 0 and unanchored.shape[0] > 0:
        raise ValueError("No remaining unanchored scaffolds!")
    # m[, .(scaffold1=scaffold, super1=super, chr1=chr, cM1=cM, size1=super_nbin,
    # d1 = pmin(bin - 1, super_nbin - bin))][xy, on="scaffold1"]->xy
    distances = membership[["scaffold_index", "super", "chr", "cM", "super_nbin", "bin"]].copy()
    distances.loc[:, "d"] = np.minimum(distances.eval("bin - 1"), distances.eval("super_nbin - bin"))
    distances = distances.loc[:, ["scaffold_index", "super", "chr", "cM", "super_nbin", "d"]]
    distances.columns = ["scaffold_index", "super", "chr", "cM", "size", "d"]
    distances = distances.set_index("scaffold_index")
    print("Distances", distances.shape[0])
    left_side = distances[:].add_suffix("1")
    left_side.index.names = ["scaffold_index1"]
    print(linkages.shape[0])
    linkages = linkages.set_index("scaffold_index1")
    linkages = left_side.merge(linkages,
                               on=["scaffold_index1"], how="right").reset_index(drop=False)
    print("Merged on scaffold_index1", linkages.shape[0])
    # m[, .(scaffold2=scaffold, super2=super, chr2=chr, cM2=cM, size2=super_nbin,
    # d2 = pmin(bin - 1, super_nbin - bin))][xy, on="scaffold2"]->xy
    left_side = distances[:].add_suffix("2")
    left_side.index.names = ["scaffold_index2"]
    linkages = left_side.merge(linkages.set_index("scaffold_index2"),
                               on=["scaffold_index2"], how="right").reset_index(drop=False)
    print("Merged on scaffold_index2", linkages.shape[0])
    print(linkages.d1.min())
    print(linkages.d2.min())
    print(linkages.size1.max())
    print(linkages.size2.max())
    print()
    print(linkages.query("chr1 == chr2").shape[0])

    new_linkages = linkages.loc[~((linkages["super2"].isna()) | (linkages["super1"].isna()))].query(
        "(super2 != super1) & (d1 == d2 == 0) & (size1 > 1) & (size2 > 1) & (chr1 == chr2)")

    print("Filtered", new_linkages.shape[0])
    # Now find and keep those scaffold links that are unambiguous. Use the index to remove the double counting.
    nscl = new_linkages.drop_duplicates(
        ["scaffold_index1", "scaffold_link", "scaffold_index2"]).query(
        "scaffold_index1 < scaffold_index2").groupby(["scaffold_link"]).size().to_frame("nscl")

    linkages = nscl.merge(new_linkages, on="scaffold_link")

    # xy[super1 < super2][, c("n", "g"):=list(.N, .GRP),
    # by=.(super1, super2)][order(-link_length)][!duplicated(g)]->zz
    keys = ["super1", "super2"]
    grouped = linkages.query("super1 < super2").groupby(keys)
    group_number = grouped.ngroup().to_frame("g")
    link_count = linkages.eval("g = @group_number").merge(grouped.size().to_frame("n"), on=keys)
    link_count = link_count.sort_values("link_length", ascending=False).drop_duplicates("g")
    if "super1" not in link_count.columns:
        # Reset the index
        link_count = link_count.reset_index(drop=False)
    _scaffold_link = link_count["scaffold_link"].values.flatten()
    _scaffold1 = link_count["scaffold_index1"].values.flatten()
    _scaffold2 = link_count["scaffold_index2"].values.flatten()
    sel = pd.DataFrame().assign(
        scaffold_index1=np.concatenate([_scaffold_link, _scaffold_link, _scaffold1, _scaffold2]),
        scaffold_index2=np.concatenate([_scaffold1, _scaffold2, _scaffold_link, _scaffold_link]))
    lower = sample_count.merge(sel, how="right", on=["scaffold_index1", "scaffold_index2"])
    print(sel.head())
    print(sel.shape)
    print(lower.head())
    print(lower.shape)
    links2 = pd.concat([links, lower]).reset_index(drop=False)
    links2.drop("index", axis=1, inplace=True, errors="ignore")
    links2.drop("cidx", axis=1, inplace=True, errors="ignore")

    out = make_super_scaffolds(links=links2, prefix=prefix, info=info, excluded=excluded,
                               ncores=ncores)
    membership = out["membership"]
    res = out["info"]
    return membership, res
