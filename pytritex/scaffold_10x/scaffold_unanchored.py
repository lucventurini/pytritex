import pandas as pd
import numpy as np
from pytritex.graph_utils.make_super_scaffolds import make_super_scaffolds
import time
from dask.distributed import Client
import dask.dataframe as dd
from dask.delayed import delayed
import dask.array as da


def _scaffold_unanchored(links: str,
                         excluded,
                         membership: str,
                         info: str,
                         sample_count: str,
                         save_dir: str,
                         client: Client,
                         ncores=1,
                         verbose=False):

    if verbose:
        print(time.ctime(), "Starting to unanchor unassigned scaffolds")

    sample_count = dd.read_parquet(sample_count, infer_divisions=True)
    links = dd.read_parquet(links, infer_divisions=True)
    membership = dd.read_parquet(membership, infer_divisions=True)
    info = dd.read_parquet(info, infer_divisions=True)

    bait1 = sample_count["popseq_chr1"].isna()
    bait2 = sample_count["popseq_chr2"].isna()

    unanchored = sample_count.loc[(bait1 & ~bait2) | (bait2 & ~bait1)][
        ["scaffold_index1", "length1", "scaffold_index2"]].copy()

    left_side = unanchored.copy()
    left_side.columns = ["scaffold_link", "link_length", "scaffold_index1"]
    left_side = left_side.set_index("scaffold_link")
    right_side = unanchored.copy()[["scaffold_index1", "scaffold_index2"]]
    right_side.columns = ["scaffold_link", "scaffold_index2"]
    right_side = right_side.set_index("scaffold_link")
    left_side = client.scatter(left_side)
    right_side = client.scatter(right_side)
    func = delayed(dd.merge)(left_side, right_side, on="scaffold_link", how="outer")
    linkages = client.compute(func).result().query(
        "scaffold_index1 != scaffold_index2").reset_index(drop=False)
    if linkages.shape[0] == 0 and unanchored.shape[0] > 0:
        raise ValueError("No remaining unanchored scaffolds!")
    # m[, .(scaffold1=scaffold, super1=super, chr1=chr, cM1=cM, size1=super_nbin,
    # d1 = pmin(bin - 1, super_nbin - bin))][xy, on="scaffold1"]->xy
    distances = membership[["super", "chr", "cM", "super_nbin", "bin"]].copy()
    chunks = tuple(distances.map_partitions(len).compute().values.tolist())

    distances["d"] = da.from_array(np.minimum(distances.eval("bin - 1").compute(),
                                              distances.eval("super_nbin - bin").compute()),
                                   chunks=chunks)
    distances = distances.loc[:, ["super", "chr", "cM", "super_nbin", "d"]]
    distances.columns = ["super", "chr", "cM", "size", "d"]
    left_side = distances.copy()
    left_side.columns = ["{}1".format(col) for col in left_side.columns]
    left_side.index = left_side.index.rename("scaffold_index1")
    left_side = client.scatter(left_side)
    linkages = client.scatter(linkages.set_index("scaffold_index1"))
    func = delayed(dd.merge)(left_side, linkages, left_index=True, right_index=True, how="right")
    linkages = client.compute(func).result().reset_index(drop=False)

    left_side = distances.copy()
    left_side.columns = ["{}2".format(col) for col in left_side.columns]
    left_side.index = left_side.index.rename("scaffold_index2")
    left_side = client.scatter(left_side)
    linkages = client.scatter(linkages.set_index("scaffold_index2"))
    func = delayed(dd.merge)(left_side, linkages, left_index=True, right_index=True, how="right")
    linkages = client.compute(func).result().reset_index(drop=False)

    new_linkages = linkages.loc[~((linkages["super2"].isna()) | (linkages["super1"].isna()))].query(
        "(super2 != super1) & (d1 == d2 == 0) & (size1 > 1) & (size2 > 1) & (chr1 == chr2)")

    # Now find and keep those scaffold links that are unambiguous. Use the index to remove the double counting.
    nscl = new_linkages.drop_duplicates(
        ["scaffold_index1", "scaffold_link", "scaffold_index2"]).query(
        "scaffold_index1 < scaffold_index2").groupby(["scaffold_link"]).size().to_frame("nscl")
    nscl = client.scatter(nscl)
    func = delayed(dd.merge)(nscl, new_linkages, on="scaffold_link")
    linkages = client.compute(func).result()
    # xy[super1 < super2][, c("n", "g"):=list(.N, .GRP),
    # by=.(super1, super2)][order(-link_length)][!duplicated(g)]->zz
    keys = ["super1", "super2"]
    grouped = linkages.query("super1 < super2").groupby(keys)
    group_number = grouped.ngroup().to_frame("g").compute()
    
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

    out = make_super_scaffolds(links=links2, save_dir=save_dir, info=info, excluded=excluded,
                               ncores=ncores, client=client)
    membership = out["membership"]
    res = out["info"]
    return membership, res
