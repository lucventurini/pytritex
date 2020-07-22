import os
import pandas as pd
import numpy as np
import dask.dataframe as dd
import dask.array as da
from dask.delayed import delayed
from dask.distributed import Client
import logging
logger = logging.getLogger("distributed.worker")
import time
from ..graph_utils.make_super_scaffolds import add_missing_scaffolds, add_statistics
# from .make_agp import make_agp


def super_position(group):
    orig_index = group.index[:]
    group = group.reset_index(drop=True)
    indexer = group.index[:]
    group = group.sort_values(["bin", "rank"])
    # return group[["length"]].values
    vals = (group["length"].shift(1, fill_value=0).cumsum() + 1).astype(np.int)
    s = vals.shape[0]
    vals = vals[indexer]
    # assert s == vals.shape[0]
    # vals.index = orig_index
    vals = np.vstack([group["scaffold_index"], vals])
    return vals


def chrom_percentage(group):
    chr_sum = group["nchr"].sum()
    return (group["nchr"] / chr_sum).values


def _create_association(membership, info, link_pos, client, max_dist_orientation):
    logger.warning("%s Getting scaffolds in super-scaffolds", time.ctime())
    m_greater_one = membership.query("super_nbin > 1")
    assert m_greater_one.shape[0].compute() > 0, membership.head()
    left = m_greater_one.loc[:, ["bin", "super"]].rename(columns={"bin": "bin1", "super": "super1"})
    left.index = left.index.rename("scaffold_index1")
    association = dd.merge(left, link_pos, on="scaffold_index1", how="inner")
    left = m_greater_one[["bin", "super"]].rename(columns={"bin": "bin2", "super": "super2"})
    left.index = left.index.rename("scaffold_index2")
    association = dd.merge(left, association, on="scaffold_index2", how="inner")
    logger.warning("%s Obtained the association table", time.ctime())
    association = association.query("(super1 == super2) & (bin1 != bin2)")
    association = association.eval("d = abs(bin2 - bin1)")
    association = association[association["d"] <= max_dist_orientation]
    association = association.rename(columns={"scaffold_index1": "scaffold_index"})
    assert association.shape[0].compute() > 0
    logger.warning("%s Grouping the association table by SI", time.ctime())
    grouped = association.groupby("scaffold_index")
    nxt = grouped.apply(lambda group: group.loc[group.eval("bin2 > bin1"), "pos1"].mean(), meta=float).to_frame("nxt")
    prev = grouped.apply(lambda group: group.loc[group.eval("bin2 < bin1"), "pos1"].mean(), meta=float).to_frame("prv")
    final_association = dd.concat([prev, nxt], axis=1)
    final_association = info[["length"]].merge(final_association, left_index=True, right_index=True)
    return final_association


def _calculate_orientation_column(membership, final_association, client):
    logger.warning("%s Calculating the orientation variable for scaffolds", time.ctime())
    idx1 = final_association.eval("(prv == prv) & (nxt == nxt)")
    index = final_association.index.compute()
    orientation = pd.Series([np.nan] * index.shape[0], index=index)
    orientation.loc[idx1.compute()] = final_association.loc[idx1].eval("prv <= nxt").compute()
    idx2 = final_association.eval("(prv != prv) & (nxt == nxt)")
    orientation.loc[idx2.compute()] = final_association.loc[idx2].eval("length - nxt <= nxt").compute()
    idx3 = final_association.eval("(prv == prv) & (nxt != nxt)")
    orientation.loc[idx3.compute()] = final_association.loc[idx3].eval("prv < length - prv").compute()
    orientation = orientation.map({True: 1, False: -1})
    chunks = tuple(final_association.map_partitions(len).compute().values.tolist())
    orientation = da.from_array(orientation, chunks=chunks)
    final_association["orientation"] = orientation
    assert isinstance(final_association, dd.DataFrame), type(final_association)
    logger.warning("%s Merging the orientation variable back into membership", time.ctime())
    membership = dd.merge(membership, final_association[["orientation"]],
                          left_index=True, right_index=True, how="left").drop_duplicates()
    logger.warning("%s Merged the orientation variable back into membership", time.ctime())
    return membership


def _calculate_oriented_column(membership):
    # Now assign an "orientation" value to each (1 if + or unknown, -1 for -)
    # and an "oriented" value (strictly boolean)
    chunks = tuple(membership.map_partitions(len).compute().values.tolist())
    orientation = membership["orientation"].compute()
    logger.warning("%s Got the orientation column", time.ctime())
    oriented = pd.Series([True] * orientation.shape[0], index=orientation.index, dtype=bool)
    logger.warning("%s Created the oriented column", time.ctime())
    idx1 = orientation.isna().values
    logger.warning("%s Calculated the idx1 variable", time.ctime())
    oriented.loc[idx1] = False
    orientation.loc[idx1] = 1
    logger.warning("%s Creating the dask arrays", time.ctime())
    orientation = da.from_array(orientation, chunks=chunks)
    oriented = da.from_array(oriented, chunks=chunks)
    membership["orientation"] = orientation
    membership["oriented"] = oriented
    logger.warning("%s Calculated the orientation/oriented columns", time.ctime())
    return membership


def orient_scaffolds(info: str, res: str,
                     membership: str,
                     link_pos: str,
                     client: Client,
                     max_dist_orientation: float,
                     save_dir: str):

    membership = dd.read_parquet(membership, infer_divisions=True, engine="auto")
    link_pos = dd.read_parquet(link_pos, infer_divisions=True, engine="auto")
    info = dd.read_parquet(info, infer_divisions=True, engine="auto")
    res = dd.read_parquet(res, infer_divisions=True, engine="auto")

    membership = membership.drop_duplicates()
    assert membership.shape[0].compute() > 0
    maxidx = membership["super"].values.max().compute()
    # #  m[super_nbin > 1, .(scaffold1=scaffold, bin1=bin, super1=super)][link_pos, on="scaffold1", nomatch=0]->a

    final_association = _create_association(membership, info, link_pos, client, max_dist_orientation)
    membership = _calculate_orientation_column(membership, final_association, client)
    membership = _calculate_oriented_column(membership)

    # Next step: assign to each scaffold a base-pair position in the super-scaffold, "super_pos"
    # This position should be determined by bin and rank.

    logger.warning("%s Calculating the positions within super-scaffolds", time.ctime())
    # excluded = membership[(membership["excluded"] == True) | (membership.super_size == 1)]
    # excluded["super_pos"] = 1
    # chunks = tuple(excluded.map_partitions(len).compute().values.tolist())
    # excluded["super_pos"] = np.nan
    non_excluded = membership[(membership["excluded"] == False) & (membership.super_size > 1)]
    grouped = non_excluded.reset_index(drop=False).groupby("super")
    super_pos = grouped.apply(super_position, meta=np.int).compute()
    indices = np.concatenate([_[0, :] for _ in super_pos.values])
    pos = np.concatenate([_[1, :] for _ in super_pos.values])
    super_pos = pd.DataFrame().assign(scaffold_index=indices, super_pos=pos).astype({"scaffold_index": np.int,
                                                                                     "super_pos": np.float}).set_index("scaffold_index")
    super_pos = dd.from_pandas(super_pos, npartitions=min(10, super_pos.shape[0]))
    assert membership.index.name == "scaffold_index", (membership.index.name,)
    membership = membership.reset_index(drop=False).astype({"scaffold_index": np.int}).set_index("scaffold_index")
    membership = dd.merge(membership, super_pos, on="scaffold_index", how="left")
    membership["super_pos"] = membership["super_pos"].fillna(1)
    # Now add back missing scaffolds
    excluded_scaffolds = membership[membership["excluded"] == True].index.values.compute()
    membership = add_missing_scaffolds(info, membership, maxidx, excluded_scaffolds, client, save_dir)
    logger.warning("%s Added the super_pos columns", time.ctime())

    # Now we have to assign the genetic positions, ie anchor
    anchored = membership[~membership["chr"].isna()]
    logger.warning("%s Calculating nchr", time.ctime())
    nchr = anchored.groupby(["chr", "super"])["length"].size().to_frame("nchr").reset_index(drop=False)
    nchr_sum = nchr.groupby("super")["nchr"].transform("sum", meta=np.int)
    logger.warning("%s Calculating pchr", time.ctime())
    nchr["pchr"] = da.from_array((nchr["nchr"] / nchr_sum).compute(),
                                 chunks=tuple(nchr.map_partitions(len).compute().values.tolist()))
    nchr["max_nchr"] = nchr.groupby("super")["nchr"].transform("max", meta=np.int).values
    logger.warning("%s Dropping max_nchr duplicates", time.ctime())
    nchr = nchr.query("nchr == max_nchr").drop_duplicates(subset="super").drop("max_nchr", axis=1)
    # chr, super, nchr, pchr - no index
    logger.warning("%s Merging into nchr_with_cms", time.ctime())
    nchr_with_cms = dd.merge(membership[~membership["cM"].isna()].reset_index(drop=False),
                             nchr, on=["chr", "super"]).set_index("scaffold_index")
    logger.warning("%s Merged into nchr_with_cms", time.ctime())
    grouped_nchr_with_cms = nchr_with_cms.groupby("super")
    logger.warning("%s Calculating the nchr_with_cms stats", time.ctime())
    # m[!is.na(chr), .(nchr=.N), key=.(chr, super)][,
    # pchr := nchr/sum(nchr), by=super][order(-nchr)][!duplicated(super)]->y
    #  m[!is.na(cM)][y, on=c("chr", "super")]->yy
    #  y[res, on="super"]->res
    #  yy[, .(cM=mean(cM), min_cM=min(cM), max_cM=max(cM)), key=super][res, on='super']->res
    #  setorder(res, -length)
    min_cM = grouped_nchr_with_cms["cM"].apply(min, meta=np.float).to_frame("min_cM").compute()
    max_cM = grouped_nchr_with_cms["cM"].apply(max, meta=np.float).to_frame("max_cM").compute()
    mean_cM = grouped_nchr_with_cms["cM"].apply(np.mean, meta=np.float).to_frame("cM").compute()
    nchr_with_cms = pd.merge(mean_cM, pd.merge(min_cM, max_cM, on="super"), on="super")
    logger.warning("%s Calculated the nchr_with_cms stats", time.ctime())
    if res.index.name is None:
        res = res.set_index("super")
    nchr = nchr.set_index("super")
    res = dd.merge(nchr, res, on="super", how="right")
    res = dd.merge(nchr_with_cms, res, on="super", how="right")

    logger.warning("%s Merged everything into res, saving", time.ctime())
    res_name = os.path.join(save_dir, "res")
    dd.to_parquet(res, res_name, compression="gzip", compute=True, engine="pyarrow")
    mem_name = os.path.join(save_dir, "membership")
    dd.to_parquet(membership, mem_name, compression="gzip", compute=True, engine="pyarrow")
    logger.warning("%s Finished orient scaffolds", time.ctime())
    return mem_name, res_name
