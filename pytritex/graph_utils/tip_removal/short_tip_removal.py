from pytritex.graph_utils.make_super_scaffolds import make_super_scaffolds
from dask.distributed import Client
from . import minimum_distance, _calculate_degree
from .bulge_removal import _remove_bulges
import dask.dataframe as dd
import numpy as np
import pandas as pd


def _remove_short_tips(links: str, excluded, membership: str, info: str,
                       client: Client, save_dir: str,
                       min_dist=minimum_distance, ncores=1):
    # #   links[!scaffold2 %in% ex][, .(degree=.N), key=.(scaffold_index=scaffold1)]-> degree
    # # inner <- m[rank == 1][, .(super=c(super, super, super), bin=c(bin, bin-1, bin+1))]
    # # middle <- m[inner, on=c("super", "bin")]
    # # a <- degree[middle, on="scaffold_index"]
    # #  add <-  a[degree == 1 & length <= 1e4]$scaffold
    # #   if(length(add) > 0){
    # #    ex <- c(ex, add)
    # #    make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
    # #    out$m -> m
    # #   }
    # #   # remove short tips/bulges of rank 1
    # #   m[rank == 1 & length <= 1e4]$scaffold -> add
    # #   ex <- c(ex, add)
    # #   if(length(add) > 0){
    # #    ex <- c(ex, add)
    # #    make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
    # #    out$m -> m
    # #   }

    dd_membership = dd.read_parquet(membership, infer_divisions=True)
    inner = dd_membership.loc[dd_membership["rank"] == 1, ["super", "bin"]].values.compute()
    inner0 = np.tile(inner[:, 0], 3)
    inner1 = np.tile(inner[:, 1], 3).reshape(3, inner.shape[0]) + np.array([0, -1, 1]).reshape(3, 1)
    inner1 = inner1.flatten()
    inner = pd.DataFrame().assign(super=inner0, bin=inner1)
    middle = dd.merge(membership, inner, how="right",
                      on=["super", "bin"]).reset_index(drop=False).set_index("scaffold_index")
    links = dd.read_parquet(links, infer_divisions=True)
    degree = _calculate_degree(links, excluded)
    a = degree.merge(middle, on="scaffold_index").reset_index(drop=False)
    add = a.query("(degree == 1) & (length <= @min_dist)")["scaffold_index"].values
    out = {"membership": membership, "info": None}
    if add.shape[0] > 0:
        excluded.update(add.tolist())
        out = make_super_scaffolds(links=links, info=info, client=client,
                                   save_dir=save_dir,
                                   membership=membership,
                                   excluded=excluded, ncores=ncores)
        membership = out["membership"]

    info = dd.read_parquet(info, infer_divisions=True)
    membership, res, excluded = _remove_bulges(links=links, excluded=excluded,
                                               membership=membership,
                                               client=client, save_dir=save_dir,
                                               info=info, min_dist=min_dist, ncores=ncores)
    if res is not None:
        out["info"] = res

    return membership, out["info"], excluded