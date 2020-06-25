from pytritex.graph_utils.make_super_scaffolds import make_super_scaffolds
from dask.distributed import Client
from . import minimum_distance, _calculate_degree
from .bulge_removal import _remove_bulges
import dask.dataframe as dd
from dask.delayed import delayed
import numpy as np
import pandas as pd


def _remove_short_tips(links: dd.DataFrame,
                       excluded,
                       membership: dd.DataFrame,
                       info: dd.DataFrame,
                       client: Client,
                       save_dir: str,
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

    inner = membership.loc[membership["rank"] == 1, ["super", "bin"]].values.compute()
    inner0 = np.tile(inner[:, 0], 3)
    inner1 = np.tile(inner[:, 1], 3).reshape(3, inner.shape[0]) + np.array([0, -1, 1]).reshape(3, 1)
    inner1 = inner1.flatten()
    inner = pd.DataFrame().assign(super=inner0, bin=inner1)
    right = client.scatter(membership.reset_index(drop=False))
    func = delayed(dd.merge)(right, inner, how="right", on=["super", "bin"])
    middle = client.compute(func).result().set_index("scaffold_index")
    degree = _calculate_degree(links, excluded)
    middle = client.scatter(middle)
    func = delayed(dd.merge)(degree, middle, left_index=True, right_index=True)
    a = client.compute(func).result()
    add = a.query("(degree == 1) & (length <= @min_dist)").index.compute().values
    out = {"membership": membership, "info": None}
    if add.shape[0] > 0:
        excluded.update(set(add.tolist()))
        out = make_super_scaffolds(links=links, info=info, client=client,
                                   save_dir=save_dir,
                                   membership=membership,
                                   excluded=excluded, ncores=ncores)
        membership = out["membership"]

    membership, res, excluded = _remove_bulges(links=links, excluded=excluded,
                                               membership=membership,
                                               client=client, save_dir=save_dir,
                                               info=info, min_dist=min_dist, ncores=ncores)
    if res is not None:
        out["info"] = res

    return membership, out["info"], excluded
