import pandas as pd
from pytritex.graph_utils.make_super_scaffolds import make_super_scaffolds
import dask.dataframe as dd
import numpy as np
from time import ctime
from dask.distributed import Client
from . import minimum_distance, _calculate_degree
from .short_tip_removal import _remove_short_tips


def remove_tips(links: str, excluded, out: dict, info: str,
                client: Client,
                save_dir: str,
                ncores=1,
                verbose=False, min_dist=minimum_distance):

    if verbose:
        print(ctime(), "Starting tip removal")
    # out = {"membership": membership, "info": info}
    membership = out["membership"]
    dd_membership = dd.read_parquet(out["membership"], infer_divisions=True)
    if dd_membership.head(npartitions=-1).shape[0] == 0:
        return out

    membership, res, excluded = _remove_short_tips(links, excluded, membership, info,
                                                   client=client, save_dir=save_dir,
                                                   min_dist=min_dist, ncores=ncores)

    if membership.head(npartitions=-1).shape[0] == 0:
        print("This set of parameters leads to lose everything.")
        return membership, info, excluded
    if res is not None:
        out["info"] = res

    membership, res, excluded = _remove_bifurcations(links, excluded, membership, info,
                                                     client=client, save_dir=save_dir,
                                                     min_dist=min_dist, ncores=ncores)

    if membership.head(npartitions=-1).shape[0] == 0:
        print("This set of parameters leads to lose everything.")
        return membership, info, excluded
    if res is not None:
        out["info"] = res

    degree = _calculate_degree(links, excluded)
    add = degree.merge(
        membership.query("rank == 1"), on="scaffold_index").query("degree == 1")
    add = add.index.values.compute()
    if add.shape[0].compute() > 0:
        excluded.update(add.tolist())
        out = make_super_scaffolds(links=links, info=info, excluded=excluded, ncores=ncores,
                                   client=client, save_dir=save_dir)
        membership = dd.read_parquet(out["membership"], infer_divisions=True)

    add = membership.query("rank > 0")["scaffold_index"]
    if add.shape[0].compute() > 0:
        excluded.update(add.tolist())
        out = make_super_scaffolds(links=links, info=info, excluded=excluded,
                                   membership=membership,
                                   ncores=ncores,
                                   client=client, save_dir=save_dir)

    return membership, out["info"], excluded
