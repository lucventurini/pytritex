import pandas as pd
from pytritex.graph_utils.make_super_scaffolds import make_super_scaffolds
import dask.dataframe as dd
import numpy as np
from dask.delayed import delayed
from time import ctime
from dask.distributed import Client
from . import minimum_distance, _calculate_degree
from .short_tip_removal import _remove_short_tips
from .bifurcation_removal import _remove_bifurcations
import os


def remove_tips(links: str, excluded, out: dict, info: str,
                client: Client,
                save_dir: str,
                ncores=1,
                verbose=False, min_dist=minimum_distance):

    if verbose:
        print(ctime(), "Starting tip removal")
    # out = {"membership": membership, "info": info}
    links = dd.read_parquet(links, infer_divisions=True)
    info = dd.read_parquet(info, infer_divisions=True)
    # membership = out["membership"]
    membership = dd.read_parquet(out["membership"], infer_divisions=True)
    if membership.head(npartitions=-1).shape[0] == 0:
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
    scattered = client.scatter(membership.query("rank == 1"))
    func = delayed(dd.merge)(scattered, degree, on="scaffold_index")
    add = client.compute(func).result().query("degree == 1")
    add = add.index.compute().values
    if add.shape[0] > 0:
        excluded.update(set(add.tolist()))
        out = make_super_scaffolds(links=links, info=info,
                                   membership=membership,
                                   excluded=excluded,
                                   ncores=ncores,
                                   client=client, save_dir=save_dir,
                                   to_parquet=False)

    add = membership.query("rank > 0").index.compute().values
    if add.shape[0] > 0:
        excluded.update(set(add.tolist()))
        out = make_super_scaffolds(links=links, info=info, excluded=excluded,
                                   membership=membership,
                                   ncores=ncores,
                                   client=client, save_dir=save_dir,
                                   to_parquet=False)
    mem_name = os.path.join(save_dir, "membership")
    dd.to_parquet(out["membership"], mem_name, compute=True, compression="gzip", engine="pyarrow")
    res_name = os.path.join(save_dir, "result")
    dd.to_parquet(out["info"], res_name, compute=True, compression="gzip", engine="pyarrow")

    return membership, out["info"], excluded
