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
import logging
dask_logger = logging.getLogger("dask")


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
    if membership.shape[0].compute() == 0:
        return out

    dask_logger.warning("%s Removing short tips", ctime())
    membership, res, excluded = _remove_short_tips(links, excluded, membership, info,
                                                   client=client, save_dir=save_dir,
                                                   min_dist=min_dist, ncores=ncores)
    dask_logger.warning("%s Removed short tips", ctime())    

    if membership.head(npartitions=-1).shape[0] == 0:
        dask_logger.warning("%s This set of parameters leads to lose everything.", ctime())
        return membership, info, excluded
    if res is not None:
        out["info"] = res

    dask_logger.warning("%s Removing bifurcations", ctime())        
    membership, res, excluded = _remove_bifurcations(links, excluded, membership, info,
                                                     client=client, save_dir=save_dir,
                                                     min_dist=min_dist, ncores=ncores)
    dask_logger.warning("%s Removed bifurcations", ctime())

    if isinstance(membership, str):
        membership = dd.read_parquet(membership, infer_divisions=True)
    else:
        assert isinstance(membership, dd.DataFrame)

    if isinstance(res, str):
        res = dd.read_parquet(res, infer_divisions=True)
    else:
        assert isinstance(res, dd.DataFrame) or res is None

    if membership.shape[0].compute() == 0:
        dask_logger.warning("%s This set of parameters leads to lose everything.", ctime())
        return membership, info, excluded
    if res is not None:
        out["info"] = res

    dask_logger.warning("%s Removing the last tips", ctime())
    degree = _calculate_degree(links, excluded)
    scattered = membership.query("rank == 1")
    add = dd.merge(scattered, degree, on="scaffold_index").query("degree == 1")
    add = add.index.compute().values
    if add.shape[0] > 0:
        excluded.update(set(add.tolist()))
        out = make_super_scaffolds(links=links, info=info,
                                   membership=membership,
                                   excluded=excluded,
                                   ncores=ncores,
                                   client=client, save_dir=save_dir,
                                   to_parquet=False)
        membership = out["membership"]
        if isinstance(membership, str):
            membership = dd.read_parquet(membership, infer_divisions=True)
        else:
            assert isinstance(membership, dd.DataFrame)

    dask_logger.warning("%s Removing any stray non-0 rank scaffolds.", ctime())            
    add = membership.query("rank > 0").index.compute().values
    if add.shape[0] > 0:
        excluded.update(set(add.tolist()))
        out = make_super_scaffolds(links=links, info=info, excluded=excluded,
                                   membership=membership,
                                   ncores=ncores,
                                   client=client, save_dir=save_dir,
                                   to_parquet=False)
        if isinstance(membership, str):
            membership = dd.read_parquet(membership, infer_divisions=True)
        else:
            assert isinstance(membership, dd.DataFrame)
        
    mem_name = os.path.join(save_dir, "membership")
    dd.to_parquet(out["membership"], mem_name, compute=True, compression="gzip", engine="pyarrow")
    res_name = os.path.join(save_dir, "result")
    dd.to_parquet(out["info"], res_name, compute=True, compression="gzip", engine="pyarrow")
    dask_logger.warning("%s Finished tip removal.", ctime())
    
    return mem_name, res_name, excluded
