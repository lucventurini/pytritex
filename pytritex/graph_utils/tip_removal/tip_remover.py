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
from typing import Union
dask_logger = logging.getLogger("dask")


def remove_tips(links: Union[str, pd.DataFrame, dd.DataFrame],
                excluded: Union[None, set], out: dict,
                info: Union[str, pd.DataFrame, dd.DataFrame],
                client: Client,
                save_dir: str,
                ncores=1,
                verbose=False, min_dist=minimum_distance):

    if verbose:
        print(ctime(), "Starting tip removal")
    # out = {"membership": membership, "info": info}
    if isinstance(links, str):
        links = dd.read_parquet(links, infer_divisions=True)
    if isinstance(info, str):
        info = dd.read_parquet(info, infer_divisions=True)
    # membership = out["membership"]
    membership = out["membership"]
    if isinstance(membership, str):
        membership = dd.read_parquet(out["membership"], infer_divisions=True)
    if membership.shape[0].compute() == 0:
        return out

    dask_logger.warning("%s Removing short tips", ctime())
    new_membership, new_res, excluded = _remove_short_tips(links, excluded, membership, info,
                                                   client=client, save_dir=save_dir,
                                                   min_dist=min_dist, ncores=ncores)

    dask_logger.warning("%s Removed short tips", ctime())    

    if new_membership.head(npartitions=-1).shape[0] == 0:
        dask_logger.warning("%s This set of parameters leads to lose everything.", ctime())
        return membership, info, excluded
    if new_res is not None:
        out["info"] = new_res
        out["membership"] = new_membership
        membership = new_membership
        info = out["info"]

    dask_logger.warning("%s Removing bifurcations", ctime())        
    new_membership, new_res, excluded = _remove_bifurcations(
        links, excluded, out["membership"], out["info"],
        client=client, save_dir=save_dir,
        min_dist=min_dist, ncores=ncores)
    dask_logger.warning("%s Removed bifurcations", ctime())
    
    if isinstance(new_membership, str):
        new_membership = dd.read_parquet(new_membership, infer_divisions=True)
    else:
        assert isinstance(new_membership, dd.DataFrame)

    if isinstance(new_res, str):
        new_res = dd.read_parquet(new_res, infer_divisions=True)
    else:
        assert isinstance(new_res, dd.DataFrame) or new_res is None

    if new_membership.shape[0].compute() == 0:
        dask_logger.warning("%s This set of parameters leads to lose everything.", ctime())
        return membership, info, excluded
    if new_res is not None:
        out["info"] = new_res
        out["membership"] = new_membership
        info = out["info"]
        membership = new_membership

    dask_logger.warning("%s Removing the last tips", ctime())
    degree = _calculate_degree(links, excluded)
    scattered = membership.query("rank == 1")
    add = dd.merge(scattered, degree, on="scaffold_index").query("degree == 1")
    add = add.index.compute().values
    membership = out["membership"].copy()
    if add.shape[0] > 0:
        excluded.update(set(add.tolist()))
        new_out = make_super_scaffolds(links=links, info=out["info"],
                                   membership=out["membership"],
                                   excluded=excluded,
                                   ncores=ncores,
                                   client=client, save_dir=save_dir,
                                   to_parquet=False)
        if isinstance(new_out["membership"], str):
            new_membership = dd.read_parquet(new_out["membership"], infer_divisions=True)
        else:
            new_membership = new_out["membership"]
            assert isinstance(new_membership, dd.DataFrame)
        if new_out["info"] is not None and new_membership.shape[0].compute() > 0:
            out["info"] = new_res
            out["membership"] = new_membership
            info = out["info"]
            membership = new_membership

    dask_logger.warning("%s Removing any stray non-0 rank scaffolds.", ctime())            
    add = membership.query("rank > 0").index.compute().values
    if add.shape[0] > 0:
        new_excluded = {}
        new_excluded.update(excluded)
        new_excluded.update(set(add.tolist()))
        new_out = make_super_scaffolds(links=links, info=out["info"], excluded=new_excluded,
                                   membership=out["membership"],
                                   ncores=ncores,
                                   client=client, save_dir=save_dir,
                                   to_parquet=False)
        if isinstance(new_out["membership"], str):
            new_membership = dd.read_parquet(new_out["membership"], infer_divisions=True)
        else:
            new_membership = new_out["membership"]
            assert isinstance(new_membership, dd.DataFrame)
        if new_out["info"] is not None and new_membership.shape[0].compute() > 0:
            out["info"] = new_res
            out["membership"] = new_membership
            info = out["info"]
            membership = new_membership
            excluded = new_excluded

    mem_name = os.path.join(save_dir, "membership")
    if isinstance(out["membership"], dd.DataFrame):  # We have removed some stuff.
        dd.to_parquet(out["membership"], mem_name, compute=True, compression="gzip", engine="pyarrow",
                      schema="infer")
    res_name = os.path.join(save_dir, "result")
    dd.to_parquet(out["info"], res_name, compute=True, compression="gzip", engine="pyarrow", schema="infer")
    dask_logger.warning("%s Finished tip removal.", ctime())
    
    return mem_name, res_name, excluded
