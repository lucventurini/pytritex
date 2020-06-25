from pytritex.graph_utils.make_super_scaffolds import make_super_scaffolds
from dask.distributed import Client
from . import minimum_distance
import dask.dataframe as dd


def _remove_bulges(links: str, excluded,
                   membership: str,
                   client: Client,
                   save_dir: str,
                   info: str, min_dist=minimum_distance, ncores=1):
    dd_membership = dd.read_parquet(membership, infer_divisions=True)
    add = dd_membership[(dd_membership["rank"] == 1) & (
            dd_membership["length"] <= min_dist)].index.values.compute()
    if add.shape[0] > 0:
        excluded.update(add.tolist())
        out = make_super_scaffolds(links=links, client=client, save_dir=save_dir,
                                   membership=membership,
                                   info=info, excluded=excluded, ncores=ncores)
        membership = out["membership"]
        res = out["info"]
    else:
        res = None
    return membership, res, excluded