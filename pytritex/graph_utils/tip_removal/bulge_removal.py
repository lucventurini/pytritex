from pytritex.graph_utils.make_super_scaffolds import make_super_scaffolds
from dask.distributed import Client
from . import minimum_distance
import dask.dataframe as dd


def _remove_bulges(links: dd.DataFrame,
                   excluded,
                   membership: dd.DataFrame,
                   client: Client,
                   save_dir: str,
                   info: dd.DataFrame, min_dist=minimum_distance, ncores=1):

    add = membership[(membership["rank"] == 1) & (membership["length"] <= min_dist)].index.compute().values
    if add.shape[0] > 0:
        excluded.update(add.tolist())
        out = make_super_scaffolds(links=links, client=client, save_dir=save_dir,
                                   membership=membership,
                                   info=info, excluded=excluded, ncores=ncores,
                                   to_parquet=False)
        membership = out["membership"]
        res = out["info"]
    else:
        res = None
    return membership, res, excluded
