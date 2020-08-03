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

    if isinstance(membership, str):
        membership = dd.read_parquet(membership, infer_divisions=True)
    else:
        assert isinstance(membership, dd.DataFrame)

    assert isinstance(min_dist, (float, int)), (min_dist, type(min_dist))
    inner = membership.loc[membership["rank"] == 1, ["super", "bin"]].values.compute()
    if inner.shape[0] > 0:
        inner0 = np.tile(inner[:, 0], 3)
        inner1 = np.tile(inner[:, 1], 3).reshape(3, inner.shape[0]) + np.array([0, -1, 1]).reshape(3, 1)
        inner1 = inner1.flatten()
        inner = pd.DataFrame().assign(super=inner0, bin=inner1).astype({"super": np.int})
        right = membership.reset_index(drop=False).astype({"super": np.int})
        middle = dd.merge(right, inner, how="right", on=["super", "bin"]).set_index("scaffold_index").compute()
        assert isinstance(middle, pd.DataFrame), (type(middle),)
        degree = _calculate_degree(links, excluded)
        a = pd.merge(degree, middle, left_index=True, right_index=True)
        add = a.loc[(a["degree"] == 1) & (a["length"] <= min_dist)]
        assert isinstance(add, pd.DataFrame), (type(add),)
        add = add.index.values
        out = {"membership": membership, "info": None}
        if add.shape[0] > 0:
            excluded.update(set(add.tolist()))
            out = make_super_scaffolds(links=links, info=info, client=client,
                                   save_dir=save_dir,
                                   membership=membership,
                                   excluded=excluded, ncores=ncores)
            membership = out["membership"]
    else:
        out = {"membership": membership, "info": info}

    membership, res, excluded = _remove_bulges(links=links, excluded=excluded,
                                               membership=membership,
                                               client=client, save_dir=save_dir,
                                               info=info, min_dist=min_dist, ncores=ncores)
    if res is not None:
        out["info"] = res

    return membership, out["info"], excluded
