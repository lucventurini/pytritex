from pytritex.graph_utils.make_super_scaffolds import make_super_scaffolds
from dask.distributed import Client
from . import minimum_distance, _calculate_degree
from .bulge_removal import _remove_bulges
import dask.dataframe as dd
import numpy as np
from dask.delayed import delayed


def _remove_bifurcations(links: dd.DataFrame,
                         excluded,
                         membership: dd.DataFrame,
                         info: dd.DataFrame,
                         client: Client, save_dir: str,
                         min_dist=minimum_distance, ncores=1):

    if isinstance(membership, str):
        membership = dd.read_parquet(membership, infer_divisions=True)
    else:
        assert isinstance(membership, dd.DataFrame)

    if isinstance(info, str):
        info = dd.read_parquet(info, infer_divisions=True)
    else:
        assert isinstance(info, dd.DataFrame)

    if membership.shape[0].compute() == 0:
        return membership, info, excluded
    keys = ["super", "super_nbin", "length", "bin"]
    bifurcated = membership.query("(rank > 0) & ((bin == 2) | (super_nbin - 1 == bin) )")
    if bifurcated.shape[0].compute() == 0:
        return membership, info, excluded
    bifurcated = bifurcated[keys]
    bifurcated = bifurcated.eval("type = (bin == 2)").rename(columns={"bin": "bin0"})
    bifurcated = bifurcated[["super", "super_nbin", "type", "length", "bin0"]]
    bifurcated = bifurcated.reset_index(drop=False).drop_duplicates()
    bifurcated = bifurcated.astype({"super": np.int})
    key = ["super", "bin"]
    # indexed = dd_membership.set_index(key)
    # m[x[type == T,.(super, bin0, bin=1)], on = c("super", "bin")],
    left = membership.reset_index(drop=False).astype({"super": np.int})
    upper = bifurcated.loc[bifurcated.type, ["super", "bin0"]].assign(bin=1)
    upper = dd.merge(left, upper, on=key)
    upper_size = upper.shape[0].compute()
    assert "scaffold_index" in upper.columns
    # m[x[type == F,.(super, bin0, bin=super_nbin)], on = c("super", "bin")]
    lower = bifurcated.loc[~bifurcated.type, ["super", "bin0", "super_nbin"]].rename(
        columns={"super_nbin": "bin"})
    lower = dd.merge(left, lower, on=key)
    lower_size = lower.shape[0].compute()
    assert max(lower_size, upper_size) > 0
    assert "scaffold_index" in lower.columns
    a = dd.concat([upper, lower])
    assert isinstance(a, dd.DataFrame)
    # assert a.shape[0].compute() > 0
    a = a.drop_duplicates().rename(columns={"length": "length2", "scaffold_index": "scaffold_index2"})
    # Now merge with "x" to find places to exclude
    a = dd.merge(a, bifurcated, on=["super", "bin0"], how="right").compute()
    add = np.where((a.length >= a.length2), a.scaffold_index2, a.scaffold_index)
    out = {"membership": membership, "info": None}
    if add.shape[0] > 0:
        excluded.update(add.tolist())
        out = make_super_scaffolds(links=links, info=info,
                                   membership=membership,
                                   excluded=excluded, ncores=ncores,
                                   client=client, save_dir=save_dir,
                                   to_parquet=False)
        membership = out["membership"]
    membership, res, excluded = _remove_bulges(links=links, excluded=excluded,
                                               membership=membership,
                                               save_dir=save_dir,
                                               client=client,
                                               info=info, min_dist=min_dist, ncores=ncores)
    if res is not None:
        out["info"] = res
    return membership, out["info"], excluded
