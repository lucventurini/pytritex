from pytritex.graph_utils.make_super_scaffolds import make_super_scaffolds
from dask.distributed import Client
from . import minimum_distance, _calculate_degree
from .bulge_removal import _remove_bulges
import dask.dataframe as dd
import numpy as np
import pandas as pd


def _remove_bifurcations(links: str,
                         excluded,
                         membership: str,
                         info: str,
                         client: Client, save_dir: str,
                         min_dist=minimum_distance, ncores=1):
    # resolve length-one-bifurcations at the ends of paths
    dd_membership = dd.read_parquet(membership, infer_divisions=True)
    if dd_membership.shape[0].compute() == 0:
        return membership, info, excluded
    keys = ["super", "super_nbin", "length", "bin"]
    bifurcated = dd_membership.query("(rank > 0) & ((bin == 2) | (super_nbin - 1 == bin) )")
    if bifurcated.shape[0].compute() == 0:
        return bifurcated, info, excluded
    bifurcated = bifurcated[keys]
    bifurcated = bifurcated.eval("type = (bin == 2)").rename(columns={"bin": "bin0"})
    bifurcated = bifurcated[["super", "super_nbin", "type", "length", "bin0"]]
    key = ["super", "bin"]
    # indexed = dd_membership.set_index(key)
    # m[x[type == T,.(super, bin0, bin=1)], on = c("super", "bin")],
    upper = dd.merge(dd_membership,
        bifurcated.loc[bifurcated.type, ["super", "bin0"]].assign(bin=1),
        on=key)
    upper = upper.persist()
    assert upper.index.name == "scaffold_index"
    # m[x[type == F,.(super, bin0, bin=super_nbin)], on = c("super", "bin")]
    lower = dd.merge(dd_membership,
        bifurcated.loc[~bifurcated.type, ["super", "bin0", "super_nbin"]].rename(
            columns={"super_nbin": "bin"}),
        on=key)
    assert lower.index.name == "scaffold_index"
    a = dd.concat([upper, lower]).drop_duplicates().rename(
        columns={"length": "length2", "scaffold_index": "scaffold_index2"})

    # Now merge with "x" to find places to exclude
    a = a.merge(bifurcated, on=["super", "bin0"], how="right")
    add = np.where((a.length >= a.length2), a.scaffold_index2, a.scaffold_index)
    out = {"membership": membership, "info": None}
    if add.shape[0] > 0:
        excluded.update(add.tolist())
        out = make_super_scaffolds(links=links, info=info,
                                   membership=membership,
                                   excluded=excluded, ncores=ncores,
                                   client=client, save_dir=save_dir)
        membership = out["membership"]
    membership, res, excluded = _remove_bulges(links=links, excluded=excluded,
                                               membership=membership,
                                               save_dir=save_dir,
                                               client=client,
                                               info=info, min_dist=min_dist, ncores=ncores)
    if res is not None:
        out["info"] = res
    return membership, out["info"], excluded
