from pytritex.graph_utils.make_super_scaffolds import make_super_scaffolds, add_missing_scaffolds, add_statistics
import time
from dask.distributed import Client
import dask.dataframe as dd
import logging
import os
from functools import partial
import numpy as np
logger = logging.getLogger("distributed.worker")


def iteration(counter, membership, excluded, links, save_dir, client, info, ncores):
    run = True
    logger.warning("%s Starting run %s , excluded: %s",
                   time.ctime(), counter, len(excluded))
    out = make_super_scaffolds(links=links, save_dir=save_dir,
                               client=client, membership=membership,
                               info=info, excluded=excluded, ncores=ncores,
                               to_parquet=False)
    membership = out["membership"]
    # dd_membership = dd.read_parquet(membership, infer_divisions=True)
    # Now we have to exclude from consideration those scaffolds
    # that are the backbone of super-scaffolds where there is at least a scaffold with
    # rank > 1.
    # Ie: remove the backbones so that the branches can be reassigned to somewhere else.
    a = membership.merge(
        membership.loc[membership["rank"] > 1, ["super", "bin"]].drop_duplicates(),
        how="inner", on=["super", "bin"])
    add = a.loc[a["rank"] == 0, :]
    if add.shape[0].compute() == 0:
        run = False
    else:
        # Drop all the links between the backbone of the "fuzzy" scaffolds and the spikes.
        links = links[
            (links["scaffold_index1"].isin(add) & (links["scaffold_index2"].isin(add))) |
            (~links["scaffold_index1"].isin(add) & (~links["scaffold_index2"].isin(add)))
        ]
        links = links.persist()
    logger.warning("Finished run %s", counter)
    return out, links, run


def _initial_branch_remover(client: Client,
                            save_dir: str,
                            links: str,
                            info: str, excluded: set, ncores):

    print(time.ctime(), "Starting the run")
    if excluded is None:
        excluded = set()

    links = dd.read_parquet(links, infer_divisions=True)
    info = dd.read_parquet(info, infer_divisions=True)
    scaffolds_to_use = np.unique(links[["scaffold_index1", "scaffold_index2"]].values.compute().flatten())
    info_to_use = info.loc[scaffolds_to_use].persist()

    _iterator = partial(iteration,
                        excluded=excluded,
                        save_dir=save_dir,
                        client=client,
                        info=info_to_use,
                        ncores=ncores)
    counter = 1
    out, links, run = _iterator(counter=counter, membership=None, excluded=excluded)
    membership = out["membership"]
    while run is True:
        counter += 1
        out, excluded, run = _iterator(counter=counter,
                                       membership=membership,
                                       links=links)
        membership = out["membership"]

    # Now we need to rejoin things
    maxidx = out["membership"]["super"].max().compute()
    out["membership"] = add_missing_scaffolds(info, out["membership"],
                                              maxidx, excluded, client)
    out["membership"], out["info"] = add_statistics(out["membership"], client)

    dd.to_parquet(out["membership"], os.path.join(save_dir, "membership"))
    # res = dd.from_pandas(res, chunksize=1000)
    dd.to_parquet(out["info"], os.path.join(save_dir, "result"))
    out = {"membership": os.path.join(save_dir, "membership"),
           "info": os.path.join(save_dir, "result")}

    return out, excluded
