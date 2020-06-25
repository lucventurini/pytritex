from pytritex.graph_utils.make_super_scaffolds import make_super_scaffolds
import time
from joblib import Memory
from dask.distributed import Client
import dask.dataframe as dd
import logging
logger = logging.getLogger("distributed.worker")


def _initial_branch_remover(client: Client,
                            save_dir: str,
                            links: str,
                            info: str, excluded: set, ncores):

    print(time.ctime(), "Starting the run")
    if excluded is None:
        excluded = set()

    membership = None

    run = True
    counter = 0
    while run is True:
        counter += 1
        # print(time.ctime(), "Starting run", counter)
        logger.warning("%s Starting run %s , excluded: %s",
                       time.ctime(), counter, len(excluded))
        out = make_super_scaffolds(links=links, save_dir=save_dir,
            client=client, membership=membership,
            info=info, excluded=excluded, ncores=ncores)
        membership = dd.read_parquet(out["membership"], infer_divisions=True)
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
            previous = len(excluded)
            excluded.update(set(add.index.values.compute().tolist()))
            if previous == len(excluded):
                # Nothing more to remove
                run = False
            assert excluded is not None
            logger.warning("%s Run %s excluding %s", time.ctime(), counter, len(excluded))
        logger.warning("Finished run %s", counter)

    return out, excluded
