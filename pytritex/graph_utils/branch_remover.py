from pytritex.graph_utils.make_super_scaffolds import make_super_scaffolds
import time
from joblib import Memory
from dask.distributed import Client
import dask.dataframe as dd
import logging
logger = logging.getLogger("distributed.worker")


def _initial_branch_remover(raw, memory: Memory, client: Client,
                            save_dir: str,
                            links: str,
                            info: str, excluded: set, ncores):

    print(time.ctime(), "Starting the run")
    if excluded is None:
        excluded = set()

    if raw is False:
        run = True
        counter = 0
        while run is True:
            counter += 1
            # print(time.ctime(), "Starting run", counter)
            logger.info("Starting run %s , excluded: %s",
                        counter, len(excluded))
            out = memory.cache(make_super_scaffolds, ignore=["ncores", "client"])(
                links=links, save_dir=save_dir,
                client=client,
                info=info, excluded=excluded, ncores=ncores)
            membership = dd.read_parquet(out["membership"], infer_divisions=True)
            # res = out["info"]
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
                print(time.ctime(), "Run", counter, ", excluding", len(excluded))
            logger.info("Finished run %s", counter)
    else:
        out = memory.cache(make_super_scaffolds, ignore=["ncores"])(
                links=links, client=client,
            save_dir=save_dir, info=info, excluded=excluded, ncores=ncores)

    return out, excluded
