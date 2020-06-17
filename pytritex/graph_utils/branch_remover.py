from pytritex.graph_utils.make_super_scaffolds import make_super_scaffolds
import time
import dask.dataframe as dd
from joblib import Memory


def _initial_branch_remover(raw, memory: Memory, save_dir: str, links: str,
                            info: str, excluded, ncores):

    print(time.ctime(), "Starting the run")
    if raw is False:
        run = True
        counter = 0
        while run is True:
            counter += 1
            print(time.ctime(), "Starting run", counter)
            out = memory.cache(make_super_scaffolds, ignore=["ncores"])(
                links=links, save_dir=save_dir, info=info, excluded=excluded, ncores=ncores)
            membership = out["membership"]
            # res = out["info"]
            a = membership.merge(
                membership.loc[membership["rank"] > 1, ["super", "bin"]].drop_duplicates(),
                how="inner", on=["super", "bin"])
            add = a.loc[a["rank"] == 0, :]
            if add.shape[0] == 0:
                run = False
            else:
                if excluded is None:
                    excluded = set(add["scaffold_index"].to_list())
                else:
                    previous = len(excluded)
                    excluded.update(set(add["scaffold_index"].to_list()))
                    if previous == len(excluded):
                        # Nothing more to remove
                        run = False
                assert excluded is not None
                print(time.ctime(), "Run", counter, ", excluding", len(excluded))
    else:
        out = memory.cache(make_super_scaffolds, ignore=["ncores"])(
                links=links, save_dir=save_dir, info=info, excluded=excluded, ncores=ncores)

    return out, excluded
