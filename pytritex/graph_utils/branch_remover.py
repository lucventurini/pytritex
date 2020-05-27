from pytritex.graph_utils.make_super_scaffolds import make_super_scaffolds
import time


def _initial_branch_remover(raw, links, info, excluded, ncores, prefix=None):

    print(time.ctime(), "Starting the run")
    if raw is False:
        run = True
        counter = 0
        while run is True:
            counter += 1
            print(time.ctime(), "Starting run", counter)
            out = make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=excluded, ncores=ncores)
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
        out = make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=excluded, ncores=ncores)
        # membership = out["membership"]
        # res = out["info"]

    return out, excluded
