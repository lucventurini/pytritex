from pytritex.graph_utils.make_super_scaffolds import make_super_scaffolds


def _initial_branch_remover(raw, links, info, excluded, ncores, prefix=None):

    if raw is False:
        run = True
        while run is True:
            out = make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=excluded, ncores=ncores)
            membership = out["membership"]
            res = out["info"]
            a = membership.merge(
                membership.loc[membership["rank"] > 1, ["super", "bin"]].drop_duplicates(),
                how="inner", on=["super", "bin"])
            add = a.loc[a["rank"] == 0, :]
            if add.shape[0] == 0:
                run = False
            else:
                if excluded is None:
                    excluded = add["scaffold_index"].to_list()
                else:
                    excluded = excluded.extend(add["scaffold_index"].to_list())
    else:
        out = make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=excluded, ncores=ncores)
        membership = out["membership"]
        res = out["info"]

    return membership, res, excluded
