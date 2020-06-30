from ..graph_utils.branch_remover import _initial_branch_remover
from .link_finder import _initial_link_finder
from .orient_scaffolds import orient_scaffolds
from pytritex.graph_utils.tip_removal.tip_remover import remove_tips
from .scaffold_unanchored import _scaffold_unanchored
import time
from joblib import Memory
import os
import hashlib
from dask.distributed import Client

sha = hashlib.sha256()


def scaffold_10x(assembly: dict, memory: Memory, save_dir: str,
                 client: Client,
                 min_npairs=5, max_dist=1e5, min_nmol=6,
                 min_nsample=2, popseq_dist=5, max_dist_orientation=5,
                 ncores=1, verbose=True, raw=False, unanchored=True):
    """Construct super-scaffold from 10X links. Use heuristics based on genetic map information
    to prune erroneous edges in the scaffold graph."""

    info, molecules, fai = assembly["info"], assembly["molecules"], assembly["fai"]
    # print(assembly.keys())
    print(time.ctime(), "Finding initial links")

    params = (min_npairs, max_dist, min_nmol, min_nsample, popseq_dist, max_dist_orientation, raw)
    sha.update(str(params).encode())
    folder = os.path.join(save_dir, "joblib", "pytritex", "scaffold_10x", sha.hexdigest())
    sample_count, links, link_pos = memory.cache(_initial_link_finder,
                                                 ignore=["client"])(
        info=info, molecules=molecules,
        save_dir=folder,
        client=client,
        fai=fai, verbose=verbose, popseq_dist=popseq_dist,
        min_npairs=min_npairs, max_dist=max_dist, min_nmol=min_nmol,
        min_nsample=min_nsample)
    print(time.ctime(), "Found initial links")
    excluded = set()
    print(time.ctime(), "Starting initial pruning")
    out, excluded = memory.cache(_initial_branch_remover, ignore=["ncores", "client"])(
        client=client, save_dir=folder, links=links, info=info, excluded=excluded, ncores=ncores)
    res = out["info"]
    membership = out["membership"]
    print(time.ctime(), "Finished initial pruning")

    if verbose is True:
        print("Finding initial super-scaffolds")

    if raw is False:

        membership, res, excluded = memory.cache(remove_tips, ignore=["client", "ncores", "verbose"])(
            links=links, excluded=excluded,
            out=out, info=info,
            client=client,
            save_dir=os.path.join(folder, "tip_removal"),
            ncores=ncores, verbose=verbose,
            min_dist=1e4)
        if popseq_dist > 0 and unanchored is True:

            membership, res = memory.cache(_scaffold_unanchored,
                                           ignore=["client", "ncores", "verbose"])(
                links,
                excluded, membership, info,
                sample_count,
                save_dir=os.path.join(folder, "unanchored"), client=client,
                ncores=1, verbose=False)

    membership, result = memory.cache(orient_scaffolds, ignore=["client"])(
        info=info, res=res, membership=membership,
        link_pos=link_pos, max_dist_orientation=max_dist_orientation,
        save_dir=os.path.join(folder, "orientation"),
        client=client)

    return membership, result
