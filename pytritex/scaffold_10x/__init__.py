from ..graph_utils.branch_remover import _initial_branch_remover
from .link_finder import _initial_link_finder
from .orient_scaffolds import orient_scaffolds
from .tip_remover import remove_tips
from .scaffold_unanchored import _scaffold_unanchored
import time


def scaffold_10x(assembly: dict, prefix="super", min_npairs=5, max_dist=1e5, min_nmol=6,
                 min_nsample=2, popseq_dist=5, max_dist_orientation=5,
                 ncores=1, verbose=True, raw=False, unanchored=True):
    """Construct super-scaffold from 10X links. Use heuristics based on genetic map information
    to prune erroneous edges in the scaffold graph."""

    info, molecules = assembly["info"], assembly["molecules"]
    # print(assembly.keys())
    print(time.ctime(), "Finding initial links")
    sample_count, links, link_pos = _initial_link_finder(info=info, molecules=molecules,
                                      verbose=verbose, popseq_dist=popseq_dist,
                                      min_npairs=min_npairs, max_dist=max_dist, min_nmol=min_nmol,
                                      min_nsample=min_nsample)
    print(time.ctime(), "Found initial links")
    excluded = set()
    print(time.ctime(), "Starting initial pruning")
    out, excluded = _initial_branch_remover(
        raw, links=links, info=info, excluded=excluded, ncores=ncores, prefix=prefix)
    res = out["info"]
    membership=out["membership"]
    print(time.ctime(), "Finished initial pruning")

    if verbose is True:
        print("Finding initial super-scaffolds")

    if raw is False:
        membership, res, excluded = remove_tips(links=links, excluded=excluded,
                                                out=out, info=info,
                                                prefix=prefix, ncores=ncores, verbose=verbose,
                                                min_dist=1e4)
        if popseq_dist > 0 and unanchored is True:
            membership, res = _scaffold_unanchored(links, excluded, membership, info, sample_count, ncores=1,
                                                   prefix=prefix, verbose=False)
    membership, result = orient_scaffolds(
        info=info, res=res, membership=membership,
        link_pos=link_pos, max_dist_orientation=max_dist_orientation, verbose=verbose)
    return membership, result
