from ..graph_utils.branch_remover import _initial_branch_remover
from .link_finder import _initial_link_finder
from .orient_scaffolds import orient_scaffolds
from .non_raw_analyser import non_raw_analyser
import time


def scaffold_10x(assembly: dict, prefix="super", min_npairs=5, max_dist=1e5, min_nmol=6,
                min_nsample=2, popseq_dist=5, max_dist_orientation=5,
                ncores=1, verbose=True, raw=False, unanchored=True):
    """Construct super-scaffold from 10X links. Use heuristics based on genetic map information
    to prune erroneous edges in the scaffold graph."""

    info, molecules = assembly["info"], assembly["molecules"]
    # print(assembly.keys())
    print(time.ctime(), "Finding initial links")
    ww2, links, link_pos = _initial_link_finder(info=info, molecules=molecules,
                                      verbose=verbose, popseq_dist=popseq_dist,
                                      min_npairs=min_npairs, max_dist=max_dist, min_nmol=min_nmol,
                                      min_nsample=min_nsample)
    print(time.ctime(), "Found initial links")
    excluded = set()
    print(time.ctime(), "Starting initial pruning")
    membership, res, excluded = _initial_branch_remover(
        raw, links=links, info=info, excluded=excluded, ncores=ncores, prefix=prefix)
    print(time.ctime(), "Finished initial pruning")

    if verbose is True:
        print("Finding initial super-scaffolds")

    if raw is False:
        membership, res, excluded = non_raw_analyser(links=links, excluded=excluded,
                         membership=membership, info=info,
                         prefix=prefix, ncores=ncores, verbose=verbose, unanchored=unanchored,
                         popseq_dist=popseq_dist)

    membership, result = orient_scaffolds(info=info, res=res, membership=membership,
                              link_pos=link_pos,
                              max_dist_orientation=max_dist_orientation,
                              verbose=verbose)
    return membership, result
