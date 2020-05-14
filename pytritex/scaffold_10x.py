import pandas as pd
import numpy as np
from .scaffold_10x_utils.branch_remover import _initial_branch_remover
from .scaffold_10x_utils.link_finder import _initial_link_finder
from .scaffold_10x_utils.orient_scaffolds import orient_scaffolds
from .scaffold_10x_utils.non_raw_analyser import non_raw_analyser


def scaffold_10x(assembly: dict, prefix="super", min_npairs=5, max_dist=1e5, min_nmol=6,
                min_nsample=2, popseq_dist=5, max_dist_orientation=5,
                ncores=1, verbose=True, raw=False, unanchored=True):
    """Construct super-scaffold from 10X links. Use heuristics based on genetic map information
    to prune erroneous edges in the scaffold graph."""

    info, molecules = assembly["info"], assembly["molecules"]
    print(assembly.keys())
    ww2, links, link_pos = _initial_link_finder(info=info, molecules=molecules,
                                      verbose=verbose, popseq_dist=popseq_dist,
                                      min_npairs=min_npairs, max_dist=max_dist, min_nmol=min_nmol,
                                      min_nsample=min_nsample)
    excluded = []
    membership, res, excluded = _initial_branch_remover(
        raw, links=links, info=info, excluded=excluded, ncores=ncores, prefix=prefix)

    if verbose is True:
        print("Finding initial super-scaffolds")

    if raw is False:
        membership, res, excluded = non_raw_analyser(links=links, excluded=excluded,
                         membership=membership, info=info,
                         prefix=prefix, ncores=ncores, verbose=verbose, unanchored=unanchored,
                         popseq_dist=popseq_dist)

    result = orient_scaffolds(info=info, res=res, membership=membership,
                              link_pos=link_pos,
                              max_dist_orientation=max_dist_orientation,
                              verbose=verbose)
    return result
