import pandas as pd
from pytritex.utils.chrnames import chrNames
from .anchor_css_alignment import anchor_css_alignment
from .add_centimorgans import add_centimorgans
from .add_hic_statistics import add_hic_statistics
from .find_wrong_assignment import find_wrong_assignments
import numpy as np


def anchor_scaffolds(assembly: dict,
                     popseq: pd.DataFrame,
                     species=None,
                     sorted_percentile=95,
                     popseq_percentile=90,
                     hic_percentile=98):
    if species is None:
        raise KeyError(
            "Parameter 'species' is NULL. Please set 'species' to one of "
            "\"wheat\", \"barley\", \"oats\", \"lolium\", \"sharonensis\" or \"rye\".")
    elif species not in ("wheat", "barley", "rye", "oats", "sharonensis", "lolium"):
        raise KeyError(
            "Parameter 'species' is not valid. Please set 'species' to one of "
            "\"wheat\", \"barley\", \"oats\", \"lolium\", \"sharonensis\" or \"rye\".")
    wheatchr = chrNames(species=species).rename(columns={"alphachr": "popseq_alphachr", "chr": "popseq_chr"})

    fai = assembly["fai"]
    initial_size = fai.shape[0]
    cssaln = assembly["cssaln"]
    if "fpairs" not in assembly:
        fpairs = None
        hic = False
    else:
        fpairs = assembly["fpairs"]
        hic = (assembly["fpairs"].shape[0] > 0)

    anchored_css = anchor_css_alignment(cssaln, fai, wheatchr)
    anchored_css.loc[:, "scaffold_index"] = anchored_css["scaffold_index"].astype(np.int)
    # anchored_css.loc[:, "popseq_index"] = anchored_css["popseq_index"].astype(np.int)
    anchored_css = add_centimorgans(cssaln, popseq, anchored_css, wheatchr)
    anchored_css.loc[:, "scaffold_index"] = anchored_css["scaffold_index"].astype(np.int)
    # anchored_css.loc[:, "popseq_index"] = anchored_css["popseq_index"].astype(np.int)

    # # Assignment of POPSEQ genetic positions
    if hic is True:
        anchored_css, anchored_hic_links = add_hic_statistics(anchored_css, fpairs)
        measure = ["popseq_chr", "hic_chr", "sorted_chr"]
    else:
        measure = ["popseq_chr", "sorted_chr"]
        anchored_hic_links = None
    anchored_css = find_wrong_assignments(anchored_css, measure,
                                          sorted_percentile=sorted_percentile, hic_percentile=hic_percentile,
                                          popseq_percentile=popseq_percentile, hic=hic)

    anchored_css.loc[:, "scaffold_index"] = anchored_css["scaffold_index"].fillna(0).astype(np.int)
    assembly["info"] = anchored_css
    assembly["popseq"] = popseq
    if hic is True:
        assembly["fpairs"] = anchored_hic_links
    print(anchored_css.columns)
    return assembly
