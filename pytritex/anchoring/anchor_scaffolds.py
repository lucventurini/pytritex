from pytritex.utils.chrnames import chrNames
from .assign_carma import assign_carma
from .assign_popseq_position import assign_popseq_position
from .add_hic_statistics import add_hic_statistics
from .find_wrong_assignment import find_wrong_assignments
import dask.dataframe as dd
from dask.distributed import Client
from dask.delayed import delayed
from typing import Union
import os
import numpy as np
np.seterr(all='raise')


def anchor_scaffolds(assembly: dict,
                     save: Union[str, None],
                     client: Client,
                     scaffolds=None,
                     species=None,
                     sorted_percentile=95,
                     popseq_percentile=90,
                     hic_percentile=98) -> dict:
    if species is None:
        raise KeyError(
            "Parameter 'species' is NULL. Please set 'species' to one of "
            "\"wheat\", \"barley\", \"oats\", \"lolium\", \"sharonensis\" or \"rye\".")
    elif species not in ("wheat", "barley", "rye", "oats", "sharonensis", "lolium"):
        raise KeyError(
            "Parameter 'species' is not valid. Please set 'species' to one of "
            "\"wheat\", \"barley\", \"oats\", \"lolium\", \"sharonensis\" or \"rye\".")

    wheatchr = chrNames(species=species)
    if isinstance(assembly["fai"], str):
        fai = dd.read_parquet(assembly["fai"], infer_divisions=True)
    else:
        fai = assembly["fai"]

    if isinstance(assembly["cssaln"], str):
        cssaln = dd.read_parquet(assembly["cssaln"], infer_divisions=True)
    else:
        cssaln = assembly["cssaln"]
    assert isinstance(cssaln, dd.DataFrame)
    popseq = dd.read_parquet(assembly["popseq"])
    assert isinstance(popseq, dd.DataFrame)
    if "fpairs" not in assembly:
        fpairs = None
        hic = False
    else:
        fpairs = assembly["fpairs"]
        if isinstance(fpairs, str):
            fpairs = dd.read_parquet(fpairs, infer_divisions=True)
        hic = (fpairs.head(5, npartitions=-1).shape[0] > 0)

    # if scaffolds is not None:
    #     anchored_css = assembly["anchored_css"]
    #
    # else:
    anchored_css = assign_carma(cssaln, fai, wheatchr, client)
    anchored_css = assign_popseq_position(
        cssaln=cssaln, popseq=popseq, client=client,
        anchored_css=anchored_css, wheatchr=wheatchr)
    if hic is True:
        anchored_css, anchored_hic_links = add_hic_statistics(anchored_css, fpairs, client=client)
        measure = ["popseq_chr", "hic_chr", "sorted_chr"]
    else:
        measure = ["popseq_chr", "sorted_chr"]
        anchored_hic_links = None
    anchored_css = find_wrong_assignments(
        anchored_css, measure, client=client,
        sorted_percentile=sorted_percentile, hic_percentile=hic_percentile,
        popseq_percentile=popseq_percentile, hic=hic)
    assert isinstance(anchored_css, dd.DataFrame), (type(anchored_css))
    assert isinstance(anchored_css, dd.DataFrame)
    if save is not None:
        dd.to_parquet(anchored_css, os.path.join(save, "anchored_css"), compute=True)
        assembly["info"] = os.path.join(save, "anchored_css")
    else:
        assembly["info"] = anchored_css
    if hic is True:
        if save is not None:
            dd.to_parquet(anchored_hic_links, os.path.join(save, "anchored_hic_links"), compute=True)
            assembly["fpairs"] = os.path.join(save, "anchored_hic_links")
        else:
            assembly["fpairs"] = anchored_hic_links

    return assembly
