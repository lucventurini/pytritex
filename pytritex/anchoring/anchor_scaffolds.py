from pytritex.utils.chrnames import chrNames
from .assign_carma import assign_carma
from .assign_popseq_position import assign_popseq_position
from .add_hic_statistics import add_hic_statistics
from .find_wrong_assignment import find_wrong_assignments
import dask.dataframe as dd
from dask.distributed import Client
from dask import delayed
import os


def anchor_scaffolds(assembly: dict,
                     save: str,
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

    wheatchr = chrNames(species=species)
    fai = dd.read_parquet(assembly["fai"])
    cssaln = dd.read_parquet(assembly["cssaln"])
    assert isinstance(cssaln, dd.DataFrame)
    popseq = dd.read_parquet(assembly["popseq"])
    assert isinstance(popseq, dd.DataFrame)
    if "fpairs" not in assembly:
        fpairs = None
        hic = False
    else:
        fpairs = dd.read_parquet(assembly["fpairs"])
        hic = (fpairs.head(5).shape[0] > 0)

    anchored_css = assign_carma(cssaln, fai, wheatchr)
    anchored_css = assign_popseq_position(
        cssaln=cssaln, popseq=popseq,
        anchored_css=anchored_css, wheatchr=wheatchr)
    if hic is True:
        anchored_css, anchored_hic_links = add_hic_statistics(anchored_css, fpairs)
        measure = ["popseq_chr", "hic_chr", "sorted_chr"]
    else:
        measure = ["popseq_chr", "sorted_chr"]
        anchored_hic_links = None
    anchored_css = find_wrong_assignments(anchored_css, measure,
        sorted_percentile=sorted_percentile, hic_percentile=hic_percentile,
        popseq_percentile=popseq_percentile, hic=hic)
    assert isinstance(anchored_css, dd.DataFrame), (type(anchored_css))
    # save_dir = os.path.join(save, "joblib", "pytritex", "anchoring")
    assert isinstance(anchored_css, dd.DataFrame)
    dd.to_parquet(anchored_css, os.path.join(save, "anchored_css"), compute=True)
    assembly["info"] = os.path.join(save, "anchored_css")
    if hic is True:
        dd.to_parquet(anchored_hic_links, os.path.join(save, "anchored_hic_links"), compute=True)
        assembly["fpairs"] = os.path.join(save, "anchored_hic_links")
    return assembly
