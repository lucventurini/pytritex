from pytritex.utils.chrnames import chrNames
from .assign_carma import assign_carma
from .assign_popseq_position import assign_popseq_position
from .add_hic_statistics import add_hic_statistics
from .find_wrong_assignment import find_wrong_assignments
import dask.dataframe as dd
from dask.distributed import Client
import os


def anchor_scaffolds(assembly: dict,
                     save: str,
                     client: Client,
                     species=None,
                     sorted_percentile=95,
                     popseq_percentile=90,
                     hic_percentile=98,
                     verbose=False):
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
    popseq = dd.read_parquet(assembly["popseq"])
    if "fpairs" not in assembly:
        fpairs = None
        hic = False
    else:
        fpairs = dd.read_parquet(assembly["fpairs"])
        hic = (fpairs.head(5).shape[0] > 0)

    scattered_css = client.scatter(cssaln)
    scattered_fai = client.scatter(fai)

    anchored_css = client.submit(assign_carma, scattered_css, scattered_fai, wheatchr)
    anchored_css = client.submit(assign_popseq_position, cssaln=scattered_css,
                                 popseq=popseq, anchored_css=anchored_css,
                                 wheatchr=wheatchr)
    anchored_css = client.gather(anchored_css)
    anchored_css = client.scatter(anchored_css)
    if hic is True:
        assigner = client.submit(add_hic_statistics, anchored_css, fpairs)
        anchored_css, anchored_hic_links = client.gather(assigner)
        anchored_css = client.scatter(anchored_css)
        measure = ["popseq_chr", "hic_chr", "sorted_chr"]
    else:
        measure = ["popseq_chr", "sorted_chr"]
        anchored_hic_links = None
    anchored_css = client.submit(find_wrong_assignments,
                                 anchored_css, measure,
                                 sorted_percentile=sorted_percentile, hic_percentile=hic_percentile,
                                 popseq_percentile=popseq_percentile, hic=hic)
    # anchored_css = client.gather(anchored_css)
    # anchored_css.loc[:, "scaffold_index"] = anchored_css["scaffold_index"].fillna(0).astype(np.int)
    # Store in parquet
    anchored_css = client.gather(anchored_css)
    save_dir = os.path.join(save, "joblib", "pytritex", "anchoring")
    parqueting = client.submit(dd.to_parquet, anchored_css, os.path.join(save_dir, "anchored_css"),
                               compute=True)
    client.gather(parqueting)
    assembly["info"] = os.path.join(save_dir, "anchored_css")
    if hic is True:
        parqueting = client.submit(dd.to_parquet,
                                   anchored_hic_links, os.path.join(save_dir, "anchored_hic_links"))
        client.gather(parqueting)
        assembly["fpairs"] = os.path.join(save_dir, "anchored_hic_links")
    return assembly
