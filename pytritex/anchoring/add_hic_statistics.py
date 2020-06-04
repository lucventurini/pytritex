import pandas as pd
import numpy as np
from ..utils import first, second
import dask.dataframe as dd


def add_hic_statistics(anchored_css: dd.DataFrame, fpairs: dd.DataFrame, verbose=False):
    """This function will add the HiC statistics to the anchored dataframe."""

    # info[!(popseq_chr != sorted_chr)][, .(scaffold, chr=popseq_chr)]->info0
    #   setnames(copy(info0), names(info0), sub("$", "1", names(info0)))[fpairs, on="scaffold1"]->tcc_pos
    #   setnames(copy(info0), names(info0), sub("$", "2", names(info0)))[tcc_pos, on="scaffold2"]->tcc_pos
    #   tcc_pos[!is.na(chr1), .N, key=.(scaffold=scaffold2, hic_chr=chr1)]->z
    #   z[order(-N)][, .(Nhic=sum(N), hic_chr=hic_chr[1], hic_N1=N[1],
    # 		   hic_chr2=hic_chr[2], hic_N2=N[2]), keyby=scaffold]->zz
    #   zz[, hic_pchr := hic_N1/Nhic]
    #   zz[, hic_p12 := hic_N2/hic_N1]
    #   zz[info, on="scaffold"]->info
    #   info[is.na(Nhic), Nhic := 0]
    #   info[is.na(hic_N1), hic_N1 := 0]
    #   info[is.na(hic_N2), hic_N2 := 0]
    anchored0 = anchored_css.query("popseq_chr == sorted_chr")[["popseq_chr"]].compute().reset_index(drop=False)
    anchored0 = anchored0.rename(columns={"popseq_chr": "chr"}, errors="raise")
    anchored_hic_links = dd.merge(anchored0.copy().rename(
        columns={"chr": "chr1",
                 "scaffold_index": "scaffold_index1"},
    errors="raise"), fpairs, on="scaffold_index1", how="right")
    anchored_hic_links = dd.merge(anchored0.rename(
        columns={"chr": "chr2",
                 "scaffold_index": "scaffold_index2"}),
        anchored_hic_links, on="scaffold_index2", how="right")
    if verbose:
        print("Anchored HiC links columns", anchored_hic_links.columns)
    hic_stats = anchored_hic_links.query("chr1 == chr1").rename(
        columns={"scaffold_index2": "scaffold_index", "chr1": "hic_chr"})
    hic_stats = hic_stats.groupby(
        ["scaffold_index", "hic_chr"]).size().to_frame("N").reset_index(drop=False).compute()
    hic_stats["scaffold_index"] = pd.to_numeric(hic_stats["scaffold_index"],
                                                downcast="integer")
    grouped_hic_stats = hic_stats.sort_values(["scaffold_index", "N"], ascending=[True, False])
    hic_stats = grouped_hic_stats.groupby(
        "scaffold_index", sort=False).agg(
        hic_chr=("hic_chr", first), hic_chr2=("hic_chr", second), hic_N1=("N", first), hic_N2=("N", second),
        Nhic=("N", "sum"))
    hic_stats.loc[:, "hic_pchr"] = pd.to_numeric(hic_stats["hic_N1"].div(hic_stats["Nhic"], fill_value=0),
                                                 downcast="float")
    hic_stats.loc[:, "hic_p12"] = pd.to_numeric(hic_stats["hic_N2"].div(hic_stats["hic_N1"], fill_value=0),
                                             downcast="float")
    hic_stats = dd.from_pandas(hic_stats, npartitions=100)
    anchored_css = dd.merge(hic_stats, anchored_css, on="scaffold_index", how="right")
    anchored_css = anchored_css.persist()
    return anchored_css, anchored_hic_links
