import pandas as pd
import numpy as np
from ..utils import first, second


def add_hic_statistics(anchored_css: pd.DataFrame, fpairs: pd.DataFrame, verbose=False):
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
    anchored0 = anchored_css.loc[anchored_css.eval("popseq_chr == sorted_chr"),
                                     ["scaffold_index", "popseq_chr"]].rename(columns={"popseq_chr": "chr"})
    anchored_hic_links = anchored0.copy().rename(
        columns={"chr": "chr1",
                 "scaffold_index": "scaffold_index1"}).merge(fpairs, on="scaffold_index1", how="right")
    anchored_hic_links.loc[:, "scaffold_index1"] = pd.to_numeric(anchored_hic_links["scaffold_index1"],
                                                                 downcast="unsigned")
    anchored_hic_links = anchored0.rename(
        columns={"chr": "chr2",
                 "scaffold_index": "scaffold_index2"}).merge(anchored_hic_links,
                                                             on="scaffold_index2", how="right")
    anchored_hic_links.loc[:, "scaffold_index2"] = pd.to_numeric(anchored_hic_links["scaffold_index2"],
                                                                 downcast="unsigned")
    if verbose:
        print("Anchored HiC links columns", anchored_hic_links.columns)
    hic_stats = anchored_hic_links[~(anchored_hic_links["chr1"].isna() |
                                     anchored_hic_links["chr1"].isna())].rename(
        columns={"scaffold_index2": "scaffold_index", "chr1": "hic_chr"}
    ).groupby(["scaffold_index", "hic_chr"]).size().to_frame("N").reset_index(drop=False).apply(
        pd.to_numeric, downcast="signed")
    N_dtype, chr_dtype = hic_stats["N"].dtype, hic_stats["hic_chr"].dtype
    grouped_hic_stats = hic_stats.sort_values(["scaffold_index", "N"], ascending=[True, False])
    hic_stats = grouped_hic_stats.groupby("scaffold_index", sort=False).agg(
        hic_chr=("hic_chr", first), hic_chr2=("hic_chr", second), hic_N1=("N", first), hic_N2=("N", second),
        Nhic=("N", "sum"))
    type_mapper = {"Nhic": N_dtype, "hic_N1": N_dtype, "hic_N2": N_dtype,
                   "hic_chr": chr_dtype, "hic_chr2": chr_dtype}
    for key in ["Nhic", "hic_N1", "hic_N2", "hic_chr", "hic_chr2"]:
        if hic_stats[key].isna().any() is True and hic_stats[key].isin([0]).any() is True:
            continue
        try:
            hic_stats.loc[:, key] = hic_stats[key].fillna(0).astype(type_mapper[key])
        except ValueError:
            raise ValueError((key, hic_stats[key].value_counts(dropna=False).sort_values().head(20)))
    hic_stats.loc[:, "hic_pchr"] = pd.to_numeric(hic_stats["hic_N1"].div(hic_stats["Nhic"], fill_value=0),
                                                 downcast="float")
    hic_stats.loc[:, "hic_p12"] = pd.to_numeric(hic_stats["hic_N2"].div(hic_stats["hic_N1"], fill_value=0),
                                                downcast="float")
    anchored_css = hic_stats.merge(anchored_css, right_on="scaffold_index", how="right", left_index=True)

    return anchored_css, anchored_hic_links
