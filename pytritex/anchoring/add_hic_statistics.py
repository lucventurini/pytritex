import pandas as pd
import numpy as np
from ..utils import first, second


def add_hic_statistics(anchored_css: pd.DataFrame, fpairs: pd.DataFrame):
    """This function will add the HiC statistics to the anchored dataframe."""

    anchored0 = anchored_css[anchored_css["popseq_chr"] == anchored_css["sorted_chr"]][
        ["scaffold_index", "popseq_chr"]].rename(columns={"popseq_chr": "chr"})
    anchored_hic_links = anchored0.copy().rename(
        columns={"chr": "chr1", "scaffold_index": "scaffold_index1"}).merge(
        fpairs, on="scaffold_index1", how="right")
    anchored_hic_links = anchored0.rename(
        columns={"chr": "chr2", "scaffold_index": "scaffold_index2"}).merge(
        anchored_hic_links, on="scaffold_index2", how="right")
    hic_stats = anchored_hic_links[~anchored_hic_links["chr1"].isna()].rename(
        columns={"scaffold_index2": "scaffold_index", "chr1": "hic_chr"}
    ).groupby(["scaffold_index", "hic_chr"]).size().to_frame("N").reset_index(drop=False)
    grouped_hic_stats = hic_stats.sort_values(["scaffold_index", "N"], ascending=[True, False])
    Nhic = grouped_hic_stats.groupby("scaffold_index", sort=False, observed=True).size().to_frame("Nhic")
    hic_stats = grouped_hic_stats[
        ["scaffold_index", "N", "hic_chr"]].drop_duplicates().groupby("scaffold_index", sort=False).agg(
        hic_chr=("hic_chr", first), hic_chr2=("hic_chr", second), hic_N1=("N", first), hic_N2=("N", second)).merge(
        Nhic, left_index=True, right_index=True)
    hic_stats["hic_pchr"] = hic_stats["hic_N1"].div(hic_stats["Nhic"], fill_value=0)
    hic_stats["hic_p12"] = hic_stats["hic_N2"].div(hic_stats["hic_N1"], fill_value=0)
    anchored_css = hic_stats.merge(anchored_css, on="scaffold_index", how="right", left_index=True)
    for col in ["Nhic", "hic_N1", "hic_N2"]:
        anchored_css.loc[:, col] = anchored_css[col].fillna(0)
    return anchored_css, anchored_hic_links
