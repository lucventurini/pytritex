import pandas as pd
import numpy as np


def find_wrong_assignments(anchored_css: pd.DataFrame, measure: list, sorted_percentile, popseq_percentile,
                           hic_percentile, hic=False):

    """This function will find those scaffolds for which the assignment smells "fishy", ie
    scaffolds where the ratio of second-best to best is over the defined percentiles."""

    melted = pd.melt(anchored_css, id_vars="scaffold_index",
                     value_vars=measure, var_name="map", value_name="chr").dropna()
    #  w[, .N, key=.(scaffold, chr)]->w
    melted = melted.groupby(["scaffold_index", "chr"]).size().to_frame("N").reset_index(drop=False)
    melted = melted.sort_values("N", ascending=False).groupby("scaffold_index").agg(
        {"N": ["sum", "size"]})
    # Nchr_ass: number of assignments for the scaffold.
    # Nchr_ass_uniq: number of unique assignments for the scaffold.
    melted.columns = ["Nchr_ass", "Nchr_ass_uniq"]
    anchored_css = melted.merge(anchored_css, right_on="scaffold_index", left_index=True, how="right")
    anchored_css.loc[:, "scaffold_index"] = anchored_css["scaffold_index"].astype(np.int)
    anchored_css.loc[: "Nchr_ass"] = anchored_css.loc[: "Nchr_ass_uniq"].fillna(0)
    anchored_css.loc[: "Nchr_ass_uniq"] = anchored_css.loc[: "Nchr_ass_uniq"].fillna(0)
    sorted_threshold = anchored_css.loc[anchored_css["Ncss"] >= 30, "sorted_p12"].quantile((sorted_percentile + 1) / 100)
    anchored_css.loc[:, "bad_sorted"] = (anchored_css["sorted_p12"] >= sorted_threshold) & (
            anchored_css["sorted_Ncss2"] >= 2)
    pop_threshold = anchored_css.loc[anchored_css["popseq_Ncss"] >= 30, "popseq_p12"].quantile(
        (popseq_percentile + 1) / 100)
    anchored_css.loc[:, "bad_popseq"] = (anchored_css["popseq_p12"] >= pop_threshold) & (
            anchored_css["popseq_Ncss2"] >= 2)
    anchored_css.loc[anchored_css["bad_sorted"].isna(), "bad_sorted"] = False
    anchored_css.loc[anchored_css["bad_popseq"].isna(), "bad_popseq"] = False
    measure = ["bad_sorted", "bad_popseq"]
    if hic is True:
        hic_threshold = anchored_css.loc[anchored_css["Nhic"] >= 30,
                                         "hic_p12"].quantile((hic_percentile + 1) / 100)
        anchored_css.loc[anchored_css["Nhic"] >= 30, "bad_hic"] = ((anchored_css["hic_p12"] >= hic_threshold) & (
                anchored_css["hic_N2"] >= 2))
        anchored_css.loc[anchored_css["bad_hic"].isna(), "bad_hic"] = False
        measure.append("bad_hic")

    # Finally, let's count how many bad assignments we found for each scaffold. This can be between 0 and 3.
    melted = pd.melt(anchored_css,
                id_vars=["scaffold_index"],
                value_vars=measure,
                value_name="bad",
                var_name="map").dropna().loc[lambda df: df["bad"] == True, :]
    melted = melted.groupby("scaffold_index").size().to_frame("Nbad").reset_index(drop=False)
    anchored_css = melted.merge(anchored_css, on="scaffold_index", how="right")
    anchored_css.loc[:, "Nbad"] = anchored_css["Nbad"].fillna(0)
    return anchored_css
