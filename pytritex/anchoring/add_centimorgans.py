import pandas as pd
import numpy as np
from scipy.stats import median_absolute_deviation
from ..utils import first, second


def add_centimorgans(cssaln: pd.DataFrame, popseq: pd.DataFrame, anchored_css: pd.DataFrame, wheatchr: pd.DataFrame):
    centimorgans = pd.merge(
        popseq.loc[~popseq["sorted_alphachr"].isna(), ["popseq_index", "popseq_alphachr", "popseq_cM"]],
        cssaln.loc[:, ["popseq_index", "scaffold_index"]],
        on="popseq_index", how="left")
    centi_grouped = centimorgans.groupby(["scaffold_index", "popseq_alphachr"], observed=True)
    centi_count = centi_grouped.size().to_frame("N")
    centi_stats = centi_grouped.agg({"popseq_cM": [np.mean, np.std, median_absolute_deviation]})
    centi_stats.columns = centi_stats.columns.to_flat_index()
    centi_stats = centi_stats.rename(
        columns={centi_stats.columns[0]: "popseq_cM",
                 centi_stats.columns[1]: "popseq_cM_sd",
                 centi_stats.columns[2]: "popseq_cM_mad"})
    centi_stats = pd.merge(centi_count, centi_stats, left_index=True, right_index=True).reset_index(drop=False)
    centi_stats.loc[:, "popseq_Ncss"] = centi_stats.groupby("scaffold_index")["N"].transform("size")
    # So far we have congregated the different statistics about the centimorgans.
    # Now we have to consider .... ?
    best_locations = centi_stats.sort_values(
        ["scaffold_index", "N"], ascending=[True, False])[
        ["scaffold_index", "N", "popseq_alphachr"]].drop_duplicates().groupby("scaffold_index").agg(
        {"N": [first, second], "popseq_alphachr": [first, second]
         })
    best_locations.columns = ["popseq_Ncss1", "popseq_Ncss2", "popseq_alphachr", "popseq_alphachr2"]
    best_locations.reset_index(drop=False, inplace=True)
    centi_stats = centi_stats.loc[:,
                  ["scaffold_index", "popseq_alphachr", "popseq_Ncss",
                    "popseq_cM", "popseq_cM_sd", "popseq_cM_mad"]].merge(
        best_locations, on=["scaffold_index", "popseq_alphachr"], how="right")
    centi_stats.loc[:, "popseq_pchr"] = centi_stats["popseq_Ncss1"].div(centi_stats["popseq_Ncss"], fill_value=0)
    centi_stats.loc[:, "popseq_p12"] = centi_stats["popseq_Ncss2"].div(centi_stats["popseq_Ncss1"], fill_value=0)
    centi_stats = wheatchr.rename(columns={"popseq_chr": "popseq_chr2",
                                  "popseq_alphachr": "popseq_alphachr2"}).merge(
        wheatchr.merge(centi_stats, on="popseq_alphachr"), on="popseq_alphachr2")
    anchored_css = centi_stats.copy().merge(anchored_css.copy(), on="scaffold_index", how="right")
    for column in ["popseq_Ncss", "popseq_Ncss1", "popseq_Ncss2"]:
        anchored_css.loc[:, column] = anchored_css[column].fillna(0)
    return anchored_css
