import pandas as pd
import numpy as np
from scipy.stats import median_absolute_deviation
from ..utils import first, second


def assign_popseq_position(cssaln: pd.DataFrame, popseq: pd.DataFrame, anchored_css: pd.DataFrame, wheatchr: pd.DataFrame):


    # # Assignment of POPSEQ genetic positions
    #  popseq[!is.na(popseq_alphachr), .(css_contig, popseq_alphachr, popseq_cM)][
    #            cssaln[, .(css_contig, scaffold)], on="css_contig", nomatch=0]->z
    #  z[, .(.N, popseq_cM=mean(popseq_cM), popseq_cM_sd=ifelse(length(popseq_cM) > 1, sd(popseq_cM), 0),
    #  popseq_cM_mad=mad(popseq_cM)), keyby=.(scaffold, popseq_alphachr)]->zz
    #  zz[, popseq_Ncss := sum(N), by=scaffold]->zz
    #  zz[order(-N)][, .(popseq_alphachr=popseq_alphachr[1], popseq_Ncss1=N[1],
    # 			 popseq_alphachr2=popseq_alphachr[2], popseq_Ncss2=N[2]), keyby=scaffold]->x
    #  zz[, .(scaffold, popseq_alphachr, popseq_Ncss, popseq_cM, popseq_cM_sd,  popseq_cM_mad)][x, on=c("scaffold", "popseq_alphachr")]->zz
    #  zz[, popseq_pchr := popseq_Ncss1/popseq_Ncss]
    #  zz[, popseq_p12 := popseq_Ncss2/popseq_Ncss1]
    #  wheatchr[zz, on="popseq_alphachr"]->zz
    #  setnames(copy(wheatchr), paste0(names(wheatchr), 2))[zz, on="popseq_alphachr2"]->zz
    #  zz[info, on="scaffold"]->info
    #  info[is.na(popseq_Ncss), popseq_Ncss := 0]
    #  info[is.na(popseq_Ncss1), popseq_Ncss1 := 0]
    #  info[is.na(popseq_Ncss2), popseq_Ncss2 := 0]

    popseq_positions = popseq.loc[~popseq["popseq_alphachr"].isna(), ["popseq_index", "popseq_alphachr", "popseq_cM"]]
    popseq_positions = popseq_positions.merge(cssaln[["popseq_index", "scaffold_index"]],
                                              on="popseq_index", how="left")
    popseq_stats = popseq_positions.groupby(["scaffold_index", "popseq_alphachr"], observed=True)
    popseq_count = popseq_stats.size().to_frame("N").astype(np.uint32)
    popseq_stats = popseq_stats.agg({"popseq_cM": [np.mean, np.std, median_absolute_deviation]})
    popseq_stats.columns = popseq_stats.columns.to_flat_index()
    popseq_stats.columns = ["popseq_cM", "popseq_cM_sd", "popseq_cM_mad"]
    popseq_stats = pd.merge(popseq_count, popseq_stats, left_index=True, right_index=True).reset_index(drop=False)
    popseq_stats.loc[:, "popseq_Ncss"] = popseq_stats.groupby("scaffold_index")["N"].transform("sum").astype(np.uint32)
    popseq_stats = popseq_stats.reset_index(drop=False)
    # So far we have congregated the different statistics about the centimorgans.
    # Now we have to consider .... ?
    best_locations = popseq_stats.sort_values("N", ascending=False).groupby("scaffold_index").agg(
        popseq_alphachr=("popseq_alphachr", first), popseq_Ncss1=("N", first),
        popseq_alphachr2=("popseq_alphachr", second), popseq_Ncss2=("N", second)).reset_index()
    popseq_stats = popseq_stats[
        ["scaffold_index", "popseq_alphachr", "popseq_Ncss", "popseq_cM", "popseq_cM_sd", "popseq_cM_mad"]].merge(
        best_locations, on=["scaffold_index", "popseq_alphachr"], how="right")
    popseq_stats.loc[:, "popseq_pchr"] = pd.to_numeric(
        popseq_stats["popseq_Ncss1"].div(popseq_stats["popseq_Ncss"], fill_value=0), downcast="float")
    popseq_stats.loc[:, "popseq_p12"] = pd.to_numeric(
        popseq_stats["popseq_Ncss2"].div(popseq_stats["popseq_Ncss1"], fill_value=0), downcast="float")
    wheatchr1 = wheatchr.copy().rename(columns={"chr": "popseq_chr", "alphachr": "popseq_alphachr"})
    popseq_stats = wheatchr1.merge(popseq_stats, on="popseq_alphachr", how="right")
    wheatchr2 = wheatchr.copy().rename(columns={"chr": "popseq_chr2", "alphachr": "popseq_alphachr2"})
    popseq_stats = wheatchr2.merge(popseq_stats, on="popseq_alphachr2", how="right")
    anchored_css = popseq_stats.merge(anchored_css, on="scaffold_index", how="right")
    for column in ["popseq_Ncss", "popseq_Ncss1", "popseq_Ncss2"]:
        anchored_css.loc[:, column] = pd.to_numeric(anchored_css[column].fillna(0),
                                                    downcast="unsigned")
    return anchored_css
