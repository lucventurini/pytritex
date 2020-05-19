import pandas as pd
from ..utils import first, second
import numpy as np


def assign_carma(cssaln: pd.DataFrame, fai: pd.DataFrame, wheatchr: pd.DataFrame):
    """
    This function will create a master table with the anchored information from the CSS alignment.
    :param cssaln: the original CSS alignment.
    :param fai: the original information on the scaffolds.
    :param species: the species we are working on (e.g. wheat).
    :return:
    """

    # # Assignment of CARMA chromosome
    #  cssaln[!is.na(sorted_alphachr), .N, keyby=.(scaffold, sorted_alphachr)]->z
    #  z[order(-N)][, .(Ncss=sum(N), sorted_alphachr=sorted_alphachr[1], sorted_Ncss1=N[1],
    # 			 sorted_alphachr2=sorted_alphachr[2], sorted_Ncss2=N[2]), keyby=scaffold]->z
    #  z[, sorted_pchr := sorted_Ncss1/Ncss]
    #  z[, sorted_p12 := sorted_Ncss2/sorted_Ncss1]
    anchored_css = cssaln[~cssaln["sorted_alphachr"].isna()].copy()
    anchored_css_grouped = anchored_css.groupby(["scaffold_index", "sorted_alphachr"], observed=True)

    combined_stats = anchored_css_grouped.size().to_frame("N").reset_index(drop=False)
    combined_stats = combined_stats.sort_values(
        ["scaffold_index", "N"], ascending=[True, False]).groupby("scaffold_index").agg(
        Ncss=("N", "sum"), sorted_Ncss1=("N", first), sorted_Ncss2=("N", second),
        sorted_alphachr=("sorted_alphachr", first),
        sorted_alphachr2=("sorted_alphachr", second)
    )
    combined_stats.loc[:, "sorted_pchr"] = combined_stats["sorted_Ncss1"].div(combined_stats["Ncss"], fill_value=0)
    combined_stats.loc[:, "sorted_p12"] = combined_stats["sorted_Ncss2"].div(
        combined_stats["sorted_Ncss1"], fill_value=0)

    #  # Assignment of CARMA chromosome arm
    #  cssaln[sorted_arm == "L", .(NL=.N), keyby=.(scaffold, sorted_alphachr)]->al
    #  cssaln[sorted_arm == "S", .(NS=.N), keyby=.(scaffold, sorted_alphachr)]->as
    #  al[z, on=c("scaffold", "sorted_alphachr")]->z
    #  as[z, on=c("scaffold", "sorted_alphachr")]->z
    #  z[is.na(NL), NL := 0]
    #  z[is.na(NS), NS := 0]
    #  z[, sorted_arm := ifelse(NS >=  NL, "S", "L")]
    #  z[sorted_alphachr %in% c("1H", "3B"), sorted_arm := NA]
    #  z[sorted_arm == "S", sorted_parm := NS/sorted_Ncss1]
    #  z[sorted_arm == "L", sorted_parm := NL/sorted_Ncss1]
    short_arm_counts = cssaln.loc[cssaln.sorted_arm == "S"].groupby(
        ["scaffold_index", "sorted_alphachr"])["sorted_arm"].transform("size").to_frame("NS")
    long_arm_counts = cssaln.loc[cssaln.sorted_arm == "L"].groupby(
        ["scaffold_index", "sorted_alphachr"])["sorted_arm"].transform("size").to_frame("NL")

    combined_stats = short_arm_counts.merge(
        long_arm_counts.merge(combined_stats, left_index=True, right_index=True, how="right"),
        left_index=True, right_index=True, how="right")
    combined_stats.loc[:, ["NL", "NS"]] = combined_stats[["NL", "NS"]].fillna(0)
    combined_stats.loc[:, "sorted_arm"] = (combined_stats["NS"] > combined_stats["NL"]).map(
        {True: "S", False: "L"}).astype(bool)
    combined_stats.loc[combined_stats["sorted_alphachr"].isin(["1H", "3B"]), "sorted_arm"] = np.nan
    combined_stats.loc[combined_stats["sorted_arm"] == "S", "sorted_parm"] = combined_stats["NS"].div(
        combined_stats["sorted_Ncss1"], fill_value=0)
    combined_stats.loc[combined_stats["sorted_arm"] == "L", "sorted_parm"] = combined_stats["NL"].div(
        combined_stats["sorted_Ncss1"], fill_value=0)
    combined_stats = combined_stats.reset_index(drop=False)

    #  setnames(copy(wheatchr), c("sorted_alphachr", "sorted_chr"))[z, on="sorted_alphachr"]->z
    #  setnames(copy(wheatchr), c("sorted_alphachr2", "sorted_chr2"))[z, on="sorted_alphachr2"]->z
    #  z[fai, on="scaffold"]->info
    #  info[is.na(Ncss), Ncss := 0]
    #  info[is.na(NS), NS := 0]
    #  info[is.na(NL), NL := 0]
    #  info[is.na(sorted_Ncss1), sorted_Ncss1 := 0]
    #  info[is.na(sorted_Ncss2), sorted_Ncss2 := 0]

    wheatchr1 = wheatchr.copy().rename(columns={"chr": "sorted_chr", "alphachr": "sorted_alphachr"})
    # combined_stats.loc[:, "sorted_alphachr"] = pd.Categorical(combined_stats["sorted_alphachr"])
    combined_stats = wheatchr1.merge(combined_stats, how="right", on="sorted_alphachr")
    wheatchr2 = wheatchr.copy().rename(columns={"chr": "sorted_chr2", "alphachr": "sorted_alphachr2"})
    # combined_stats.loc[:, "sorted_alphachr2"] = pd.Categorical(combined_stats["sorted_alphachr2"])
    combined_stats = wheatchr2.merge(combined_stats, how="right", on="sorted_alphachr2")
    info = pd.merge(combined_stats, fai, on="scaffold_index", how="right").drop("scaffold", axis=1)
    for col in ["Ncss", "NS", "NL", "sorted_Ncss1", "sorted_Ncss2"]:
        info.loc[:, col] = info[col].fillna(0)
    return info
