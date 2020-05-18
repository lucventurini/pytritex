import pandas as pd
from ..utils import first, second
import numpy as np


def anchor_css_alignment(cssaln: pd.DataFrame, fai: pd.DataFrame, wheatchr: pd.DataFrame):
    """
    This function will create a master table with the anchored information from the CSS alignment.
    :param cssaln: the original CSS alignment.
    :param fai: the original information on the scaffolds.
    :param species: the species we are working on (e.g. wheat).
    :return:
    """

    anchored_css = cssaln[~cssaln["sorted_alphachr"].isna()].copy()
    anchored_css = anchored_css.merge(
        anchored_css.groupby(["scaffold_index", "sorted_alphachr"], observed=True).size().to_frame("N"),
        how="left", on=["scaffold_index", "sorted_alphachr"])
    anchored_css = anchored_css.sort_values(
        ["scaffold_index", "N"], ascending=[True, False])[
        ["scaffold_index", "N", "sorted_alphachr"]].drop_duplicates().groupby("scaffold_index").agg(
        Ncss=("N", "count"), sorted_Ncss1=("N", first), sorted_Ncss2=("N", second),
        sorted_alphachr=("sorted_alphachr", first),
        sorted_alphachr2=("sorted_alphachr", second)
    ).reset_index(drop=False)
    anchored_css.loc[:, "sorted_pchr"] = anchored_css["sorted_Ncss1"] / anchored_css["Ncss"]
    anchored_css.loc[:, "sorted_p12"] = anchored_css["sorted_Ncss2"] / anchored_css["sorted_Ncss1"]

    # Assignment of CARMA chromosome arm
    c_al = cssaln.loc[cssaln["sorted_arm"] == "L"].groupby(
        ["scaffold_index", "sorted_alphachr"]).size().to_frame("NL")
    # cssaln[sorted_arm == "S",.(NS=.N), keyby =.(scaffold, sorted_alphachr)]->as
    c_as = cssaln.loc[cssaln["sorted_arm"] == "S"].groupby(
        ["scaffold_index", "sorted_alphachr"]).size().to_frame("NS")
    anchored_css = pd.merge(
        c_as,
        pd.merge(c_al, anchored_css.copy().set_index(
            ["scaffold_index", "sorted_alphachr"]), left_index=True, right_index=True),
        left_index=True, right_index=True).reset_index(drop=False)
    anchored_css.loc[anchored_css["NS"].isna(), "NS"] = 0
    anchored_css.loc[anchored_css["NL"].isna(), "NL"] = 0
    anchored_css.loc[:, "sorted_arm"] = anchored_css[["NS", "NL"]].apply(
        lambda row: "S" if row.NS > row.NL else "L", axis=1)
    anchored_css.loc[anchored_css["sorted_alphachr"].isin(("1H", "3B")), "sorted_arm"] = np.nan
    anchored_css.loc[anchored_css["sorted_arm"] == "S",
                     "sorted_parm"] = anchored_css["NS"] / anchored_css["sorted_Ncss1"]
    anchored_css.loc[anchored_css["sorted_arm"] == "L",
                     "sorted_parm"] = anchored_css["NL"] / anchored_css["sorted_Ncss1"]
    anchored_css = pd.merge(
        pd.merge(anchored_css, wheatchr.rename(columns={"popseq_alphachr": "sorted_alphachr",
                                                        "popseq_chr": "sorted_chr"}),
                 left_on="sorted_alphachr",
                 right_on="sorted_alphachr"),
        wheatchr.rename(
            columns={"popseq_alphachr": "sorted_alphachr2",
                     "popseq_chr": "sorted_chr2"}),
        left_on="sorted_alphachr2", right_on="sorted_alphachr2")
    anchored_css = pd.merge(anchored_css, fai, on="scaffold_index", how="right")
    for col in ["Ncss", "NS", "NL", "sorted_Ncss1", "sorted_Ncss2"]:
        anchored_css.loc[:, col] = anchored_css[col].fillna(0)
    anchored_css.drop("scaffold", axis=1, errors="ignore", inplace=True)
    return anchored_css
