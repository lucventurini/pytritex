import pandas as pd
from ..utils import first, second_agg
import numpy as np
import dask.dataframe as dd
from dask.distributed import Client
from dask.delayed import delayed


def assign_carma(cssaln: dd.DataFrame, fai: dd.DataFrame, wheatchr: pd.DataFrame, client: Client):
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
    anchored_css = cssaln[~cssaln["sorted_alphachr"].isna()]
    combined_stats = anchored_css.groupby(by=["scaffold_index", "sorted_alphachr"]).size().to_frame("N")
    # anchored_css_grouped = anchored_css.reset_index(drop=False).groupby(["scaffold_index", "sorted_alphachr"])
    # combined_stats = anchored_css_grouped.compute().reset_index(level=1)
    combined_stats = combined_stats.reset_index(drop=False)
    # Now only keep the first two by "N"
    nsum = combined_stats.groupby("scaffold_index")["N"].sum().to_frame("Ncss")

    combined_stats = combined_stats[["scaffold_index", "N", "sorted_alphachr"]].groupby("scaffold_index").apply(
        lambda df: df.nlargest(2, "N"), meta=({"scaffold_index": int, "N": int, "sorted_alphachr": int})
    ).reset_index(drop=True)

    # Now assign first, second
    combined_stats = combined_stats.groupby("scaffold_index").agg(
        {"N": ["first", second_agg],
         "sorted_alphachr": ["first", second_agg]}
    )
    combined_stats.columns = ["sorted_Ncss1", "sorted_Ncss2", "sorted_alphachr", "sorted_alphachr2"]
    combined_stats["sorted_alphachr2"] = combined_stats["sorted_alphachr2"].astype(object)
    combined_stats = dd.merge(combined_stats, nsum, on="scaffold_index")

    combined_stats["sorted_pchr"] = combined_stats["sorted_Ncss1"].div(combined_stats["Ncss"], fill_value=0)
    combined_stats["sorted_p12"] = combined_stats["sorted_Ncss2"].div(combined_stats["sorted_Ncss1"], fill_value=0)
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

    short_arm_counts = anchored_css.query("(sorted_arm == 'S')").groupby(
            ["scaffold_index", "sorted_alphachr"])["sorted_arm"].size().to_frame("NS").reset_index(drop=False)

    long_arm_counts = anchored_css.query("(sorted_arm == 'L')").groupby(
        ["scaffold_index", "sorted_alphachr"])["sorted_arm"].size().to_frame("NL").reset_index(drop=False)

    func = delayed(dd.merge)(client.scatter(combined_stats),
                             long_arm_counts, on=["scaffold_index", "sorted_alphachr"],
                             how="left")
    func2 = delayed(dd.merge)(func, short_arm_counts, on=["scaffold_index", "sorted_alphachr"],
                             how="left")
    combined_stats = client.compute(func2).result()

    # combined_stats = short_arm_counts.merge(
    #     long_arm_counts.merge(
    #         combined_stats.reset_index(drop=False), left_index=True,
    #         right_on=long_arm_counts.index.names, how="right"),
    #     left_index=True, right_on=short_arm_counts.index.names, how="right")
    # combined_stats.loc[:, "NL"] = pd.to_numeric(combined_stats["NL"].fillna(0), downcast="signed")
    # combined_stats.loc[:, "NS"] = pd.to_numeric(combined_stats["NS"].fillna(0), downcast="signed")
    combined_stats["sorted_arm"] = (combined_stats["NS"] > combined_stats["NL"])
    combined_stats["sorted_arm"] = combined_stats["sorted_arm"].map({True: "S", False: "L"})

    combined_stats["sorted_arm"] = combined_stats["sorted_arm"].mask(
        combined_stats["sorted_arm"].isin(["1H", "3B"]), np.nan)

    combined_stats["sorted_parm"] = combined_stats["NL"].mask(
        combined_stats.sorted_arm == "L", combined_stats.NS).div(
        combined_stats["sorted_Ncss1"], fill_value=0).mask(combined_stats.sorted_arm.isna(), np.nan)

    #  setnames(copy(wheatchr), c("sorted_alphachr", "sorted_chr"))[z, on="sorted_alphachr"]->z
    #  setnames(copy(wheatchr), c("sorted_alphachr2", "sorted_chr2"))[z, on="sorted_alphachr2"]->z
    #  z[fai, on="scaffold"]->info
    #  info[is.na(Ncss), Ncss := 0]
    #  info[is.na(NS), NS := 0]
    #  info[is.na(NL), NL := 0]
    #  info[is.na(sorted_Ncss1), sorted_Ncss1 := 0]
    #  info[is.na(sorted_Ncss2), sorted_Ncss2 := 0]

    assert "sorted_alphachr" in combined_stats.columns
    assert "sorted_alphachr2" in combined_stats.columns

    info1 = delayed(dd.merge)(
        client.scatter(combined_stats),
        wheatchr.copy().rename(columns={"chr": "sorted_chr", "alphachr": "sorted_alphachr"}),
        how="left", on="sorted_alphachr")

    info2 = delayed(dd.merge)(info1,
        wheatchr.copy().rename(columns={"chr": "sorted_chr2", "alphachr": "sorted_alphachr2"}),
        how="left", on="sorted_alphachr2")

    info3 = delayed(dd.merge)(
        info2,
        client.scatter(fai), on="scaffold_index", how="right")
    info = client.compute(info3).result().drop("scaffold", axis=1)
    if info.index.name != "scaffold_index":
        info = info.set_index("scaffold_index")

    return info
