import pandas as pd
from ..utils import first, second_agg, second
import numpy as np
import dask.dataframe as dd
from dask.distributed import Client
from dask.delayed import delayed
import time
import logging
dask_logger = logging.getLogger("dask")


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
    anchored_css = cssaln[~cssaln["sorted_chr"].isna()]
    dask_logger.debug("%s Assigning CARMA - calculating N", time.ctime())
    combined_stats = anchored_css.astype({"sorted_chr": int}).groupby(
        by=["scaffold_index", "sorted_chr"]).size().to_frame("N")
    combined_stats = client.persist(combined_stats)
    # anchored_css_grouped = anchored_css.reset_index(drop=False).groupby(["scaffold_index", "sorted_alphachr"])
    # combined_stats = anchored_css_grouped.compute().reset_index(level=1)
    combined_stats = combined_stats.reset_index(drop=False)
    # Now only keep the first two by "N"
    dask_logger.debug("%s Assigning CARMA - calculating Ncss", time.ctime())    
    nsum = combined_stats.groupby("scaffold_index")["N"].sum().to_frame("Ncss")
    dask_logger.debug("%s Assigning CARMA - calculating Ncss1, Ncss2", time.ctime())
    combined_stats = combined_stats[["scaffold_index", "N", "sorted_chr"]].compute()
    combined_stats = combined_stats.sort_values("N", ascending=False).groupby("scaffold_index").head(2)
    combined_stats = combined_stats.groupby("scaffold_index").agg({"N": ["first", second], "sorted_chr": ["first", second]})
    # Back to dask
    combined_stats = dd.from_pandas(combined_stats, npartitions=anchored_css.npartitions)
    # Now assign first, second
    combined_stats.columns = ["sorted_Ncss1", "sorted_Ncss2", "sorted_chr", "sorted_chr2"]
    func = delayed(dd.merge)(client.scatter(combined_stats), nsum, on="scaffold_index")
    combined_stats = client.compute(func).result()
    dask_logger.debug("%s Assigning CARMA - percentages pchr, p12", time.ctime())    
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

    dask_logger.debug("%s Assigning CARMA - calculating arm counts", time.ctime())
    short_arm_counts = anchored_css.query("(sorted_arm == 'S')").groupby(
            ["scaffold_index", "sorted_chr"])["sorted_arm"].size().to_frame("NS").reset_index(drop=False)

    long_arm_counts = anchored_css.query("(sorted_arm == 'L')").groupby(
        ["scaffold_index", "sorted_chr"])["sorted_arm"].size().to_frame("NL").reset_index(drop=False)

    dask_logger.debug("%s Assigning CARMA - merging arm counts back", time.ctime())
    func = delayed(dd.merge)(client.scatter(combined_stats),
                             long_arm_counts, on=["scaffold_index", "sorted_chr"], how="left")
    func2 = delayed(dd.merge)(func, short_arm_counts, on=["scaffold_index", "sorted_chr"], how="left")
    combined_stats = client.compute(func2).result()

    dask_logger.debug("%s Assigning CARMA - assigning arm", time.ctime())
    combined_stats["sorted_arm"] = (combined_stats["NS"] > combined_stats["NL"])
    combined_stats["sorted_arm"] = combined_stats["sorted_arm"].map({True: "S", False: "L"})

    h1_b3 = wheatchr.loc[wheatchr.alphachr.isin(["3B", "1H"]), "chr"].values
    combined_stats["sorted_arm"] = combined_stats["sorted_arm"].mask(combined_stats["sorted_chr"].isin(h1_b3), np.nan)

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

    assert "sorted_chr" in combined_stats.columns
    assert "sorted_chr2" in combined_stats.columns

    dask_logger.debug("%s Assigning CARMA - merging combined stats with wheatchr and FAI", time.ctime())
    combined_stats = combined_stats.astype({"sorted_chr": int, "sorted_chr2": float, "scaffold_index": int})
    combined_stats.index = combined_stats.index.astype(int)

    dask_logger.debug("%s Combined stats: %s; wheatchr: %s", time.ctime(),
                        combined_stats.shape[0].compute(), wheatchr.shape[0])
    info1 = delayed(dd.merge)(
        client.scatter(combined_stats),
        wheatchr.copy().rename(columns={"chr": "sorted_chr", "alphachr": "sorted_alphachr"}),
        how="left", on="sorted_chr")

    info2 = delayed(dd.merge)(info1,
        wheatchr.copy().rename(columns={"chr": "sorted_chr2", "alphachr": "sorted_alphachr2"}),
        how="left", on="sorted_chr2")

    info2 = client.compute(info2).result()
    assert info2.scaffold_index.dtype == fai.index.dtype
    info2 = info2.set_index("scaffold_index")
    
    dask_logger.debug("%s Finished with info2 (%s), now merging with FAI (%s)",
                        time.ctime(), info2.shape[0].compute(), fai.shape[0].compute())

    info = dd.merge(info2, fai, on="scaffold_index", how="right").drop("scaffold", axis=1)
    info = info.reset_index(drop=False).set_index("scaffold_index")
    dask_logger.debug("%s Finished with info (%s)", time.ctime(), info.shape[0].compute())
    assert isinstance(info, dd.DataFrame)
    if info.index.name != "scaffold_index":
        dask_logger.debug("%s Assigning CARMA - resetting the index to scaffold_index. Npartitions: %s", time.ctime(), info.npartitions)
        info = info.set_index("scaffold_index")
        dask_logger.debug("%s Assigning CARMA - reset the index", time.ctime())
    dask_logger.debug("%s Assigning CARMA - finished", time.ctime())
        
    return info
