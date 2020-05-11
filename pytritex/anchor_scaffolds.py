import pandas as pd
import numpy as np
from .utils.chrnames import chrNames
from scipy.stats import median_absolute_deviation


def anchor_scaffolds(assembly: dict,
                     popseq,
                     species=None,
                     sorted_percentile=95,
                     popseq_percentile=90,
                     hic_percentile=98):
    if species is None:
        raise KeyError(
            "Parameter 'species' is NULL. Please set 'species' to one of "
            "\"wheat\", \"barley\", \"oats\", \"lolium\", \"sharonensis\" or \"rye\".")
    elif species not in ("wheat", "barley", "rye", "oats", "sharonensis", "lolium"):
        raise KeyError(
            "Parameter 'species' is not valid. Please set 'species' to one of "
            "\"wheat\", \"barley\", \"oats\", \"lolium\", \"sharonensis\" or \"rye\".")

    fai = assembly["info"]
    cssaln = assembly["cssaln"]
    if "fpairs" not in assembly:
        fpairs = None
        hic = False
    else:
        fpairs = assembly["fpairs"]
        hic = (assembly["fpairs"].shape[0] > 0)
    # cssaln[! is.na(sorted_alphachr),.N, keyby =.(scaffold, sorted_alphachr)]->z
    # For each scaffold in the assembly, for each chromosomal bin: how many CSS contigs confirm each possible bin?
    z = cssaln[~cssaln["sorted_alphachr"].isna()]
    z = z.merge(
        z.groupby(["scaffold", "sorted_alphachr"]).size().to_frame("N").reset_index(drop=False),
        on=["scaffold", "sorted_alphachr"])
    zgrouped = z.sort_values(["scaffold", "N"], ascending=[True, False]).groupby("scaffold")
    z = zgrouped.agg(
        {"N": [np.sum, lambda series: series.loc[series.index[0]], lambda series: series.loc[series.index[1]]],
         "sorted_alphachr": [lambda series: series.loc[series.index[0]], lambda series: series.loc[series.index[1]]]
         }).reset_index(drop=False)
    z.columns = ["scaffold", "Ncss", "sorted_Ncss1", "sorted_Ncss2", "sorted_alphachr", "sorted_alphachr2"]
    z.loc[:, "sorted_pchr"] = z["sorted_Ncss1"] / z["Ncss"]
    z.loc[:, "sorted_p12"] = z["sorted_Ncss2"] / z["sorted_Ncss1"]

    # Assignment of CARMA chromosome arm
    c_al = cssaln.loc[cssaln["sorted_arm"] == "L"].groupby(["scaffold", "sorted_alphachr"]).size().to_frame("NL")
    # cssaln[sorted_arm == "S",.(NS=.N), keyby =.(scaffold, sorted_alphachr)]->as
    c_as = cssaln.loc[cssaln["sorted_arm"] == "S"].groupby(["scaffold", "sorted_alphachr"]).size().to_frame("NS")
    z = pd.merge(c_as,
             pd.merge(c_al, z.set_index(["scaffold", "sorted_alphachr"]), left_index=True, right_index=True),
             left_index=True, right_index=True).reset_index(drop=False)
    z.loc[z["NS"].isna(), "NS"] = 0
    z.loc[z["NL"].isna(), "NL"] = 0
    z.loc[:, "sorted_arm"] = z[["NS", "NL"]].apply(lambda row: "S" if row.NS > row.NL else "L", axis=1)
    z.loc[z["sorted_alphachr"].isin(("1H", "3B")), "sorted_arm"] = np.nan
    z.loc[z["sorted_arm"] == "S", "sorted_parm"] = z["NS"] / z["sorted_Ncss1"]
    z.loc[z["sorted_arm"] == "L", "sorted_parm"] = z["NL"] / z["sorted_Ncss1"]
    wheatchr = chrNames(species=species).rename(
        columns={"alphachr": "popseq_alphachr", "chr": "popseq_chr"}
    )
    z = pd.merge(
        pd.merge(z, wheatchr.rename(columns={"popseq_alphachr": "sorted_alphachr",
                                             "popseq_chr": "sorted_chr"}),
                 left_on="sorted_alphachr",
                 right_on="sorted_alphachr"),
        wheatchr.rename(
            columns={"popseq_alphachr": "sorted_alphachr2",
                     "popseq_chr": "sorted_chr2"}),
        left_on="sorted_alphachr2", right_on="sorted_alphachr2")
    info = pd.merge(z, fai, on="scaffold")
    for column in ["Ncss", "NS", "NL", "sorted_Ncss1", "sorted_Ncss2"]:
        info.loc[info[column].isna(), column] = 0
    z = pd.merge(popseq.loc[~popseq["sorted_alphachr"].isna(), ["css_contig", "popseq_alphachr", "popseq_cM"]],
                 cssaln.loc[:, ["css_contig", "scaffold"]],
                 on="css_contig", how="left")
    z.loc[z["scaffold"].isna(), "scaffold"] = 0
    zgrouped = z.groupby(["scaffold", "popseq_alphachr"])
    zN = zgrouped.size().to_frame("N")
    zother = zgrouped.agg({"popseq_cM": [np.mean, np.std, median_absolute_deviation]})
    zother.columns = zother.columns.to_flat_index()
    zother = zother.rename(columns={zother.columns[0]: "popseq_cM", ("popseq_cM", "std"): "popseq_cM_sd",
                                    zother.columns[2]: "popseq_cM_mad"})
    zz = pd.merge(zN, zother, left_index=True, right_index=True).reset_index(drop=False)
    zz.loc[:, "popseq_Ncss"] = zz.groupby("scaffold")["N"].transform(np.sum)
    # Create a different dataframe now, we'll merge soon
    x = zz.sort_values(["scaffold", "N"], ascending=[True, False]).groupby("scaffold").agg(
        {"N": [lambda series: series.iloc[0],
               lambda series: 0 if series.shape[0] == 1 else series.iloc[1]],
         "popseq_alphachr": [lambda series: series.iloc[0],
                             lambda series: 0 if series.shape[0] == 1 else series.iloc[1]]
         })
    x.columns = ["popseq_Ncss1", "popseq_Ncss2", "popseq_alphachr", "popseq_alphachr2"]
    x.reset_index(drop=False, inplace=True)
    zz = zz.loc[:, ["scaffold", "popseq_alphachr", "popseq_Ncss",
               "popseq_cM", "popseq_cM_sd", "popseq_cM_mad"]].merge(x, on=["scaffold", "popseq_alphachr"],
                                                                    how="right")
    zz.loc[:, "popseq_pchr"] = zz["popseq_Ncss1"] / zz["popseq_Ncss"]
    zz.loc[:, "popseq_p12"] = zz["popseq_Ncss2"] / zz["popseq_Ncss1"]
    zz = wheatchr.rename(columns={"popseq_chr": "popseq_chr2",
                                  "popseq_alphachr": "popseq_alphachr2"}).merge(
        wheatchr.merge(zz, on="popseq_alphachr"), on="popseq_alphachr2")
    info = zz.merge(info, on="scaffold", how="right")
    info.loc[info["popseq_Ncss"].isna(), "popseq_Ncss"] = 0
    info.loc[info["popseq_Ncss1"].isna(), "popseq_Ncss1"] = 0
    info.loc[info["popseq_Ncss2"].isna(), "popseq_Ncss2"] = 0
    # # Assignment of POPSEQ genetic positions
    if hic:
        info0 = info[info["popseq_chr"] == info["sorted_chr"]][["scaffold", "popseq_chr"]].rename(
            columns={"popseq_chr": "chr"})
        # TODO: complete the renaming
        tcc_pos = info0.rename(columns={"chr": "chr1", "scaffold": "scaffold1"}).merge(fpairs, on="scaffold1")
        tcc_pos = info0.rename(columns={"chr": "chr2", "scaffold": "scaffold2"}).merge(tcc_pos, on="scaffold2")
        # tcc_pos[!is.na(chr1), .N, key=.(scaffold=scaffold2, hic_chr=chr1)]->z
        z = tcc_pos[~tcc_pos["chr1"].isna()].rename(
            columns={"scaffold2": "scaffold", "chr1": "hic_chr"}
        ).groupby(["scaffold", "hic_chr"]).size().to_frame("N").reset_index(drop=False)
        zgrouped = z.sort_values(["scaffold", "N"], ascending=[True, False]).groupby("scaffold")
        zz = zgrouped.agg(
            hic_chr=("hic_chr", lambda series: series.iloc[0]),
            hic_chr2=("hic_chr", lambda series: np.nan if series.shape[0] == 1 else series.iloc[1]),
            hic_N1=("N", lambda series: series.iloc[0]),
            hic_N2=("N", lambda series: np.nan if series.shape[0] == 1 else series.iloc[1]),
            Nhic=("N", "sum")
        )
        zz["hic_pchr"] = zz["hic_N1"] / zz["Nhic"]
        zz["hic_p12"] = zz["hic_N2"] / zz["hic_N1"]
        info = zz.reset_index(drop=False).merge(info, on="scaffold", how="right")
        info.loc[info["Nhic"].isna(), "Nhic"] = 0
        info.loc[info["hic_N1"].isna(), "hic_N1"] = 0
        info.loc[info["hic_N2"].isna(), "hic_N2"] = 0
        measure = ["popseq_chr", "hic_chr", "sorted_chr"]
    else:
        measure = ["popseq_chr", "sorted_chr"]

    # Now melting.
    w = pd.melt(info,
        id_vars=["scaffold"],
        value_vars=measure,
        var_name="map",
        value_name="chr",
    )
    w = w.groupby(["scaffold", "chr"]).size().to_frame("N").reset_index(drop=False)
    # w[, .N, key=.(scaffold, chr)]->w
    wgrouped = w.sort_values("N", ascending=False).groupby("scaffold")
    __temp = wgrouped.agg({"N": [np.sum, lambda df: df.shape[0]]})
    __temp.columns = ["Nchr_ass", "Nchr_ass_uniq"]
    w = pd.merge(w.set_index("scaffold"), __temp, left_index=True, right_index=True, how="left").reset_index(drop=False)
    info = w.merge(info, on="scaffold", how="right")
    info.loc[info["Nchr_ass"].isna(), "Nchr_ass"] = 0
    info.loc[info["Nchr_ass_uniq"].isna(), "Nchr_ass_uniq"] = 0
    # Get this parameter, "x"?
    x = info.loc[info["Ncss"] >= 30, "sorted_p12"].quantile((sorted_percentile + 1) / 100)
    info.loc[:, "bad_sorted"] = (info["sorted_p12"] >= x) & (info["sorted_Ncss2"] >= 2)
    x = info.loc[info["popseq_Ncss"] >= 30, "popseq_p12"].quantile((popseq_percentile + 1) / 100)
    info.loc[:, "bad_popseq"] = (info["popseq_p12"] >= x) & (info["popseq_Ncss2"] >= 2)
    info.loc[info["bad_sorted"].isna(), "bad_sorted"] = False
    info.loc[info["bad_popseq"].isna(), "bad_popseq"] = False
    measure = ["bad_sorted", "bad_popseq"]

    if hic is not None:
        x = info.loc[info["Nhic"] >= 30, "hic_p12"].quantile((hic_percentile + 1) / 100)
        info.loc[info["Nhic"] >= 30, "bad_hic"] = ((info["hic_p12"] >=x ) & (info["hic_N2"] >= 2))
        info.loc[info["bad_hic"].isna(), "bad_hic"] = False
        measure.append("bad_hic")

    w = pd.melt(info,
                id_vars=["scaffold"],
                value_vars=measure,
                value_name="bad",
                var_name = "map").dropna()
    w = w.loc[w["bad"] == True, :]
    w = w.merge(w.groupby("scaffold").size().to_frame("Nbad").reset_index(drop=False),
                on="scaffold", how="left")
    info = w.merge(info, on="scaffold", how="right")
    info.loc[info["Nbad"].isna(), "Nbad"] = 0
    assembly["info"] = info
    assembly["popseq"] = popseq
    if hic is not None:
        assembly["fpairs"] = tcc_pos
    return assembly
