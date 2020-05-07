import pandas as pd
import numpy as np
from .utils.chrnames import chrNames
from scipy.stats import median_absolute_deviation


def make_super_scaffolds(links, info: pd.DataFrame, excluded=(), ncores=1, prefix=None):
    info2 = info[["scaffold", "popseq_chr", "popseq_cM", "length"]][:].rename(columns={
        "popseq_chr": "chr", "popseq_cM": "cM"})
    excluded_scaffolds = excluded
    input_df = info2[info2.scaffold.isin(excluded_scaffolds)].rename(
        columns={"scaffold": "cluster"}
    )
    hl = links.copy().rename(columns={"scaffold1": "cluster1", "scaffold2": "cluster2"})
    s = make_super(hl=hl, cluster_info=input_df, verbose=False, prefix=prefix, cores=ncores,
                   paths=True, path_max=0, known_ends=False, maxiter=100)
    m = s["membership"][:].rename(columns={"cluster": "scaffold"})
    maxidx = s["super_info"]["super"].astype(str).str.replace("^{}_".format(prefix), "", regex=True).astype(int).max()
    binder= info[~info["scaffold"].isin(m["scaffold"])][:][["scaffold", "length" "popseq_chr", "popseq_cM"]].rename(
        columns={"popseq_chr": "chr", "popseq_cM": "cM"})
    binder["bin"], binder["rank"], binder["backbone"] = 1, 0, True
    binder["excuded"] = binder["scaffold"].isin(excluded_scaffolds)
    binder["super"] = pd.Series(list(range(maxidx + 1, maxidx + 1 + binder.shape[0]))).astype(str).replace(
        "^(.)", "{}_\\1".format(prefix), regex=True)
    m = m.append(binder)
    res = m[:].groupby("super").agg({
    "bin": ["count", "max"], "rank": "max", "length": "sum"}).rename(
        columns={("bin", "count"): "n", ("bin", "max"): "nbin", "rank": "max_rank"}
    ).reset_index()
    m = pd.merge(m, res[["super", "n", "nbin"]].rename(columns={"n": "super_size", "nbin": "super_nbin"}),
                 left_on="super", right_on="super")
    return {"membership": m, "info": res}


def anchor_scaffolds(assembly: dict,
                     popseq,
                     species=None,
                     sorted_percentile=95,
                     popseq_percentile=90,
                     hic_percentile=98):
    if species is None:
        raise KeyError(
            "Parameter 'species' is NULL. Please set 'species' to one of"
"\"wheat\", \"barley\", \"oats\", \"lolium\", \"sharonensis\" or \"rye\".")
    elif species not in ("wheat", "barley", "rye", "oats", "sharonensis", "lolium"):
        raise KeyError(
            "Parameter 'species' is not valid. Please set 'species' to one of"
"\"wheat\", \"barley\", \"oats\", \"lolium\", \"sharonensis\" or \"rye\"."
        )

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
               "popseq_cM", "popseq_cM_sd", "popseq_cM_mad"]].merge(x, on=["scaffold", "popseq_alphachr"])
    zz.loc[:, "popseq_pchr"] = zz["popseq_Ncss1"] / zz["popseq_Ncss"]
    zz.loc[:, "popseq_p12"] = zz["popseq_Ncss2"] / zz["popseq_Ncss1"]
    zz = wheatchr.rename(columns={"popseq_chr": "popseq_chr2",
                                  "popseq_alphachr": "popseq_alphachr2"}).merge(
        wheatchr.merge(zz, on="popseq_alphachr"), on="popseq_alphachr2")
    info = zz.merge(info, on="scaffold")
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
        z = tcc_pos[~tcc_pos["chr1"].isna()].rename(
            {"scaffold2": "scaffold", "chr1": "hic_chr"}
        ).groupby(["scaffold", "hic_chr"]).size().to_frame("N").reset_index(drop=False)
        zgrouped = z.sort_values(["scaffold", "N"], ascending=[True, False]).groupby("scaffold")
        zNsum = zgrouped["N"].sum()
        z_hic_chr = zgrouped.agg({
            "hic_chr": [lambda series: series.iloc[0], lambda series: np.nan if series.shape[0] == 1 else series.iloc[1]],
            "N": [lambda series: series.iloc[0], lambda series: 0 if series.shape[0] == 1 else series.iloc[1]]
        })
        zz = pd.merge(zNsum, z_hic_chr, left_index=True, right_index=True)
        zz.columns = ["Nhic", "hic_chr", "hic_chr2", "hic_N1", "hic_N2"]
        zz["hic_pchr"] = zz["hic_N1"] / zz["Nhic"]
        zz["hic_p12"] = zz["hic_N2"] / zz["hic_N1"]
        info = zz.reset_index(drop=False).merge(info, on="scaffold")
        info.loc[info["Nhic"].isna(), "Nhic"] = 0
        info.loc[info["hic_N1"], "hic_N1"] = 0
        info.loc[info["hic_N2"], "hic_N2"] = 0

        measure = ["popseq_chr", "hic_chr", "sorted_chr"]
    else:
        measure = ["popseq_chr", "sorted_chr"]

    # Now melting. Let's hope something remains ...
    w = pd.melt(
        info, id_vars=["scaffold"],
        value_vars=measure,
        var_name = "map",
        value_name = "chr",
    )
    w = w.groupby(["scaffold", "chr"]).size().to_frame("N").reset_index(drop=False)


melt(info, id.vars="scaffold", measure.vars=measure, variable.factor=F, variable.name="map", na.rm=T, value.name="chr")->w
w[, .N, key=.(scaffold, chr)]->w
w[order(-N), .(Nchr_ass = sum(N), Nchr_ass_uniq =.N), keyby=scaffold]->w
w[info, on="scaffold"]->info
info[is.na(Nchr_ass), Nchr_ass := 0]
info[is.na(Nchr_ass_uniq), Nchr_ass_uniq := 0]

x < -info[Ncss >= 30, quantile(na.omit(sorted_p12), 0:100 / 100)][sorted_percentile + 1]
info[, bad_sorted := (sorted_p12 >= x & sorted_Ncss2 >= 2)]
x < -info[popseq_Ncss >= 30, quantile(na.omit(popseq_p12), 0:100 / 100)][popseq_percentile + 1]
info[, bad_popseq := (popseq_p12 >= x & popseq_Ncss2 >= 2)]

info[ is.na(bad_sorted), bad_sorted := F]
info[ is.na(bad_popseq), bad_popseq := F]

if (hic){
x < -info[Nhic >= 30, quantile(na.omit(hic_p12), 0:100 / 100)][hic_percentile + 1]
info[Nhic >= 30, bad_hic := hic_p12 >= x & hic_N2 >= 2]
info[ is.na(bad_hic), bad_hic := F]
}

melt(info, id.vars = "scaffold", measure.vars = grep(value=T, "bad_", names(
    info)), variable.factor = F, variable.name = "map", na.rm = T, value.name = "bad")[bad == T]->w
w[,.(Nbad=.N), key = scaffold]->w
w[info, on = "scaffold"]->info
info[ is.na(Nbad), Nbad := 0]

assembly$info < - info
assembly$popseq < - popseq
if (hic)
{
    assembly$fpairs < - tcc_pos
}
assembly
}