import pandas as pd
import numpy as np


def find_wrong_assignments(anchored_css: pd.DataFrame, measure: list, sorted_percentile, popseq_percentile,
                           hic_percentile, hic=False):

    """This function will find those scaffolds for which the assignment smells "fishy", ie
    scaffolds where the ratio of second-best to best is over the defined percentiles."""

    # melt(info, id.vars="scaffold", measure.vars=measure, variable.factor=F, variable.name="map", na.rm=T, value.name="chr")->w
    #  w[, .N, key=.(scaffold, chr)]->w
    #  w[order(-N), .(Nchr_ass = sum(N), Nchr_ass_uniq = .N), keyby=scaffold]->w
    #  w[info, on="scaffold"]->info
    #  info[is.na(Nchr_ass), Nchr_ass := 0]
    #  info[is.na(Nchr_ass_uniq), Nchr_ass_uniq := 0]
    #
    #  x<-info[Ncss >= 30, quantile(na.omit(sorted_p12), 0:100/100)][sorted_percentile+1]
    #  info[, bad_sorted := (sorted_p12 >= x & sorted_Ncss2 >= 2)]
    #  x<-info[popseq_Ncss >= 30, quantile(na.omit(popseq_p12), 0:100/100)][popseq_percentile+1]
    #  info[, bad_popseq := (popseq_p12 >= x & popseq_Ncss2 >= 2)]
    #
    #  info[is.na(bad_sorted), bad_sorted := F]
    #  info[is.na(bad_popseq), bad_popseq:= F]
    #
    #  if(hic){
    #   x<-info[Nhic >= 30, quantile(na.omit(hic_p12), 0:100/100)][hic_percentile+1]
    #   info[Nhic >= 30, bad_hic := hic_p12 >= x & hic_N2 >= 2]
    #   info[is.na(bad_hic), bad_hic := F]
    #  }
    #
    #  melt(info, id.vars="scaffold", measure.vars=grep(value=T, "bad_", names(info)), variable.factor=F, variable.name="map", na.rm=T, value.name="bad")[bad == T]->w
    #  w[, .(Nbad=.N), key=scaffold]->w
    #  w[info, on="scaffold"]->info
    #  info[is.na(Nbad), Nbad := 0]

    melted = pd.melt(anchored_css, id_vars="scaffold_index",
                     value_vars=measure, var_name="map", value_name="chr").dropna()
    #  w[, .N, key=.(scaffold, chr)]->w
    melted = melted.groupby(["scaffold_index", "chr"]).size().to_frame("N").reset_index(drop=False)
    # Nchr_ass: number of assignments for the scaffold.
    # Nchr_ass_uniq: number of unique assignments for the scaffold.
    melted = melted.sort_values("N", ascending=False).groupby("scaffold_index").agg(
        Nchr_ass=("N", "sum"), Nchr_ass_uniq=("N", "size"))
    anchored_css = melted.merge(anchored_css, right_on="scaffold_index", left_index=True, how="right")
    anchored_css.loc[:, "scaffold_index"] = pd.to_numeric(anchored_css["scaffold_index"].fillna(0),
                                                          downcast="signed")
    anchored_css.loc[:, "Nchr_ass"] = pd.to_numeric(anchored_css["Nchr_ass"].fillna(0),
                                                    downcast="signed")
    anchored_css.loc[:, "Nchr_ass_uniq"] = pd.to_numeric(anchored_css["Nchr_ass_uniq"].fillna(0),
                                                         downcast="signed")
    sorted_threshold = anchored_css.loc[
        anchored_css["Ncss"] >= 30, "sorted_p12"].quantile((sorted_percentile + 1) / 100)
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
                anchored_css["hic_N2"] >= 2)).astype(bool)
        anchored_css.loc[:, "bad_hic"] = anchored_css["bad_hic"].fillna(False).astype(bool)
        measure.append("bad_hic")

    # Finally, let's count how many bad assignments we found for each scaffold. This can be between 0 and 3.
    melted = pd.melt(anchored_css,
                id_vars=["scaffold_index"],
                value_vars=measure,
                value_name="bad",
                var_name="map").dropna().loc[lambda df: df["bad"] == True, :]
    melted = melted.groupby("scaffold_index").size().to_frame("Nbad").reset_index(drop=False)
    anchored_css = melted.merge(anchored_css, on="scaffold_index", how="right")
    anchored_css.loc[:, "Nbad"] = pd.to_numeric(anchored_css["Nbad"].fillna(0),
                                                downcast="signed")
    return anchored_css
