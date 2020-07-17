import dask.dataframe as dd
import time
import logging
dask_logger = logging.getLogger("dask")


def find_wrong_assignments(anchored_css: dd.DataFrame, measure: list,
                           sorted_percentile, popseq_percentile,
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

    # scaffold_index is the index
    dask_logger.debug("%s Find wrong assignments - starting", time.ctime())
    
    to_melt = anchored_css[measure].reset_index(drop=False)
    melted = dd.melt(to_melt, id_vars="scaffold_index",
                     value_vars=measure, var_name="map", value_name="chr").dropna()
    #  w[, .N, key=.(scaffold, chr)]->w
    assert "scaffold_index" in melted.columns
    assert "chr" in melted.columns
    melted = melted.groupby(["scaffold_index", "chr"]).size().to_frame("N").reset_index(drop=False)
    dask_logger.debug("\n\n")
    dask_logger.debug("%s Find wrong assignments - melted and grouped by, aggregating", time.ctime())    
    # Nchr_ass: number of assignments for the scaffold.
    # Nchr_ass_uniq: number of unique assignments for the scaffold.
    melted = melted.groupby("scaffold_index").agg({"N": ["sum", "size"]})
    melted.columns = ["Nchr_ass", "Nchr_ass_uniq"]
    dask_logger.debug("\n\n")
    assert melted.index.name == anchored_css.index.name == "scaffold_index", (
        melted.index.name, anchored_css.index.name)
    anchored_css = dd.merge(anchored_css, melted, on="scaffold_index", how="left")
    assert isinstance(anchored_css, dd.DataFrame), (type(anchored_css))
    dask_logger.debug("%s Find wrong assignments - finding the \"sorted threshold\"", time.ctime())
    assert "sorted_p12" in anchored_css.columns, anchored_css.columns
    assert "Ncss" in anchored_css.columns, anchored_css.columns
    sorted_threshold = anchored_css.query("Ncss >= 30").compute()
    sorted_threshold = sorted_threshold["sorted_p12"]
    sorted_threshold = sorted_threshold.quantile((sorted_percentile + 1) / 100)
    dask_logger.debug("%s Find wrong assignments - \"Sorted threshold\": %s", time.ctime(),
                        sorted_threshold)
    keys = ["Ncss", "sorted_p12", "sorted_Ncss2"]

    dask_logger.debug("%s Find wrong assignments - assigning the bad_sorted column", time.ctime())
    bad_sorted = anchored_css[["Ncss", "sorted_p12", "sorted_Ncss2"]]
    bad_sorted = bad_sorted.eval("Ncss >= 30 & sorted_p12 >= @sorted_threshold & sorted_Ncss2 >= 2",
                                 local_dict={"sorted_threshold": sorted_threshold})
    anchored_css["bad_sorted"] = bad_sorted.to_dask_array()
    try:
        pop_threshold = anchored_css.query("popseq_Ncss >= 30")["popseq_p12"].compute().quantile(
            (popseq_percentile + 1) / 100)
    except KeyError as kerror:
        raise KeyError((kerror, anchored_css.columns))
    keys = ["popseq_Ncss", "popseq_p12", "popseq_Ncss2"]
    dask_logger.debug("%s Find wrong assignments - assigning the bad_popseq column", time.ctime())
    bad_popseq = anchored_css[keys].eval("popseq_Ncss >= 30 & popseq_p12 >= @pop_threshold & popseq_Ncss2 >= 2",
                                         local_dict={"pop_threshold": pop_threshold})
    
    anchored_css["bad_popseq"] = bad_popseq.to_dask_array()
    measure = ["bad_sorted", "bad_popseq"]
    if hic is True:
        keys = ["Nhic", "hic_p12", "hic_N2"]
        dask_logger.debug("%s Find wrong assignments - assigning the bad_hic column", time.ctime())
        hic_threshold = anchored_css.query("Nhic >= 30")["hic_p12"].compute().quantile((hic_percentile + 1) / 100)
        bad_hic = anchored_css[keys].eval("Nhic >= 30 & hic_p12 >= @hic_threshold & hic_N2 >= 2",
                                          local_dict={"hic_threshold": hic_threshold})
        anchored_css["bad_hic"] = bad_hic.to_dask_array()
        measure.append("bad_hic")

    # Finally, let's count how many bad assignments we found for each scaffold. This can be between 0 and 3.
    # Melt only the columns we are interested in
    dask_logger.debug("%s Find wrong assignments - melting anchored_css and taking the measure", time.ctime())
    to_melt = anchored_css[measure].reset_index(drop=False)
    melted = dd.melt(to_melt,
                id_vars=["scaffold_index"],
                value_vars=measure,
                value_name="bad",
                var_name="map").dropna().query("bad == True")
    dask_logger.debug("%s Find wrong assignments - merging back onto anchored_css", time.ctime())    
    melted = melted.groupby("scaffold_index").size().to_frame("Nbad")
    anchored_css = dd.merge(anchored_css, melted, on="scaffold_index", how="left")
    dask_logger.debug("%s Find wrong assignments - finished", time.ctime())
    return anchored_css
