import pandas as pd
import re
import dask.dataframe as dd
import joblib


def read_paf(fname, primary_only=True, save=False):
    names = dict(_ for _ in enumerate(["query", "query_length", "query_start", "query_end",
             "orientation", "reference", "reference_length",
             "reference_start", "reference_end", "matches", "alnlen", "mapq", "type"]))
    paf = dd.read_csv(fname, header=None, delimiter="\t").loc[:, :len(names)].rename(columns=names)
    paf["orientation"] = paf["orientation"].map({"+": 1, "-": -1})
    paf[["query_start", "query_end", "reference_start", "reference_end"]] += 1  # ???
    paf["type"] = paf["type"].str.replace("tp:A:", "")
    if primary_only is True:
        paf = paf[paf["type"] == "P"]
    if save:
        joblib.dump(paf, re.sub(r"\.paf(.{0, 0}|\.gz)$", "pkl", fname))
    return paf


# def summarize_paf(paf: pd.DataFrame):
#
#
#
# # Aggregate alignments for each (reference, query) pair
# summarize_paf<-function(paf){
#  setnames(dcast(
#      paf[, .(l=sum(alnlen)), key=.(query, reference, orientation)], query + reference ~ orientation, value.var="l", fill=0
# ), c("-1", "1"), c("lenrev", "lenfw"))->fr
#  paf[, .(query_start=min(query_start), query_end=max(query_end),
#        reference_start=min(reference_start), reference_end=max(reference_end),
#        matches=sum(matches), alnlen=sum(alnlen), naln=.N),
#     key=.(query, query_length, reference, reference_length)]->paf_summary
#  fr[paf_summary, on=c("query", "reference")]->paf_summary
#  setorder(paf_summary, -alnlen)
#  paf_summary[, idx := paste(1:.N)][]
# }