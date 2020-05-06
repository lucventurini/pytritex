import pandas as pd
import numpy as np


def _condition_checker(i: int, df: pd.DataFrame, el: pd.DataFrame, maxiter: int):
    i += 1
    df = df.sort_values("bin")
    x = pd.DataFrame().assign(
        old_node1=df.loc[range(df.shape[0] - 1)].reset_index()["cluster"],
        cluster=df.loc[range(1, df.shape[0])].reset_index()["cluster"],
        old_node2=df.loc[range(2, df.shape[0] + 1)].reset_index()["cluster"]
    )
    ee = el[["cluster1", "cluster2", "weight"]].set_index(["cluster1", "cluster2"])
    key = ["old_node1", "old_node2"]
    x = pd.merge(
        ee.rename(columns={"cluster1": "old_node1", "cluster2": "old_node2"}), x,
        left_on=key, right_on=key).rename(columns={"weight": "new_edge3"})
    key = ["cluster", "old_node1"]
    x = pd.merge(
        ee.rename(columns={"old_node1": "cluster", "old_node2": "old_node1"}), x,
        left_on=key, right_on=key).rename(columns={"weight": "old_edge1"})
    x = x[~x["new_edge3"].isna()]
    t = pd.merge(el[["cluster1", "cluster2"]],
                 df, left_on="cluster2", right_on="cluster").sort_values(["cluster1", "bin"])
    def index_creator(group):
        if group.shape[0] > 1:
            group["dist"] = np.append(group["bin"][1:].values - group["bin"][:-1], np.nan)
        else:
            group["dist"] = np.nan
        return group
    t = t.groupby("cluster1", as_index=False).apply(index_creator)["dist"]
    idx = t[t["dist"] == 1].index
    t = pd.DataFrame().assign(
        cluster=t.loc[idx, "cluster"],
        new_node1=t.loc[idx, "cluster2"],
        new_node2=t.loc[idx + 1, "cluster2"]
    )
    ee = el[["cluster1", "cluster2", "weight"]].sort_values(["cluster1", "cluster2"])
    ee = ee.rename(columns={"cluster1": "cluster", "cluster2": "new_node1"})
    t = pd.merge(ee.set_index(["cluster", "new_node1"]),
             t.set_index(["cluster", "new_node1"]),
             left_index=True, right_index=True).reset_index(drop=False).rename(columns={"weight": "new_edge1"})
    ee = ee.rename(columns={"new_node1": "new_node2"})
    t = pd.merge(
        ee, t, left_on=["cluster", "new_node2"], right_on=["cluster", "new_node2"]
    ).rename(columns={"weight": "new_edge2"})
    t = pd.merge(
        ee.rename(columns={"cluster": "new_node1", "new_node2": "new_node1"}),
        t, left_on=["new_node1", "new_node2"], right_on=["new_node1", "new_node2"]
    ).rename({"weight": "old_edge3"})
    m = pd.merge(x.set_index("cluster"), t.set_index("cluster"), left_index=True, right_index=True)
    m = m[~m["old_node1"].isna()]
    m["diff"] = m[["new_edge_1", "new_edge2", "new_edge3"]].sum(axis=1) - m[
        ["old_edge_1", "old_edge2", "old_edge3"]].sum(axis=1)
    m = m[m["diff"] < 0].sort_values("diff")
    condition = (m.shape[0] > 0) and (i <= maxiter)
    return i, m, condition

def node_relocation(df, el, maxiter=100, verbose=True):
    i = 0
    if df.shape[0] > 2:
        i, m, condition = _condition_checker(i, df, el, maxiter=maxiter)
        while condition is True:
            ne = m.head(1)
            idx = pd.DataFrame().assign(cluster=df["cluster"])
            


# node_relocation<-function(df, el, maxiter=100, verbose=T){
#  i<-0
#  if(nrow(df) > 2){
#   while() {
#    ne<-head(data.frame(m), 1)
#    idx<-data.frame(cluster=df$cluster, idx=1:nrow(df))
#    invisible(lapply(c("cluster", "new_node1", "new_node2", "old_node1", "old_node2"), function(i) {
#     merge(ne, by.x=i, by.y="cluster", idx)->>ne
#     colnames(ne)[which(colnames(ne) == "idx")]<<-paste(sep="_", "idx", i)
#    }))
#    df[with(ne, {
#     min_new <- min(idx_new_node1, idx_new_node2)
#     max_new <- max(idx_new_node1, idx_new_node2)
#     min_old <- min(idx_old_node1, idx_old_node2)
#     max_old <- max(idx_old_node1, idx_old_node2)
#     if(min_old < min_new) {
#      c(1:min_old, max_old:min_new, idx_cluster, max_new:nrow(df))
#     } else {
#      c(1:min_new, idx_cluster, max_new:min_old, max_old:nrow(df))
#     }
#    }),]->df
#    df$bin<-1:nrow(df)
#    df<-data.table(df)
#   }
#   if(verbose){
#    cat(paste0("Node relocation steps: ", i-1, "\n"))
#   }
#  }
#  df
# }