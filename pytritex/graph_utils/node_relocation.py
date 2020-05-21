import pandas as pd
import networkit
import numpy as np


# # Traveling salesman heuristic for Hi-C mapping construction
# node_relocation<-function(df, el, maxiter=100, verbose=T){
#  i<-0
#  if(nrow(df) > 2){
#   while({
#    i<-i+1
#    df[order(bin)]->df
#    x<-data.table(old_node1=df[1:(nrow(df)-2)]$cluster, cluster=df[2:(nrow(df)-1)]$cluster, old_node2=df[3:nrow(df)]$cluster)
#    setkeyv(el[, .(cluster1, cluster2, weight)], c("cluster1", "cluster2"))->ee
#    setnames(setnames(ee, c("cluster1", "cluster2"), c("old_node1", "old_node2"))[setkeyv(x, c("old_node1", "old_node2"))], "weight", "new_edge3")->x
#    setnames(setnames(ee, c("old_node1", "old_node2"), c("cluster", "old_node1"))[setkeyv(x, c("cluster", "old_node1"))], "weight", "old_edge1")->x
#    setnames(setnames(ee, "old_node1", "old_node2")[setkeyv(x, c("cluster", "old_node2"))], "weight", "old_edge2")->x
#    x[!is.na(new_edge3) ]->x
#
#    setkey(el[, .(cluster1, cluster2)], "cluster2")[setkey(copy(df), "cluster")][order(cluster1, bin)]->t
#    which(t[, dist:={if(.N == 1) {as.integer(NA)} else { as.integer(c(bin[2:.N],NA)-bin)}}, by=cluster1]$dist == 1)->idx
#
#    data.table(cluster=t$cluster1[idx], new_node1=t$cluster2[idx], new_node2=t$cluster2[idx+1])->t
#    setkeyv(el[, .(cluster1, cluster2, weight)], c("cluster1", "cluster2"))->ee
#    setnames(setnames(ee, c("cluster1", "cluster2"), c("cluster", "new_node1"))[setkeyv(t, c("cluster", "new_node1"))], "weight", "new_edge1")->t
#    setnames(setnames(ee, "new_node1", "new_node2")[setkeyv(t, c("cluster", "new_node2"))], "weight", "new_edge2")->t
#    setnames(setnames(ee, c("cluster", "new_node2"), c("new_node1", "new_node2"))[setkeyv(t, c("new_node1", "new_node2"))], "weight", "old_edge3")->t
#    setkey(x, "cluster")[setkey(t, "cluster")][!is.na(old_node1)]->m
#    m[ ,diff := new_edge1+new_edge2+new_edge3 - old_edge1 - old_edge2 - old_edge3]
#    nrow(m<-m[diff < 0][order(diff)]) > 0 & i <= maxiter
#    }) {
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


def checker(df, edge_list):
    x = pd.DataFrame().assign(old_node1=df["cluster"][:-2],
                              cluster=df["cluster"][1:-1],
                              old_node2=df["cluster"][2:])
    # Extract the edge weights, indexed by ege
    ee = edge_list[["cluster1", "cluster2", "weight"]]
    # setnames(ee, c("cluster1", "cluster2"), c("old_node1", "old_node2"))
    __left = ee.copy().rename(columns={"cluster1": "old_node1", "cluster2": "old_node2"})
    # setnames(ee[setkeyv(x, c("old_node1", "old_node2"))], "weight", "new_edge3")->x
    x = __left.merge(x, on=["old_node1", "old_node2"], how="right").rename(columns={"weight": "new_edge3"})
    # setnames(ee, c("old_node1", "old_node2"), c("cluster", "old_node1"))
    __left = ee.rename(columns={"cluster1": "cluster", "cluster2": "old_node1"})
    # setnames(ee[setkeyv(x, c("cluster", "old_node1"))], "weight", "old_edge1")->x
    x = __left.merge(x, on=["cluster", "old_node1"], how="right").rename(columns={"weight": "old_edge1"})
    # setnames(ee, "old_node1", "old_node2")
    __left = ee.rename(columns={"cluster1": "cluster", "cluster2": "old_node2"})
    # setnames(ee[setkeyv(x, c("cluster", "old_node2"))], "weight", "old_edge2")->x
    x = __left.merge(x, on=["cluster", "old_node2"], how="right").rename(columns={"weight": "old_edge2"})
    # x[!is.na(new_edge3) ]->x
    # This is a hack to exclude NaN values.
    x = x.query("new_edge3 == new_edge3")
    t = edge_list[["cluster1", "cluster2"]].merge(df.copy(), left_on="cluster2", right_on="cluster", how="right")
    t.loc[:, "dist"] = t.groupby("cluster1")["bin"].transform(lambda s: 0 if s.shape[0] == 0 else
                                                              np.concatenate([s[1:].values - s[:-1].values, [0]]))
    idx = (t["dist"] == 1)
    try:
        t = pd.DataFrame().assign(
            cluster=t.loc[idx, "cluster"].values, new_node1=t.loc[idx, "cluster2"].values,
            new_node2=t.loc[idx.shift().fillna(False), "cluster2"].values)
    except KeyError as exc:
        print(t.head(10))
        print(t.columns)
        raise KeyError(exc)
    #    setkeyv(el[, .(cluster1, cluster2, weight)], c("cluster1", "cluster2"))->ee
    ee = edge_list[["cluster1", "cluster2", "weight"]]
    # setnames(ee, c("cluster1", "cluster2"), c("cluster", "new_node1"))
    __left = ee.rename(columns={"cluster1": "cluster", "cluster2": "new_node1"})
    # setnames(ee[setkeyv(t, c("cluster", "new_node1"))], "weight", "new_edge1")->t
    t = __left.merge(t, on=["cluster", "new_node1"], how="right").rename(columns={"weight": "new_edge1"})
    # setnames(ee, "new_node1", "new_node2")
    __left = ee.rename(columns={"cluster1": "cluster", "cluster2": "new_node2"})
    # setnames(ee[setkeyv(t, c("cluster", "new_node2"))], "weight", "new_edge2")->t
    t = __left.merge(t, on=["cluster", "new_node2"]).rename(columns={"weight": "new_edge2"})
    # setnames(ee, c("cluster", "new_node2"), c("new_node1", "new_node2"))
    __left = ee.rename(columns={"cluster1": "new_node1", "cluster2": "new_node2"})
    #    setnames(ee[setkeyv(t, c("new_node1", "new_node2"))], "weight", "old_edge3")->t
    t = __left.merge(t, on=["new_node1", "new_node2"]).rename(columns={"weight": "old_edge3"})
    #    setkey(x, "cluster")[setkey(t, "cluster")][!is.na(old_node1)]->m
    # Again the "x == x" trick to exclude NaN values
    m = x.merge(t, on=["cluster"]).query("old_node1 == old_node1")
    m.eval("diff = new_edge1 + new_edge2 + new_edge3 - old_edge1 - old_edge2 - old_edge3", inplace=True)
    m = m.query("diff < 0").sort_values("diff")
    return m


def node_relocation(df, edge_list, maxiter=100, verbose=True):
    """Traveling salesman heuristic for Hi-C mapping construction."""
    i = 0
    if df.shape[0] > 2:
        m = checker(df, edge_list)
        df = df.sort_values("bin")
        while m.shape[0] > 0 and i <= maxiter:
            i += 1
            ne = m.head(1)
            idx = pd.DataFrame().assign(cluster=df["cluster"],
                                        idx=np.arange(1, df.shape[0] + 1, dtype=int))
            for column in ["cluster", "new_node1", "new_node2", "old_node1", "old_node2"]:
                ne = pd.merge(ne, idx, left_on=column, right_on="cluster").rename(
                    columns={"idx": "idx_{column}".format(column=column)})

            # Now we have to do the reordering of the dataframe df
            # min_new <- min(idx_new_node1, idx_new_node2)
            # max_new <- max(idx_new_node1, idx_new_node2)
            # min_old <- min(idx_old_node1, idx_old_node2)
            # max_old <- max(idx_old_node1, idx_old_node2)
            min_new, max_new = np.sort([ne["idx_new_node1"], ne["idx_new_node2"]])
            min_old, max_old = np.sort([ne["idx_old_node1"], ne["idx_old_node2"]])
            if min_old < min_new:
                df = pd.concat([df.iloc[:min_old],
                                df.iloc[max_old:min_new],
                                df.iloc[ne["idx_cluster"]],
                                df.iloc[max_new:]])
            else:
                df = pd.concat([df.iloc[:min_new],
                                df.iloc[ne["idx_cluster"]],
                                df.iloc[max_new:min_old],
                                df.iloc[max_old:]])

            df.loc[:, "bin"] = np.arange(df.shape[0], dtype=int)
            if verbose:
                print("Node relocation steps:", i - 1)
            m = checker(df, edge_list)
            continue

    return df
