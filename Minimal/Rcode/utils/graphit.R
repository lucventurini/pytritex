library(data.table)
library(parallel)
library(igraph)
library(stringi)
library(zoo)


# Insert nodes not in MST backbone to where they fit best
insert_node<-function(df, el){
 while({
  setkey(el[cluster2 %in% df$cluster & cluster1 %in% unique(el[!cluster1 %in% df$cluster & cluster2 %in% df$cluster]$cluster1)], "cluster2")[setkey(setnames(df[, .(cluster, bin)], c("cluster2", "bin")), "cluster2"), allow.cartesian=T][!is.na(cluster1)][
   order(cluster1, bin)][, dist:={if(.N == 1) {as.integer(NA)} else { as.integer(c(bin[2:.N],NA)-bin)}}, by=cluster1]->y
   y[order(cluster1, bin)]->y
   length(which(y$dist == 1)->idx) > 0
  }) {

   setkeyv(setnames(el[, .(cluster1, cluster2, weight)], c("path_node1", "path_node2", "old_path")), c("path_node1", "path_node2"))[
    setkeyv(data.table(cluster=y[idx, cluster1], path_node1=y[idx, cluster2], path_node2=y[idx+1, cluster2], 
	      weight1=y[idx, weight], weight2=y[idx+1, weight], bin=y[idx, bin]), c("path_node1", "path_node2"))]->z
    head(z[ ,diff:= weight1 + weight2 - old_path][order(diff)],1)->z
   
   m=z$bin
   data.table(rbind(df[1:m], data.table(cluster=z$cluster, bin=m+1), df[(m+1):nrow(df),][, bin:=bin+1][]))->df
 }
 df
}

# Traveling salesman heuristic for Hi-C mapping construction
node_relocation<-function(df, el, maxiter=100, verbose=T){
 i<-0
 if(nrow(df) > 2){
  while({
   i<-i+1
   df[order(bin)]->df
   x<-data.table(old_node1=df[1:(nrow(df)-2)]$cluster, cluster=df[2:(nrow(df)-1)]$cluster, old_node2=df[3:nrow(df)]$cluster)
   setkeyv(el[, .(cluster1, cluster2, weight)], c("cluster1", "cluster2"))->ee
   setnames(setnames(ee, c("cluster1", "cluster2"), c("old_node1", "old_node2"))[setkeyv(x, c("old_node1", "old_node2"))], "weight", "new_edge3")->x
   setnames(setnames(ee, c("old_node1", "old_node2"), c("cluster", "old_node1"))[setkeyv(x, c("cluster", "old_node1"))], "weight", "old_edge1")->x
   setnames(setnames(ee, "old_node1", "old_node2")[setkeyv(x, c("cluster", "old_node2"))], "weight", "old_edge2")->x
   x[!is.na(new_edge3) ]->x

   setkey(el[, .(cluster1, cluster2)], "cluster2")[setkey(copy(df), "cluster")][order(cluster1, bin)]->t
   which(t[, dist:={if(.N == 1) {as.integer(NA)} else { as.integer(c(bin[2:.N],NA)-bin)}}, by=cluster1]$dist == 1)->idx

   data.table(cluster=t$cluster1[idx], new_node1=t$cluster2[idx], new_node2=t$cluster2[idx+1])->t
   setkeyv(el[, .(cluster1, cluster2, weight)], c("cluster1", "cluster2"))->ee
   setnames(setnames(ee, c("cluster1", "cluster2"), c("cluster", "new_node1"))[setkeyv(t, c("cluster", "new_node1"))], "weight", "new_edge1")->t
   setnames(setnames(ee, "new_node1", "new_node2")[setkeyv(t, c("cluster", "new_node2"))], "weight", "new_edge2")->t
   setnames(setnames(ee, c("cluster", "new_node2"), c("new_node1", "new_node2"))[setkeyv(t, c("new_node1", "new_node2"))], "weight", "old_edge3")->t
   setkey(x, "cluster")[setkey(t, "cluster")][!is.na(old_node1)]->m
   m[ ,diff := new_edge1+new_edge2+new_edge3 - old_edge1 - old_edge2 - old_edge3]
   nrow(m<-m[diff < 0][order(diff)]) > 0 & i <= maxiter
   }) {
   ne<-head(data.frame(m), 1)
   idx<-data.frame(cluster=df$cluster, idx=1:nrow(df))
   invisible(lapply(c("cluster", "new_node1", "new_node2", "old_node1", "old_node2"), function(i) {
    merge(ne, by.x=i, by.y="cluster", idx)->>ne
    colnames(ne)[which(colnames(ne) == "idx")]<<-paste(sep="_", "idx", i)
   }))
   df[with(ne, {
    min_new <- min(idx_new_node1, idx_new_node2)
    max_new <- max(idx_new_node1, idx_new_node2)
    min_old <- min(idx_old_node1, idx_old_node2)
    max_old <- max(idx_old_node1, idx_old_node2)
    if(min_old < min_new) {
     c(1:min_old, max_old:min_new, idx_cluster, max_new:nrow(df))
    } else {
     c(1:min_new, idx_cluster, max_new:min_old, max_old:nrow(df))
    }
   }),]->df
   df$bin<-1:nrow(df)
   df<-data.table(df)
  }
  if(verbose){
   cat(paste0("Node relocation steps: ", i-1, "\n"))
  }
 }
 df
}

# Traveling salesman heuristic for Hi-C mapping construction
kopt2<-function(df, el){
 setkeyv(el, c("cluster1", "cluster2"))->el
 while({
  df[, .(cluster1=cluster[1:(nrow(df)-1)], cluster2=cluster[2:nrow(df)])]->m
  el[setkeyv(m, c("cluster1", "cluster2"))]->m
  m[ ,.(cluster1, cluster2, weight)]->m
  setnames(m, "weight", "weight12")->m
  m[, dummy:=1]
  setkey(setnames(copy(df), c("cluster", "bin"), c("cluster1", "bin1")), "cluster1")[setkey(m, "cluster1")]->m
  setkey(setnames(copy(df), c("cluster", "bin"), c("cluster2", "bin2")), "cluster2")[setkey(m, "cluster2")]->m
  copy(m)->n
  setnames(n, c("cluster1", "cluster2", "bin1", "bin2", "weight12"), c("cluster3", "cluster4", "bin3", "bin4", "weight34"))
  setkey(m, "dummy")[setkey(n, "dummy"), allow.cartesian=T]->mn
  mn[, dummy:=NULL]
  mn[bin1 < bin3]->mn
  o<-el[, .(cluster1, cluster2, weight)]
  setkeyv(setnames(copy(o), c("cluster1", "cluster3", "weight13")), c("cluster1", "cluster3"))[setkeyv(mn, c("cluster1", "cluster3"))]->mn
  setkeyv(setnames(copy(o), c("cluster2", "cluster4", "weight24")), c("cluster2", "cluster4"))[setkeyv(mn, c("cluster2", "cluster4"))]->mn
  mn[, old:=weight12+weight34]
  mn[, new:=weight13+weight24]
  mn[, diff:=old-new]
  mn[order(-diff)]->mn
  nrow(mn[diff > 0]) > 0 
 }) {
  head(mn, 1)->x
  bin1<-df[cluster == x$cluster1]$bin
  bin2<-df[cluster == x$cluster2]$bin
  bin3<-df[cluster == x$cluster3]$bin
  bin4<-df[cluster == x$cluster4]$bin
  df[c(1:bin1, bin3:bin2, bin4:nrow(df))]->df
  df[, bin:=1:nrow(df)]->df
 }
 df
}