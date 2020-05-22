insert_node<-function(df, el) {
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