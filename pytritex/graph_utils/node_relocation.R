# Title     : TODO
# Objective : TODO
# Traveling salesman heuristic for Hi-C mapping construction

library(data.table)

checker <- function(df, el, i){
   i<-i+1
   df <- df[order(bin)]
   x<-data.table(old_node1=df[1:(nrow(df)-2)]$cluster,
                 cluster=df[2:(nrow(df)-1)]$cluster,
                 old_node2=df[3:nrow(df)]$cluster)
   # Take up the edge list ..
   ee <- setkeyv(el[, .(cluster1, cluster2, weight)], c("cluster1", "cluster2"))
   # Now we are getting the weight of the edge between the two outer nodes
   # And calling it "new_edge3"
   left <- setnames(ee, c("cluster1", "cluster2"), c("old_node1", "old_node2"))
   right <- setkeyv(x, c("old_node1", "old_node2"))
   x <- setnames(left[right], "weight", "new_edge3")
   # Here we are getting the cost of going from the leftmost to the one in the middle
   left <- setnames(ee, c("old_node1", "old_node2"), c("cluster", "old_node1"))
   right <- setkeyv(x, c("cluster", "old_node1"))
   x <- setnames(left[right], "weight", "old_edge1")
   # Here we are getting the cost between the mid node and the right most
   left <- setnames(ee, "old_node1", "old_node2")
   right <- setkeyv(x, c("cluster", "old_node2"))
   x <- setnames(left[right], "weight", "old_edge2")
   # Select only places where the leftmost and rightmost node ARE connected
   x[!is.na(new_edge3) ]->x
   # Now order the edges according to the righmost node
   left <- setkey(el[, .(cluster1, cluster2)], "cluster2")
   right <- setkey(copy(df), "cluster")
   # and merge with the current path
   t <- left[right][order(cluster1, bin)]
   # Now get the difference between the possible edge costs for each node in the path
   dist <- t[, dist:={if(.N == 1) {as.integer(NA)} else { as.integer(c(bin[2:.N],NA)-bin)}}, by=cluster1]
   # Now get those positions where the distance is minimal
   idx <- which(dist$dist == 1)
   # Now create a new data table composed of clusters where the distance is 1, the subsequent node, and the node after that
   t <- data.table(cluster=t$cluster1[idx], new_node1=t$cluster2[idx], new_node2=t$cluster2[idx+1])
   # Go back to the edge list ...
   ee <- setkeyv(el[, .(cluster1, cluster2, weight)], c("cluster1", "cluster2"))
   # and get the weight of all these potential new edges ..
   left <- setnames(ee, c("cluster1", "cluster2"), c("cluster", "new_node1"))
   right <- setkeyv(t, c("cluster", "new_node1"))
   t <- setnames(left[right], "weight", "new_edge1")
   left <- setnames(ee, "new_node1", "new_node2")
   right <- setkeyv(t, c("cluster", "new_node2"))
   t <- setnames(left[right], "weight", "new_edge2")
   left <- setnames(ee, c("cluster", "new_node2"), c("new_node1", "new_node2"))
   right <- setkeyv(t, c("new_node1", "new_node2"))
   t <- setnames(left[right], "weight", "old_edge3")
   m <- setkey(x, "cluster")[setkey(t, "cluster")][!is.na(old_node1)]
   # And now check the difference of the old path with the new path.
   m[ ,diff := new_edge1 + new_edge2 + new_edge3 - old_edge1 - old_edge2 - old_edge3]
  }


node_relocation<-function(df, el, maxiter=100, verbose=T){
  i <- 0

 if(nrow(df) > 2){
  while({
    i <- i + 1
    # While there is at least one node that, if relocated, would yield a decrease in weight in the path ...
    m <- checker(df, el, i)
    nrow(m<-m[diff < 0][order(diff)]) > 0 & i <= maxiter
  }) {
    # Select the best relocated node ...
    ne<-head(data.frame(m), 1)
    idx<-data.frame(cluster=df$cluster, idx=1:nrow(df))
   # Get the other columns for this row
   invisible(
     lapply(
       c("cluster", "new_node1", "new_node2", "old_node1", "old_node2"),
       function(i) {
             merge(ne, by.x=i, by.y="cluster", idx)->>ne
             colnames(ne)[which(colnames(ne) == "idx")]<<-paste(sep="_", "idx", i)
   }))
   # Get the indices to reorder / subset the data frame
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
