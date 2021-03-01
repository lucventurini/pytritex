insert_node<-function(df, el) {
 # Inputs: DataFrame DF with "cluster", "bin" -> the path.
 # Inputs: Edge list of the form "cluster1", "cluster2", "weight"

 y_calculator <- function(df, el){
  # This function will find out whether there are nodes to insert in the backbone

  # First check: the *left* link must *not* be in the path
  # The right bit *must* be in the path. Select the left
  internal_bait <- unique(el[!cluster1 %in% df$cluster & cluster2 %in% df$cluster]$cluster1);
  # Now select edges where the *right* bit is the path and the left *is not*.
  left <- el[cluster2 %in% df$cluster & cluster1 %in% internal_bait]
  setkey(left, "cluster2")
  # Now merge with the dataframe. Use the *right* bit as a bait.
  right <- setnames(df[, .(cluster, bin)], c("cluster2", "bin"))
  setkey(right, "cluster2")
  # Do a cartesian product.
  merged <- left[right, allow.cartesian=T][!is.na(cluster1)]
  # Group by cluster, sort by bin.
  merged[order(cluster1, bin)]
  #
  y <- merged[, dist:={if(.N == 1) {as.integer(NA)} else { as.integer(c(bin[2:.N],NA)-bin)}}, by=cluster1]
  y[order(cluster1, bin)]
  return(y)
 }

 idxer <- function(y) {
  idx <- which(y$dist == 1)
  return(idx)
 }

 y <- y_calculator(df, el)
 idx <- idxer(y)
 while(length(idx) > 0) {
   left <- el[, .(cluster1, cluster2, weight)]
   setnames(left, c("path_node1", "path_node2", "old_path"))
   setkeyv(left, c("path_node1", "path_node2"))
   right <- data.table(cluster=y[idx, cluster1],
                       path_node1=y[idx, cluster2],
                       path_node2=y[idx+1, cluster2],
	                   weight1=y[idx, weight],
                       weight2=y[idx+1, weight],
                       bin=y[idx, bin])
   setkeyv(right, c("path_node1", "path_node2"))
   z <- left[right]
   z <- head(z[ ,diff:= weight1 + weight2 - old_path][order(diff)], 1)
   m <- z$bin
   df <- data.table(rbind(
   df[1:m],
   data.table(cluster=z$cluster, bin=m+1),
   df[(m+1):nrow(df),][, bin:=bin+1][]))
   y <- y_calculator(df, el)
   idx <- idxer(y)
 }
 df
}


make_super_path<-function(super, idx=NULL, start=NULL, end=NULL, maxiter=100, verbose=T){

 # Get backbone from minimum spanning tree (MST)
 # Take the those clusters that are in the super_scaffold marked as "idx"
 submem<-super$mem[super == idx]
 # Get the edges AND the weights, put in "el"
 super$edges[super == idx, .(cluster1, cluster2, weight)]->el
 # Get the minimum spanning tree
 minimum.spanning.tree(induced.subgraph(super$graph, submem$cluster))->mst

 # Set the weights to one, derive the shortest path between start and end
 # IF start and end are not known, calculate the diameter (ie the longest shortest path between two vertices)
 # and relative path.
  E(mst)$weight <- 1

 if(is.null(start) | is.null(end)){
  V(mst)[get.diameter(mst)]$name->dia
 } else {
  V(mst)[get.shortest.paths(mst, from=start, to=end)$vpath[[1]]]$name->dia
 }

 #Create a dataframe "df" of the form (in parallel, here *transposed*):
 # node1, node2, node3 .... nodeN
 # 1, 2, 3, ... N
 data.table(cluster=dia, bin=1:length(dia))->df

 # Insert the nodes absent from the derived path.
 df<-insert_node(df, el)
 # Traveling salesman heuristics
 df<-kopt2(df, el)
 # Perform the node relocation heuristic
 df<-node_relocation(df, el, maxiter=maxiter, verbose=verbose)

 # Block optimise heuristic
 data.frame(df)->df
 data.frame(el)->el
 data.frame(cluster=df$cluster, rank = 0)->ranks
 r <- 0

 n <- unique(subset(el, !cluster1 %in% df$cluster & cluster2 %in% df$cluster)$cluster1);

 while(length(n > 0)) {
  r <- r+1
  subset(el, cluster2 %in% df$cluster & cluster1 %in% n)->tmp
  tmp[!duplicated(tmp$cluster1),]->tmp
  rbind(ranks, data.frame(cluster=tmp$cluster1, rank=r))->ranks
  merge(tmp, df[c("cluster", "bin")], by.x="cluster2", by.y="cluster")->x
  rbind(df, data.frame(cluster=x$cluster1, bin=x$bin))->df
  n <- unique(subset(el, !cluster1 %in% df$cluster & cluster2 %in% df$cluster)$cluster1);
 }

 merge(df, submem)->df
 df$bin<-as.integer(df$bin)
 flip<-with(df, suppressWarnings(cor(bin, cM))) < 0
 if((!is.na(flip)) & flip) {
  with(df, max(bin) - bin + 1)->df$bin
 }
 ranks$rank<-as.numeric(ranks$rank)
 merge(df, ranks)->df
 df$backbone <- df$cluster %in% dia

 data.table(df)[, .(cluster, bin, rank, backbone)]
}