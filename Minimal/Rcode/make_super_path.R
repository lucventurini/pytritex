library(data.table)
library(parallel)
library(igraph)
library(stringi)
library(zoo)


make_super_path<-function(super, idx=NULL, start=NULL, end=NULL, maxiter=100, verbose=T){

 # Get backbone from minimum spanning tree (MST)
 submem<-super$mem[super == idx]
 super$edges[super == idx, .(cluster1, cluster2, weight)]->el
 minimum.spanning.tree(induced.subgraph(super$graph, submem$cluster))->mst

 E(mst)$weight <- 1

 if(is.null(start) | is.null(end)){
  V(mst)[get.diameter(mst)]$name->dia
 } else {
  V(mst)[get.shortest.paths(mst, from=start, to=end)$vpath[[1]]]$name->dia
 }

data.table(cluster=dia, bin=1:length(dia))->df

 df<-insert_node(df, el)
 # Traveling salesman heuristics
 df<-kopt2(df, el)
 df<-node_relocation(df, el, maxiter=maxiter, verbose=verbose)

 data.frame(df)->df
 data.frame(el)->el
 data.frame(cluster=df$cluster, rank = 0)->ranks
 r=0

 while(length(n<-unique(subset(el, !cluster1 %in% df$cluster & cluster2 %in% df$cluster)$cluster1)) > 0) {
  r = r+1
  subset(el, cluster2 %in% df$cluster & cluster1 %in% n)->tmp
  tmp[!duplicated(tmp$cluster1),]->tmp
  rbind(ranks, data.frame(cluster=tmp$cluster1, rank=r))->ranks
  merge(tmp, df[c("cluster", "bin")], by.x="cluster2", by.y="cluster")->x
  rbind(df, data.frame(cluster=x$cluster1, bin=x$bin))->df
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