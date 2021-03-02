insert_node<-function(df, el) {
 y_calculator <- function(df, el){
   bait <- unique(el[!cluster1 %in% df$cluster & cluster2 %in% df$cluster]$cluster1);
   left <- el[cluster2 %in% df$cluster & cluster1 %in% bait]
   setkey(left, "cluster2");
   right <- df[, .(cluster, bin)]
   setnames(right, c("cluster2", "bin"))
   setkey(right, "cluster2")
   merged <- left[right, allow.cartesian=T][!is.na(cluster1)]
   merged <- merged[order(cluster1, bin)]
   y <- merged[, dist:={if(.N == 1) {as.integer(NA)} else { as.integer(c(bin[2:.N],NA)-bin)}}, by=cluster1]
   y <- y[order(cluster1, bin)]
   return(y)
  }

  y <- y_calculator(df, el)
  idx <- which(y$dist == 1)
  while (length(idx) > 0){
     left <- el[, .(cluster1, cluster2, weight)]
     setnames(left, c("path_node1", "path_node2", "old_path"))
     setkeyv(left, c("path_node1", "path_node2"))
     right <- data.table(cluster=y[idx, cluster1],
                         path_node1=y[idx, cluster2],
                         path_node2=y[idx+1, cluster2],
                         weight1=y[idx, weight],
                         weight2=y[idx+1, weight], bin=y[idx, bin])
     setkeyv(right, c("path_node1", "path_node2"))
     z <- left[right]
     z <- z[ ,diff:= weight1 + weight2 - old_path][order(diff)]
     z <- head(z,1)
     m=z$bin  # This is the index
     df <- data.table(
       rbind(df[1:m],
             data.table(cluster=z$cluster, bin=m+1),
             df[(m+1):nrow(df),][, bin:=bin+1][])
     )
     y <- y_calculator(df, el)
     idx <- which(y$dist == 1)
 }
 df
}