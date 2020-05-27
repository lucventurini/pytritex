 if(!raw){
   if(verbose){
    cat("Tip removal.\n")
   }
   # remove short tips of rank 1
   links[!scaffold2 %in% ex][, .(degree=.N), key=.(scaffold=scaffold1)]->b
   a <- b[
     m[m[rank == 1][, .(super=c(super, super, super), bin=c(bin, bin-1, bin+1))], on=c("super", "bin")],
          on="scaffold_index"]
   a[degree == 1 & length <= 1e4]$scaffold->add
   if(length(add) > 0){
    ex <- c(ex, add)
    make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
    out$m -> m
   }

   # remove short tips/bulges of rank 1
   m[rank == 1 & length <= 1e4]$scaffold -> add
   ex <- c(ex, add)
   if(length(add) > 0){
    ex <- c(ex, add)
    make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
    out$m -> m
   }

   # resolve length-one-bifurcations at the ends of paths
   m[rank > 0][bin == 2 | super_nbin - 1 == bin ][, .(super, super_nbin, type = bin == 2, scaffold, length, bin0=bin)]->x
   unique(rbind(
   m[x[type == T, .(super, bin0, bin=1)], on=c("super", "bin")],
   m[x[type == F, .(super, bin0, bin=super_nbin)], on=c("super", "bin")]
   ))->a
   a[, .(super, bin0, scaffold2=scaffold, length2=length)][x, on=c("super", "bin0")][, ex := ifelse(length >= length2, scaffold2, scaffold)]$ex -> add

   if(length(add) > 0){
    ex <- c(ex, add)
    make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
    out$membership -> m
   }

   # remove short tips/bulges of rank 1
   m[rank == 1 & length <= 1e4]$scaffold -> add
   if(length(add) > 0){
    ex <- c(ex, add)
    make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
    out$m -> m
   }

   # remove tips of rank 1
   links[!scaffold2 %in% ex][, .(degree=.N), key=.(scaffold=scaffold1)]->b
   b[m[rank == 1], on="scaffold_index"][degree == 1]$scaffold -> add
   if(length(add) > 0){
    ex <- c(ex, add)
    make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
    out$m -> m
   }

   # remove remaining nodes of rank > 0
   m[rank > 0]$scaffold -> add
   if(length(add) > 0){
    ex <- c(ex, add)
    make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
    out$m -> m
   }

   if(popseq_dist > 0 & unanchored == T){
    if(verbose){
     cat("Including unanchored scaffolds.\n")
    }
    # use unanchored scaffolds to link super-scaffolds
    ww2[is.na(popseq_chr1), .(scaffold_link=scaffold1, link_length=length1, scaffold1=scaffold2)]->x
    ww2[is.na(popseq_chr1), .(scaffold_link=scaffold1, scaffold2=scaffold2)]->y
    x[y, on="scaffold_link", allow.cartesian=T][scaffold1 != scaffold2]->xy

    m[, .(scaffold1=scaffold, super1=super, chr1=chr, cM1=cM, size1=super_nbin, d1 = pmin(bin - 1, super_nbin - bin))][xy, on="scaffold_index1"]->xy
    m[, .(scaffold2=scaffold, super2=super, chr2=chr, cM2=cM, size2=super_nbin, d2 = pmin(bin - 1, super_nbin - bin))][xy, on="scaffold_index2"]->xy
    xy[super2 != super1 & d1 == 0 & d2 == 0 & size1 > 1 & size2 > 1 & chr1 == chr2]->xy
    xy[scaffold1 < scaffold2, .(nscl=.N), scaffold_link][xy, on="scaffold_link"]->xy
    xy[nscl == 1] -> xy
    xy[super1 < super2][, c("n", "g"):=list(.N, .GRP), by=.(super1, super2)][order(-link_length)][!duplicated(g)]->zz

    sel <- zz[, .(scaffold1=c(scaffold_link, scaffold_link, scaffold1, scaffold2),
 	  scaffold2=c(scaffold1, scaffold2, scaffold_link, scaffold_link))]
    rbind(links, ww2[sel, on=c("scaffold_index1", "scaffold_index2")])->links2

    make_super_scaffolds(links=links2, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
    out$m -> m

    #resolve branches
    m[rank > 0][bin == 2 | super_nbin - 1 == bin ][, .(super, super_nbin, type = bin == 2, scaffold, length, bin0=bin)]->x
    unique(rbind(
    m[x[type == T, .(super, bin0, bin=1)], on=c("super", "bin")],
    m[x[type == F, .(super, bin0, bin=super_nbin)], on=c("super", "bin")]
    ))->a
    a[, .(super, bin0, scaffold2=scaffold, length2=length)][x, on=c("super", "bin0")][, ex := ifelse(length >= length2, scaffold2, scaffold)]$ex -> add

    if(length(add) > 0){
     ex <- c(ex, add)
     make_super_scaffolds(links=links2, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
     out$membership -> m
    }

    # remove remaining nodes of rank > 0
    m[rank > 0]$scaffold -> add
    if(length(add) > 0){
     ex <- c(ex, add)
     make_super_scaffolds(links=links, prefix=prefix, info=info, excluded=ex, ncores=ncores) -> out
    }
   }

   out$m -> m
   out$info -> res
  }
