# Calculate N50 (and similar statistics) from a vector scaffold lengths
n50<-function(l, p=0.5){
 l[order(l)] -> l
 l[head(which(cumsum(as.numeric(l)) >= (1 - p) * sum(as.numeric(l))), n=1)]
}
