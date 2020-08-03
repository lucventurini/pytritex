# Define mapping table between numeric and proper chromosome names (21 -> 7D in Triticum aestivum)
chrNames<-function(agp=F, species="wheat") {
 if(species == "wheat"){
  data.table(alphachr=apply(expand.grid(1:7, c("A", "B", "D"), stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:21)->z
 } else if (species == "barley"){
  data.table(alphachr=apply(expand.grid(1:7, "H", stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:7)->z
 }
 else if (species == "rye"){
  data.table(alphachr=apply(expand.grid(1:7, "R", stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:7)->z
 }
 else if (species == "lolium"){
  data.table(alphachr=apply(expand.grid(1:7, "", stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:7)->z
 }
 else if (species == "sharonensis"){
  data.table(alphachr=apply(expand.grid(1:7, "S", stringsAsFactors=F), 1, function(x) paste(x, collapse="")), chr=1:7)->z
 }
 else if (species == "oats"){
  data.table(alphachr=sub(" ", "", apply(expand.grid(1:21, "M", stringsAsFactors=F), 1, function(x) paste(x, collapse=""))), chr=1:21)->z
 }
 if(agp){
  rbind(z, data.table(alphachr="Un", chr=0))[, agp_chr := paste0("chr", alphachr)]->z
 }
 z[]
}

# alias for backwards compatibility
wheatchr <- chrNames