import pandas as pd
import numpy as np


def hic_cov_psmol(hic_map, binsize=1e3, binsize2=1e5, maxdist=1e6, cores=1):
    """calculate Hi-C physical coverage for pseudomolecules
 fpairs <- copy(hic_map$links)
 info <- hic_map$chrlen
 setnames(fpairs, c("start1", "start2"), c("pos1", "pos2"))

 fpairs[chr1 == chr2 & pos1 < pos2][, .(chr = chr1, bin1 = pos1 %/% binsize * binsize, bin2 =pos2 %/% binsize * binsize)]->f
 f[bin2 - bin1 > 2*binsize & bin2 - bin1 <= maxdist]->f
 f[, i := 1:.N]
 f[, b := paste0(chr, ":", bin1 %/% binsize2)]
 setkey(f, b)
 rbindlist(mclapply(mc.cores=cores, unique(f$b), function(j){
  f[j][, .(chr, bin=seq(bin1+binsize, bin2-binsize, binsize)), key=i][, .(n=.N), key=.(chr, bin)]
 }))->ff

 ff[, .(n=sum(n)), key=.(chr, bin)]->ff
 info[, .(chr, length)][ff, on="chr"]->ff
 ff[, d := pmin(bin, (length-bin) %/% binsize * binsize)]
 ff[, nbin := .N, key="chr"]
 ff[, mn := mean(n), key=d]
 ff[, r := log2(n/mn)][]
}
    """
    fpairs = hic_map["links"].rename(columns={"start1": "pos1", "start2": "pos2"})
    info = hic_map["chrlen"]
    # Extract those positions where the pairs match to the same chromosome? and pos1 < pos2
    bait = (fpairs.chr1 == fpairs.chr2) & (fpairs["pos1"] < fpairs["pos2"])
    f = pd.DataFrame(
        {"chr": fpairs.loc[bait]["chr1"],
         "bin1": (fpairs["pos1"] // binsize) * binsize,
         "bin2": (fpairs["pos2"] // binsize) * binsize
         }
    )
    f = f[(f["bin2"] - f["bin1"] > 2 * binsize) & (f["bin2"] - f["bin1"] <= maxdist)]
    f["i"] = list(range(1, f.shape[0] + 1))
    f["b"] = f["chr"].astype(str) + ":" + (f["bin1"] // binsize2).astype(str)
    # In the following section we have to achieve the following:
    # for each discrete bin, count how many internal bins there are.
    def aggregator(f, j):
        f[f.b == j].groupby("chr")

    rbindlist(mclapply(mc.cores = cores, unique(f$b), function(j)
    {
        # Extract those rows where b == j
        f[j]
        #
        [,.(chr, bin=seq(bin1 + binsize, bin2 - binsize, binsize)), key = i]
    [,.(n=.N), key =.(chr, bin)]
    }))->ff

    ff["n"]

    ff[,.(n=sum(n)), key =.(chr, bin)]->ff
    info[,.(chr, length)][ff, on = "chr"]->ff
    ff[, d := pmin(bin, (length - bin) % / % binsize * binsize)]
    ff[, nbin :=.N, key = "chr"]
    ff[, mn := mean(n), key = d]
    ff["mn"] = ff.groupby(["d"]).
    # ff[, r := log2(n / mn)][]
    ff["r"] = np.log2(ff["n"] / ff["mn"])
    return ff
