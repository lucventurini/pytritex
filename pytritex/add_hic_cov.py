import pandas as pd
import numpy as np
import multiprocessing as mp

# Calculate physical coverage with Hi-C links in sliding windows along the scaffolds

def add_hic_cov(assembly, scaffolds=None, binsize=1e3, binsize2=1e5, minNbin=50, innerDist=1e5, cores=1):
    """
 fpairs[scaffold1 == scaffold2 & pos1 < pos2][, .(scaffold = scaffold1, bin1 = pos1 %/% binsize * binsize, bin2 =pos2 %/% binsize * binsize)]->f
 f[bin2 - bin1 > 2*binsize]->f
 f[, i := 1:.N]
 f[, b := paste0(scaffold, ":", bin1 %/% binsize2)]
 setkey(f, b)

 rbindlist(mclapply(mc.cores=cores, unique(f$b), function(j){
  f[j][, .(scaffold=scaffold, bin=seq(bin1+binsize, bin2-binsize, binsize)), key=i][, .(n=.N), key=.(scaffold, bin)]
 }))->ff

 if(nrow(ff) > 0){
  ff[, .(n=sum(n)), key=.(scaffold, bin)]->ff
  info[, .(scaffold, length)][ff, on="scaffold"]->ff
  ff[, d := pmin(bin, (length-bin) %/% binsize * binsize)]
  ff[, nbin := .N, key="scaffold"]
  ff[, mn := mean(n), key=d]
  ff[, r := log2(n/mn)]
  ff[nbin >= minNbin, .(mr=suppressWarnings(min(r))), key=scaffold][order(mr)]->z
  ff[nbin > minNbin & d >= innerDist, .(mri=suppressWarnings(min(r))), key=scaffold]->zi
  z[ff, on="scaffold"]->ff
  zi[ff, on="scaffold"]->ff
  z[info, on="scaffold"]->info_mr
  zi[info_mr, on="scaffold"]->info_mr
 } else {
  copy(info) -> info_mr
  info_mr[, c("mri", "mr") := list(NA, NA)]
 }

 if(null){
  assembly$info=info_mr
  assembly$cov=ff

  assembly$binsize <- binsize
  assembly$minNbin <- minNbin
  assembly$innerDist <- innerDist

  assembly
 } else {
  list(info=info_mr, cov=ff)
 }
}

    :param assembly:
    :param scaffolds:
    :param binsize:
    :param binsize2:
    :param minNbin:
    :param innerDist:
    :param cores:
    :return:
    """

    info = assembly["info"]
    if "mr" in info.columns or "mri" in info.columns:
        raise KeyError("Assembly[info] already has mr and/or mri columns; aborting.")
    fpairs = assembly["fpairs"]
    if scaffolds is None:
        scaffolds = info["scaffolds"]
        null = True
    else:
        info = info[info.scaffold.isin(scaffolds)]
        fpairs = fpairs[fpairs.scaffold1.isin(scaffolds)]
        null = False

    bait = (fpairs.scaffold1 == fpairs.scaffold2) & (fpairs.pos1 < fpairs.pos2)

    f = pd.DataFrame(
        {"scaffold": fpairs.loc[bait]["scaffold1"],
         "bin1": (fpairs.loc[bait]["pos1"] // binsize) * binsize,
         "bin2": (fpairs.loc[bait]["pos2"] // binsize) * binsize,
         }
    )
    f = f[f.bin2 - f.bin1 > 2 * binsize]
    f["i"] = range(1, f.shape[0] + 1)  # Not sure about this yet
    f["b"] = f["scaffold"].astype(str) + ":" + (f["bin1"] // binsize2).astype(str)
    grouped = f.groupby("b")
    pool = mp.Pool(cores)

    def aggregator(group, binsize):
        # {
        #     f[j][,.(scaffold=scaffold, bin=seq(bin1 + binsize, bin2 - binsize, binsize)), key = i][,.(n=.N), key =.(
        # scaffold, bin)]
        # }))->ff
        return pd.concat(
            [pd.DataFrame({"bin": list(range(row['bin1'] + binsize, row['bin2'] - binsize, binsize)),
                           "scaffold": row.scaffold, "i": row.i}) for
             index, row in group.iterrows()]).groupby(
            ["scaffold", "bin"])["scaffold"].count().to_frame("n").reset_index()

    results = [pool.apply_async(aggregator, (grouped.get_group(group), binsize))
               for group in grouped.groups.keys()]
    ff = pd.concat([res.get() for res in results])

    if ff.shape[0] > 0:
        # Sum again as we had clustered by "b"
        ff = ff.groupby(["scaffold", "bin"]).n.sum().to_frame("n")
        ff = pd.merge(info[["scaffold", "length"]], ff, left_on="scaffold", right_on="scaffold")
        ff["d"] = ff.apply(lambda row: np.min(row.bin, ((row.length - row.bin) // binsize) * binsize))
        ff = pd.merge(ff, ff.groupby("scaffold")["scaffold"].count().to_frame("nbin"),
                      left_on=["scaffold"], right_index=True)
        ff = pd.merge(ff, ff[["d", "n"]].groupby("d")["n"].mean().to_frame("mn"),
                      left_on=["d"], right_index=True)
        ff["r"] = ff.apply(lambda row: np.log2(row.n / row.mn), axis=1)
        bait = (ff.nbin >= minNbin)
        z = ff.loc[bait, ["scaffold", "r"]].groupby("scaffold")["r"].min().to_frame("mr").sort_values("mr")
        bait = (ff.nbin > minNbin) & (ff.d >= innerDist)
        zi = ff.loc[bait, ["scaffold", "r"]].groupby("scaffold")["r"].min().to_frame("mri")
        ff = pd.merge(pd.merge(ff, z, left_on=["scaffold"], right_index=True),
                      zi, left_on=["scaffold"], right_index=True)
        info_mr = pd.merge(pd.merge(info, z, left_on="scaffold", right_index=True),
                           zi, left_on="scaffold", right_index=True)
    else:
        info_mr = info.copy()
        info_mr["mri"], info["mr"] = np.nan, np.nan

    return info_mr
