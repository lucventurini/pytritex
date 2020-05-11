import pandas as pd
import matplotlib.pyplot as plot
import numpy as np
import tempfile
import subprocess as sp
import multiprocessing as mp


def _plotter(index, row, height, res, width):
    # file <- bad[i, f]
    #   pdf <- bad[i, p]
    #   s<-bad[i,scaffold]
    #   cat(paste0(i, " ", bad[i, scaffold], "\n"))
    #   png(file, height=height, res=res, width=width)
    #   plot_popseq_carma_tcc(s, breaks=breaks, page=i, info, cssaln, tcc_pos, ff[d >= mindist], molcov, species)
    #   dev.off()
    #   system(paste("convert", file, pdf))
    #   unlink(file)
    pdf = row["pdf"]
    scaffold = row["scaffold"]
    figure = plot.figure(figsize=(height, width), dpi=res)
    plot_popseq_carma_tcc(scaffold, breaks=breaks, page=index + 1, info=info, cssaln=cssaln, tcc_pos=tcc_pos,
                          ?=ff[ff["d"] >= mindist], ?=molcov, ?=species)
    figure.savefig(format="pdf")


def plot_chimeras(scaffolds, assembly, filename: str, breaks=None, mindist=0, cores=20,
                  species="wheat", mbscale=False, autobreaks=True):
    """

 bad <- copy(scaffolds)
 out <- file
 bad[, f:=tempfile(fileext=".png"), by=scaffold]
 bad[, p:=tempfile(fileext=".pdf"), by=scaffold]
 mclapply(mc.cores=cores, 1:nrow(bad), function(i){

 })
 system(paste0("gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=", out, " ", paste(bad$p, collapse=' ')))
 unlink(bad$p)
    """
    info = assembly["info"]
    ff = assembly["cov"]
    tcc_pos = assembly["fpairs"]
    cssaln = assembly["cssaln"]
    molcov = assembly["molecule_cov"]
    width, height, res = 2500, 1000, 150
    if tcc_pos is not None and tcc_pos.shape[0] > 0:
        height += 1000
    if molcov is not None and molcov.shape[0] > 0:
        height += 1000
    scaffolds = scaffolds.loc[~scaffolds.scaffold.duplicated(), "scaffold"]
    bad = scaffolds.copy()
    out = filename
    for group, item in bad.groupby("scaffold").groups.items():
        bad.loc[item, "png"] = tempfile.NamedTemporaryFile(suffix=".png").name
    pool = mp.Pool(cores)
    pool.starmap_async(

    )

    sp.call("gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile={out} {pdfs}".format(
        out=out, pdfs=" ".join(bad["pdf"].astype(str))),
        shell=False)
