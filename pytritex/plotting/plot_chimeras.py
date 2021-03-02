import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import tempfile
import subprocess as sp
import multiprocessing as mp
from .plot_popseq_carma_tcc import plot_popseq_carma_tcc


def _plotter(index, row, height, res, width, breaks, info, cssaln, tcc_pos, ff, species, molcov, mindist, autobreaks):
    pdf = row["pdf"]
    scaffold = row["scaffold"]
    figure = plt.figure(figsize=(height, width), dpi=res)
    plot_popseq_carma_tcc(figure=figure,
                          scaffold=scaffold,
                          breaks=breaks, page=index + 1, info=info, cssaln=cssaln, tcc_pos=tcc_pos,
                          span=ff[ff["d"] >= mindist], molcov=molcov, species=species, autobreaks=autobreaks)
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
