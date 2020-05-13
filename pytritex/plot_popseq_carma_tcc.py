import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from .utils.chrnames import chrNames


def plot_popseq_carma_tcc(figure: plt.Figure,
                          scaffold: str, page: int,
                          info: pd.DataFrame, cssaln: pd.DataFrame,
                          tcc_pos=None,
                          span=None,
                          autobreaks=True,
                          molcov=None, species="wheat", breaks=None):

    wheatchr = chrNames(species=species)
    i = scaffold
    hic = not (tcc_pos is None or tcc_pos.shape[0] == 0)
    tenex = not (molcov is None or molcov.shape[0] == 0)
    nrow = 1 + int(hic) + int(tenex)
    ax, subplots = figure.subplots(nrow, 3)
    if hic is True:
        plt.subplots_adjust(left=1, bottom=4, right=4, top=1, wspace=0, hspace=0)
    else:
        plt.subplots_adjust(left=4, bottom=4, right=4, top=1, wspace=0, hspace=0)
    # par(cex=1) ???

    if breaks is not None:
        br = breaks.loc[breaks["scaffold"] == i, "br"] / 1e6
    elif hic and autobreaks:
        br = span.loc[span["scaffold"] == i, [["r", "bin"]]].sort_values("r").iloc[0]["bin"] / 1e6
    else:
        br = None

    if hic is True:
        ## TODO: ??
        #    info[scaffold == i, paste0(", bad Hi-C: ", bad_hic)] -> badhic
        badhic = info.loc[info["scaffold"] == i, ]
        xlab = ""
    else:
        badhic = ""
        xlab = "Position in scaffold (Mb)"

    z = cssaln.merge(scaffold, on="scaffold")
    l = info.loc[i, "length"] / 1e6

    ymin = 8 if species in ("rye", "lolium", "barley") else 24
    ylab = "PopSeq chromosome" if mbscale is False else "Morex V2 chromosome"

    # TODO
    #   plot(xlim=c(0,l), 0, ylim=c(ymin,1), bty='l', col=0, yaxt='n', pch=20, xlab=xlab, ylab=ylab)
    plt.plot()

    if br is not None:
        plt.vlines(br, colors="blue", linestyles="dashed", linewidths=3)

    if mbscale is False:
        plt.title("POPSEQ")
    else:
        plt.title("Morex V2")

    #   z[, points(pos/1e6, popseq_chr, pch=20)]
    z.assign(pos=z["pos"] / 1e6).plot.scatter("pos", "popseq_chr", 20)
    # TODO: ??
    #   axis(2, las=1, wheatchr$chr, wheatchr$alphachr)

    # TODO: ????
    #   if(mbscale == F){
    #    info[scaffold == i, title(outer=T, line=0, paste("page: ", page, ", ", scaffold, sep="", " (", round(l,1), " Mb),\n",
    # 			   "bad POPSEQ: ", bad_popseq, ", bad CARMA: ", bad_sorted, badhic))]
    #   } else {
    #    info[scaffold == i, title(outer=T, line=0, paste("page: ", page, "\n", scaffold, sep="", " (", round(l,1), " Mb)"))]
    #   }

    w = info.loc[info["scaffold"] == i, "popseq_chr"]
    if w.shape[0] == 0:
        plt.scatter()

#   w<-info[i, on="scaffold"]$popseq_chr
#   if(is.na(w)){
#    plot(xlab='', ylab='', type='n', 0, axes=F)
#   } else {
#    if(mbscale == F){
#     ylab = "genetic position (cM)"
#    } else {
#     ylab = "pseudomolecule position (Mb)"
#    }
#    z[popseq_chr == w][, plot(pch=20, xlab=xlab, ylab=ylab, xlim=c(0,l), bty='l', las=1, pos/1e6, popseq_cM)]
#    if(mbscale == F){
#     title(main="Genetic positions of CSS contigs\nfrom major chromosome")
#    } else {
#     title(main="Physical positions of 1 kb single-copy tags\nfrom major chromosome")
#    }
#    if(hic){
#     abline(v=br, lwd=3, col="blue")
#    }
#   }
#
#   if(species != 'lolium'){
#    if(mbscale == F){
#     ylab="flow-sorting chromosome"
#     ttt="flow-sorting CARMA"
#    } else {
#     ylab="chromsome arm"
#     ttt="Hi-C based chromosome arm assignment"
#    }
#    plot(xlim=c(0,l), 0, bty='l', ylim=c(ymin,1), yaxt='n', pch=20, col=0, xlab=xlab, ylab=ylab)
#    title(ttt, line=1)
#    if(!is.null(br)){
#     abline(v=br, lwd=3, col="blue")
#    }
#    if(species == "wheat"){
#     z[, points(pos/1e6, sorted_chr, col=ifelse(sorted_alphachr == "3B" | sorted_arm == "S", 1, 2), pch=20)]
#     legend(horiz=T, "bottomleft", pch=19, col=1:2, bg="white", legend=c("short arm / 3B", "long arm"))
#    }
#    if(species == "rye"){
#     z[, points(pos/1e6, sorted_chr, col=1, pch=20)]
#    }
#    if(species == "barley"){
#     z[, points(pos/1e6, sorted_chr, col=ifelse(sorted_alphachr == "1H" | sorted_arm == "S", 1, 2), pch=20)]
#     if(mbscale == F){
#      llx =  "short arm / 1H"
#     } else {
#      llx = "short arm"
#     }
#     legend(horiz=T, "bottomleft", pch=19, col=1:2, bg="white", legend=c(llx, "long arm"))
#    }
#    axis(2, las=1, wheatchr$chr, wheatchr$alphachr)
#   } else {
#    plot(0, type='n', xlab="", ylab="", axes=F)
#   }
#
#   if(hic){
#    par(mar=c(4,4,4,1))
#    plot(0, xlim=c(0,l), col=0, bty='l', ylab='Hi-C chromosome', yaxt='n', ylim=c(ymin,1), xlab="position in scaffold (Mb)")
#    title("Interchromosomal Hi-C links", line=1)
#    tcc_pos[i, on="scaffold1", points(pos1/1e6, col="#00000003", chr2, pch=20)]
#    abline(v=br, lwd=3, col="blue")
#    axis(2, las=1, wheatchr$chr, wheatchr$alphachr)
#
#    span[i, on="scaffold"]->w
#    if(nrow(w[!is.na(n)]) == 0){
#     plot(xlab='', ylab='', type='n', 0, axes=F)
#    } else {
#     w[order(bin), plot(bty='l', bin/1e6, n, log='y', xlim=c(0,l), col=0, ylab='coverage', las=1, xlab="position in scaffold (Mb)")]
#     title("Intrascaffold physical Hi-C coverage", line=1)
#     abline(v=br, lwd=3, col="blue")
#     w[order(bin), points(bin/1e6, n, xlim=c(0,l), type='l', lwd=3, col=1)]
#    }
#
#    if(nrow(w[!is.na(n)]) == 0){
#     plot(xlab='', ylab='', type='n', 0, axes=F)
#    } else {
#     w[order(bin), plot(bty='l', bin/1e6, r, xlim=c(0,l), col=0, ylab='log2(observed/expected ratio)', las=1, xlab="position in scaffold (Mb)")]
#     title("Hi-C expected vs. observed coverage", line=1)
#     abline(v=br, lwd=3, col="blue")
#     w[order(bin), points(bin/1e6, r, xlim=c(0,l), type='l', lwd=3, col=1)]
#    }
#   }
#
#   if(tenex){
#    molcov[i, on="scaffold"]->w
#    if(nrow(w) == 0){
#     plot(xlab='', ylab='', type='n', 0, axes=F)
#    } else {
#     w[order(bin), plot(bty='l', bin/1e6, n, log='y', xlim=c(0,l), col=0, ylab='coverage', las=1, xlab="position in scaffold (Mb)")]
#     title("10X molecule coverage", line=1)
#     if(!is.null(br)){
#      abline(v=br, lwd=3, col="blue")
#     }
#     w[order(bin), points(bin/1e6, n, xlim=c(0,l), type='l', lwd=3, col=1)]
#    }
#
#    if(nrow(w) == 0){
#     plot(xlab='', ylab='', type='n', 0, axes=F)
#    } else {
#     w[order(bin), plot(bty='l', bin/1e6, r, xlim=c(0,l), col=0, ylab='log2(observed/expected ratio)', las=1, xlab="position in scaffold (Mb)")]
#     title("10X expected vs. observed coverage", line=1)
#     if(!is.null(br)){
#      abline(v=br, lwd=3, col="blue")
#     }
#     w[order(bin), points(bin/1e6, r, xlim=c(0,l), type='l', lwd=3, col=1)]
#    }
#   }
#  }