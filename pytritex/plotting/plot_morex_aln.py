import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def _plot_function(data: dict, membership: pd.DataFrame, aln, minlen=5e3):

    xx = aln.merge(data["super"], on=["super"])
    xx = xx.loc[(xx["agp_chr"] == xx["super_chr"]) & (xx["alnlen"] >= xx["minlen"])]
    if xx.shape[0] > 0:
        ylim = xx["agp_start"].quantile(q=np.arange(0, 1 + 1/50, 1/50)).iloc[[1, -2]] / 1e6
        xlim = np.array([0, data["length"]]) / 1e6
        # TODO: type='n' means no plotting. Ie, WHT??
        # # xx[, plot(0, xlim=xlim, type='n', ylim=ylim, las=1, xlab="10X super-scaffold position (Mb)",
        #     #           ylab="Morex AGP position (Mb)", bty='l')]
        # data[, title(main=paste0(sub("super", "super_scaffold", super), ", ", chr, "H, ", round(length / 1e6, 1), " Mb"))]
        if membership is not None:
            # abline(v=mem[data$super, on="super"]$super_pos / 1e6, col="gray")
            plt.axvline
        xx.loc[:, "idx"] = np.range(1, xx.shape[0] + 1, dtype=np.int)



    # xx[, idx := 1:.
    #     N]
    # xx[orientation == 1, lines(c(super_start / 1e6, super_end / 1e6), c(agp_start / 1e6, agp_end / 1e6)), by = idx]
    # xx[orientation == -1, lines(c(super_start / 1e6, super_end / 1e6), c(agp_end / 1e6, agp_start / 1e6)), by = idx]
    # } else {
    #     plot(0, type='n', axes=F, xlab="", ylab="")
    # data[, title(main=paste0(sub("super", "super_scaffold", super), ", ", chr, "H, ", round(length / 1e6, 1), " Mb",
    #                          "\nNo alignments to Morex pseudomolecules."))]
    # }
    # }



def plot_morex_aln(data: pd.DataFrame, morex_aln: pd.DataFrame, filename: str, ncores=1, minlen=5e3):

    data = data.copy().assign(idx=range(1, data.shape[0] + 1))
    paralle


plot_morex_aln<-function(data, morex_aln, file, ncores=1, minlen=5e3){
 copy(data)[, idx := 1:.N] -> data


 parallel_plot(data=data, group="super", cores=ncores, aln=morex_aln, minlen=minlen,
		file=file, height=700, width=700, res=150, plot_function=plotfu)
}