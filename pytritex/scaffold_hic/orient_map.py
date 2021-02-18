import dask.dataframe as dd
import pandas as pd
import numpy as np
import os
from scipy.stats import spearmanr
from typing import Union


def orient_hic_map(hic_map_bin: Union[pd.DataFrame, dd.DataFrame]):

    w = hic_map_bin.query("hic_bin == hic_bin")[["pos", "chr", "hic_bin"]]

    # w<-hic_map_bin[!is.na(hic_bin), .(id=scaffold, scaffold=sub(":.*$", "", scaffold), pos=as.integer(sub("^.*:", "", scaffold)), chr, hic_bin)]
    #     w<-w[, .(gbin=mean(na.omit(hic_bin)),
    # 	     hic_cor=as.numeric(suppressWarnings(cor(method='s', hic_bin, pos, use='p')))), keyby=scaffold][!is.na(hic_cor)]
    #     hic_map[!is.na(hic_bin) & scaffold %in% w$scaffold][order(chr, hic_bin)]->z0
    #     z0[,.(scaffold1=scaffold[1:(.N-2)], scaffold2=scaffold[2:(.N-1)], scaffold3=scaffold[3:(.N)]), by=chr]->z
    #     z0[, data.table(key="scaffold1", scaffold1=scaffold, hic_bin1=hic_bin)][setkey(z, "scaffold1")]->z
    #     z0[, data.table(key="scaffold2", scaffold2=scaffold, hic_bin2=hic_bin)][setkey(z, "scaffold2")]->z
    #     z0[, data.table(key="scaffold3", scaffold3=scaffold, hic_bin3=hic_bin)][setkey(z, "scaffold3")]->z
    #     w[, data.table(key="scaffold1", scaffold1=scaffold, gbin1=gbin)][setkey(z, "scaffold1")]->z
    #     w[, data.table(key="scaffold2", scaffold2=scaffold, gbin2=gbin)][setkey(z, "scaffold2")]->z
    #     w[, data.table(key="scaffold3", scaffold3=scaffold, gbin3=gbin)][setkey(z, "scaffold3")]->z
    #     z[, cc:= apply(z[, .(hic_bin1, hic_bin2, hic_bin3, gbin1, gbin2, gbin3)],1,function(x) {
    # 		    suppressWarnings(cor(x[1:3], x[4:6]))
    # 		    })]
    #     z[, data.table(key="scaffold", scaffold=scaffold2, cc=ifelse(cc > 0, 1, -1))]->ccor
    #     ccor[w]->m
    #     m[, hic_orientation:=ifelse(hic_cor > 0, 1 * cc, -1 * cc)]
    #     m[, .(scaffold, hic_cor, hic_invert=cc, hic_orientation)][hic_map, on="scaffold"]->hic_map_oriented
    #     hic_map_oriented[is.na(hic_orientation), hic_orientation := ifelse(hic_cor > 0, 1, -1)]
    #     setnames(hic_map_oriented, "chr", "consensus_chr")
    #     setnames(hic_map_oriented, "cM", "consensus_cM")
    #     hic_map_oriented[, consensus_orientation := hic_orientation]