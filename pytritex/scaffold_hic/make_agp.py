import numpy as np
import pandas as pd
import dask.dataframe as dd
from ..utils.chrnames import chrNames
import dask.array as da
from typing import Union


# make_agp<-function(hic_map_oriented, gap_size=100, species){
#
#  hic_map_oriented[, .(scaffold, chr = consensus_chr,
# 	  popseq_cM=ifelse(consensus_chr == popseq_chr | is.na(consensus_chr), popseq_cM, NA),
# 	  scaffold_length = length, hic_bin, orientation=consensus_orientation)]->z
#  chrNames(agp=T, species=species)[z, on="chr"]->z
#
#  z[, agp_chr := "chrUn"]
#  z[!is.na(hic_bin), agp_chr := sub("NA", "Un", paste0("chr", alphachr))]
#  z[, alphachr := NULL]
#  z[order(agp_chr, hic_bin, chr, popseq_cM, -scaffold_length)]->z
#  z[, index := 2*1:.N-1]
#  z[, gap := F]
#  rbind(z, data.table(scaffold="gap", gap=T, chr=NA, popseq_cM=NA, scaffold_length = gap_size, hic_bin = NA, orientation = NA, agp_chr=z$agp_chr, index=z$index+1))->z
#  z[order(index)][, head(.SD, .N-1), by=agp_chr]->z
#  z[, agp_start := cumsum(c(0, scaffold_length[1:(.N-1)]))+1, by = agp_chr]
#  z[, agp_end := cumsum(scaffold_length), by = agp_chr]
#
#  z[, .(scaffold=scaffold, bed_start=0, bed_end=scaffold_length, name=scaffold, score=1, strand=ifelse(is.na(orientation) | orientation == 1, "+", "-"), agp_chr=agp_chr)]->agp_bed
#
#  list(agp=z, agp_bed=agp_bed)
# }


def make_agp(hic_map_oriented: Union[pd.DataFrame, dd.DataFrame], gap_size=100, species=None):
    names = chrNames(agp=True, species=species)

    initial = hic_map_oriented.reset_index(drop=False)[
        ["scaffold_index", "consensus_chr", "popseq_cM", "length", "hic_bin", "consensus_orientation"]]

    initial = initial.rename(
        columns={"consensus_chr": "chr", "length": "scaffold_length",
                 "consensus_orientation": "orientation"})
    initial["popseq_cM"] = initial["popseq_cM"].mask(
        (initial["chr"] != initial["popseq_chr"]) & ~initial["chr"].isna(), np.nan)
    initial = initial.merge(names, on="chr", how="left")
    initial["index"] = da.from_array(np.arange(1, initial.shape[0].compute() + 1) * 2 - 1,
                                     chunks=initial.map_partitions(len).compute().values.tolist())
    initial = initial.compute()
    initial = initial.sort_values(["chr", "hic_bin", "popseq_cM", "scaffold_length"],
                                  ascending=[True, True, True, False])
    initial.loc[:, "agp_index"] = np.arange(1, initial.shape[0].compute() + 1, dtype=int) * 2 - 1
    initial.loc[:, "gap"] = False
    initial.loc[:, "start"] = 1
    gaps = pd.DataFrame().assign(chr=initial["chr"].values.compute(),
                                 hic_bin=np.nan,
                                 start=1,  # This is the "gap_type" column
                                 length=gap_size,
                                 orientation=np.nan,
                                 index=initial.index.values.compute() + 1,
                                 alphachr=np.nan,
                                 cM=np.nan,
                                 popseq_cM=np.nan,
                                 scaffold_index=-1,
                                 gap=True)
    with_gaps = (initial.groupby("chr").size() > 1).compute().to_frame("with_gaps")
    gaps = gaps.merge(with_gaps, on="chr", how="left")
    gaps = gaps.query("with_gaps == True")[:].drop("with_gaps", axis=1).set_index("index")
    assert (gaps["scaffold_index"].unique() == [-1]).all()
    data = dd.concat([initial, gaps]).reset_index(drop=False).set_index("index")
    data["chr_start"] = (data.groupby("chr")["length"].cumsum().shift(1).fillna(0) + 1).astype(int)
    data["chr_end"] = data.groupby("chr")["length"].cumsum().compute()
    data["chr_index"] = (data.groupby("chr").cumcount() + 1)

    # Now we have to get the names for the scaffolds
