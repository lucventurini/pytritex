import pandas as pd
import numpy as np
from pytritex.graph_utils.make_super_scaffolds import make_super_scaffolds
import time
from dask.distributed import Client
import dask.dataframe as dd
from dask.delayed import delayed
import dask.array as da


def _scaffold_unanchored(links: str,
                         excluded,
                         membership: str,
                         info: str,
                         sample_count: str,
                         save_dir: str,
                         client: Client,
                         ncores=1,
                         verbose=False):

    if verbose:
        print(time.ctime(), "Starting to unanchor unassigned scaffolds")

    # # use unanchored scaffolds to link super-scaffolds
    # ww2[is.na(popseq_chr1), .(scaffold_link=scaffold1, link_length=length1, scaffold1=scaffold2)]->x
    # ww2[is.na(popseq_chr1), .(scaffold_link=scaffold1, scaffold2=scaffold2)]->y
    # x[y, on="scaffold_link", allow.cartesian=T][scaffold1 != scaffold2]->xy
    #
    # m[, .(scaffold1=scaffold, super1=super, chr1=chr, cM1=cM, size1=super_nbin, d1 = pmin(bin - 1, super_nbin - bin))][xy, on="scaffold1"]->xy
    # m[, .(scaffold2=scaffold, super2=super, chr2=chr, cM2=cM, size2=super_nbin, d2 = pmin(bin - 1, super_nbin - bin))][xy, on="scaffold2"]->xy
    # xy[super2 != super1 & d1 == 0 & d2 == 0 & size1 > 1 & size2 > 1 & chr1 == chr2]->xy
    # xy[scaffold1 < scaffold2, .(nscl=.N), scaffold_link][xy, on="scaffold_link"]->xy
    # xy[nscl == 1] -> xy
    # xy[super1 < super2][, c("n", "g"):=list(.N, .GRP), by=.(super1, super2)][order(-link_length)][!duplicated(g)]->zz


    sample_count = dd.read_parquet(sample_count, infer_divisions=True)
    links = dd.read_parquet(links, infer_divisions=True)
    membership = dd.read_parquet(membership, infer_divisions=True)
    info = dd.read_parquet(info, infer_divisions=True)

    left_linkage = sample_count[sample_count["popseq_chr1"].isna()][["scaffold_index1", "length1", "scaffold_index2"]]
    left_linkage = left_linkage.rename(
        columns={"scaffold_index1": "scaffold_link",
                 "length1": "link_length",
                 "scaffold_index2": "scaffold_index1"})
    right_linkage = sample_count[sample_count["popseq_chr1"].isna()][
        ["scaffold_index1", "scaffold_index2"]].rename(
        columns={"scaffold_index1": "scaffold_link"})
    linkages = dd.merge(left_linkage, right_linkage, on="scaffold_link")

    distances = membership[["super", "chr", "cM", "super_nbin", "bin"]]
    distances["d"] = da.from_array(
        np.minimum(distances["bin"].compute() - 1, (distances["super_nbin"] - distances["bin"]).compute()),
        chunks=tuple(distances.map_partitions(len).compute().values.tolist()))
    left1 = distances.copy().rename(columns=dict((col, col + "1") for col in distances.columns))
    left1.index = left1.index.rename(distances.index.name + "1")
    linkages = dd.merge(linkages.set_index("scaffold_index1"), left1, on="scaffold_index1")
    left2 = distances.copy().rename(columns=dict((col, col + "2") for col in distances.columns))
    left2.index = left2.index.rename(distances.index.name + "2")
    linkages = dd.merge(linkages.reset_index(drop=False).set_index("scaffold_index2"),
                        left2, on="scaffold_index2").reset_index(drop=False)

    linkages = linkages.query(
        "super2 != super1 & d1 == 0 & d2 == 0 & super_nbin1 > 1 & super_nbin2 > 1 & chr1 == chr2")
    # Only keep unambiguous links
    linkages = linkages.merge(
        linkages.query("scaffold_index1 < scaffold_index2").groupby("scaffold_link").size().to_frame("nscl"),
        on="scaffold_link").query("nscl == 1 & super1 < super2")

    # xy[super1 < super2][, c("n", "g"):=list(.N, .GRP), by=.(super1, super2)][order(-link_length)][!duplicated(g)]->zz
    # sel <- zz[, .(scaffold1=c(scaffold_link, scaffold_link, scaffold1, scaffold2),
    # scaffold2=c(scaffold1, scaffold2, scaffold_link, scaffold_link))]
    # rbind(links, ww2[sel, on=c("scaffold1", "scaffold2")])->links2
    sizes = linkages.groupby(["super1", "super2"]).size().to_frame("n")
    sizes["g"] = range(1, sizes.shape[0].compute() + 1)



