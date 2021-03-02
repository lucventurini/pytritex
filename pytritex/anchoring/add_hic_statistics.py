from ..utils import second
import dask.dataframe as dd
from dask.distributed import Client
from dask.delayed import delayed
import numpy as np
import dask.array as da
import time
import logging
dask_logger = logging.getLogger("dask")


def add_hic_statistics(anchored_css: dd.DataFrame, fpairs: dd.DataFrame):
    """This function will add the HiC statistics to the anchored dataframe."""

    # info[!(popseq_chr != sorted_chr)][, .(scaffold, chr=popseq_chr)]->info0
    #   setnames(copy(info0), names(info0), sub("$", "1", names(info0)))[fpairs, on="scaffold1"]->tcc_pos
    #   setnames(copy(info0), names(info0), sub("$", "2", names(info0)))[tcc_pos, on="scaffold2"]->tcc_pos
    #   tcc_pos[!is.na(chr1), .N, key=.(scaffold=scaffold2, hic_chr=chr1)]->z
    #   z[order(-N)][, .(Nhic=sum(N), hic_chr=hic_chr[1], hic_N1=N[1],
    # 		   hic_chr2=hic_chr[2], hic_N2=N[2]), keyby=scaffold]->zz
    #   zz[, hic_pchr := hic_N1/Nhic]
    #   zz[, hic_p12 := hic_N2/hic_N1]
    #   zz[info, on="scaffold"]->info
    #   info[is.na(Nhic), Nhic := 0]
    #   info[is.na(hic_N1), hic_N1 := 0]
    #   info[is.na(hic_N2), hic_N2 := 0]
    dask_logger.debug("%s Starting to add the HiC statistics", time.ctime())
    anchored0 = anchored_css.query("popseq_chr == sorted_chr")[
        ["popseq_chr"]].compute().reset_index(drop=False)
    anchored0 = anchored0.rename(columns={"popseq_chr": "chr"}).astype({"chr": str})
    anchored1 = anchored0[:].rename(columns={"chr": "chr1"})
    anchored1.index = anchored1.index.rename("scaffold_index1")
    dask_logger.debug("%s Merging anchored_css to fpairs", time.ctime())
    nparts = fpairs.npartitions
    fpairs1 = dd.merge(fpairs, anchored1, on="scaffold_index1", how="left", npartitions=nparts)
    anchored2 = anchored0[:].rename(columns={"chr": "chr2"})
    anchored2.index = anchored2.index.rename("scaffold_index2")
    anchored_hic_links = dd.merge(fpairs1, anchored2, on="scaffold_index2", how="left", npartitions=nparts)
    dask_logger.debug("%s Merged anchored_css to fpairs", time.ctime())
    assert isinstance(anchored_hic_links, dd.DataFrame)

    dask_logger.debug("%s Creating the index columns", time.ctime())
    anchored_hic_links["hic_index"] = da.from_array(
        np.arange(1, anchored_hic_links.shape[0].compute() + 1),
        chunks=tuple(anchored_hic_links.map_partitions(len).compute().values.tolist()))
    dask_logger.debug("%s Indexing", time.ctime())
    anchored_hic_links = anchored_hic_links.set_index("hic_index")
    dask_logger.debug("%s Indexed", time.ctime())

    hic_stats = anchored_hic_links.query("chr1 == chr1").rename(
        columns={"scaffold_index2": "scaffold_index", "chr1": "hic_chr"})
    assert "scaffold_index" in hic_stats.columns
    dask_logger.debug("%s Calculating hic_stats", time.ctime())
    hic_stats = hic_stats.groupby(["scaffold_index", "hic_chr"]).size().to_frame("N")
    hic_stats = hic_stats.reset_index(drop=False).compute()
    # grouped_hic_stats = hic_stats.sort_values(["scaffold_index", "N"], ascending=[True, False])
    hic_ns = hic_stats.groupby("scaffold_index")["N"].sum()
    # meta = dict(hic_stats.dtypes)
    hic_stats = hic_stats.sort_values("N", ascending=False).groupby("scaffold_index").head(2)
    hic_stats = hic_stats.groupby("scaffold_index").agg({"hic_chr": ["first", second], "N": ["first", second]})
    hic_stats.columns = ["hic_chr", "hic_chr2", "hic_N1", "hic_N2"]
    hic_stats["Nhic"] = hic_ns    
    # Now merge back to anchored_css
    hic_stats["hic_pchr"] = hic_stats["hic_N1"].div(hic_stats["Nhic"], fill_value=0)
    hic_stats["hic_p12"] = hic_stats["hic_N2"].div(hic_stats["hic_N1"], fill_value=0)
    dask_logger.debug("%s Calculated hic_stats", time.ctime())
    anchored_css = dd.merge(anchored_css, hic_stats, on="scaffold_index", how="left",
                            npartitions=anchored_css.npartitions)
    dask_logger.debug("%s Merged hic_stats back into anchored_css", time.ctime())
    assert isinstance(anchored_css, dd.DataFrame)
    return anchored_css, anchored_hic_links
