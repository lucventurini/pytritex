from ..utils import second_agg
import dask.dataframe as dd
from dask.distributed import Client
from dask.delayed import delayed


def add_hic_statistics(anchored_css: dd.DataFrame, fpairs: dd.DataFrame, client: Client, verbose=False):
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
    anchored0 = anchored_css.query("popseq_chr == sorted_chr")[["popseq_chr"]]  # .compute().reset_index(drop=False)
    anchored0 = anchored0.rename(columns={"popseq_chr": "chr"}).astype({"chr": str})
    anchored1 = anchored0[:].rename(columns={"chr": "chr1"})
    anchored1.index = anchored1.index.rename("scaffold_index1")
    fpairs1 = delayed(dd.merge)(
        client.scatter(fpairs), client.scatter(anchored1), on="scaffold_index1", how="left")
    anchored2 = anchored0[:].rename(columns={"chr": "chr2"})
    anchored2.index = anchored2.index.rename("scaffold_index2")
    fpairs2 = delayed(dd.merge)(fpairs1, client.scatter(anchored2), on="scaffold_index2", how="left")
    anchored_hic_links = client.compute(fpairs2).result()
    assert isinstance(anchored_hic_links, dd.DataFrame)

    hic_stats = anchored_hic_links.query("chr1 == chr1").rename(
        columns={"scaffold_index2": "scaffold_index", "chr1": "hic_chr"})
    assert "scaffold_index" in hic_stats.columns
    hic_stats = hic_stats.groupby(
        ["scaffold_index", "hic_chr"]).size().to_frame("N").reset_index(drop=False)
    # grouped_hic_stats = hic_stats.sort_values(["scaffold_index", "N"], ascending=[True, False])
    hic_ns = hic_stats.groupby("scaffold_index")["N"].sum(split_out=hic_stats.npartitions).to_dask_array()
    meta = dict(hic_stats.dtypes)
    grouped_hic_stats = hic_stats.groupby("scaffold_index").apply(
        lambda df: df.nlargest(2, "N"), meta=meta).reset_index(drop=True).groupby("scaffold_index")
    hic_stats = grouped_hic_stats.agg({"hic_chr": ["first", second_agg], "N": ["first", second_agg]},
                                      split_out=hic_stats.npartitions)
    hic_stats.columns = ["hic_chr", "hic_chr2", "hic_N1", "hic_N2"]
    hic_stats["Nhic"] = hic_ns
    # Now merge back to anchored_css
    hic_stats["hic_pchr"] = hic_stats["hic_N1"].div(hic_stats["Nhic"], fill_value=0)
    hic_stats["hic_p12"] = hic_stats["hic_N2"].div(hic_stats["hic_N1"], fill_value=0)
    func = delayed(dd.merge)(client.scatter(hic_stats),
                             client.scatter(anchored_css), on="scaffold_index", how="right")
    anchored_css = client.compute(func).result()
    return anchored_css, anchored_hic_links
