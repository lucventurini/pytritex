import pandas as pd
import numpy as np
from scipy.stats import median_absolute_deviation
from ..utils import second_agg
import dask.dataframe as dd
import itertools as it
from dask.distributed import Client
from dask.delayed import delayed
import time
import logging
dask_logger = logging.getLogger("dask")


def _mad1(chunks):
    c = chunks.apply(list)
    return c


def _mad2(grouped):
    # First create a list for each of the grouped chunks
    def internal(c):
        if (c != c).all():
            return [np.nan]
        f = [_ for _ in c if _ == _]
        f = [_ if isinstance(_, list) else [_] for _ in f]
        return list(it.chain.from_iterable(f))
    chunks = grouped.apply(internal)
    return chunks


def _mad3(chunks):
    chunks = chunks.apply(lambda s: np.nan if len(s) == 0 else
                          median_absolute_deviation(s, nan_policy="omit"))
    return chunks


mad = dd.Aggregation("mad", chunk=_mad1, agg=_mad2, finalize=_mad3)


def assign_popseq_position(cssaln: pd.DataFrame, popseq: pd.DataFrame,
                           anchored_css: dd.DataFrame, wheatchr: pd.DataFrame,
                           client: Client):


    # # Assignment of POPSEQ genetic positions
    #  popseq[!is.na(popseq_alphachr), .(css_contig, popseq_alphachr, popseq_cM)][
    #            cssaln[, .(css_contig, scaffold)], on="css_contig", nomatch=0]->z
    #  z[, .(.N, popseq_cM=mean(popseq_cM), popseq_cM_sd=ifelse(length(popseq_cM) > 1, sd(popseq_cM), 0),
    #  popseq_cM_mad=mad(popseq_cM)), keyby=.(scaffold, popseq_alphachr)]->zz
    #  zz[, popseq_Ncss := sum(N), by=scaffold]->zz
    #  zz[order(-N)][, .(popseq_alphachr=popseq_alphachr[1], popseq_Ncss1=N[1],
    # 			 popseq_alphachr2=popseq_alphachr[2], popseq_Ncss2=N[2]), keyby=scaffold]->x
    #  zz[, .(scaffold, popseq_alphachr, popseq_Ncss, popseq_cM, popseq_cM_sd,  popseq_cM_mad)][x, on=c("scaffold", "popseq_alphachr")]->zz
    #  zz[, popseq_pchr := popseq_Ncss1/popseq_Ncss]
    #  zz[, popseq_p12 := popseq_Ncss2/popseq_Ncss1]
    #  wheatchr[zz, on="popseq_alphachr"]->zz
    #  setnames(copy(wheatchr), paste0(names(wheatchr), 2))[zz, on="popseq_alphachr2"]->zz
    #  zz[info, on="scaffold"]->info
    #  info[is.na(popseq_Ncss), popseq_Ncss := 0]
    #  info[is.na(popseq_Ncss1), popseq_Ncss1 := 0]
    #  info[is.na(popseq_Ncss2), popseq_Ncss2 := 0]

    assert isinstance(popseq, dd.DataFrame)
    popseq_positions = popseq.loc[
                           ~popseq["popseq_alphachr"].isna(),
                           ["popseq_index", "popseq_chr", "popseq_cM"]
                       ][:].set_index("popseq_index").astype({"popseq_chr": int})
    assert isinstance(popseq_positions, dd.DataFrame)
    # This will be ["scaffold_index", "popseq_index"]
    dask_logger.debug("%s Assigning PopSeq - calculating right", time.ctime())
    right = cssaln[["popseq_index"]].reset_index(drop=False)
    dask_logger.debug("%s Assigning PopSeq - merging with the popseq_index", time.ctime())
    popseq_positions = dd.merge(popseq_positions, right, on="popseq_index", how="right",
                                npartitions=right.npartitions)
    # Group by "scaffold_index", "popseq_alphachr"
    dask_logger.debug("%s Assigning PopSeq (nparts: %s) - starting with popseq_stats", time.ctime(),
                        popseq_positions.npartitions)
    popseq_stats = popseq_positions.groupby(["scaffold_index", "popseq_chr"])
    popseq_count = popseq_stats.size(split_out=popseq_positions.npartitions).to_dask_array()
    dask_logger.debug("%s Assigning PopSeq (nparts %s) - calculating mean, std, MAD", time.ctime(),
                        popseq_positions.npartitions)
    popseq_stats = popseq_stats.agg({"popseq_cM": [np.mean, np.std, mad]}, split_out=popseq_positions.npartitions)
    popseq_stats.columns = ["popseq_cM", "popseq_cM_sd", "popseq_cM_mad"]
    popseq_stats["N"] = popseq_count
    dask_logger.debug("%s Assigning PopSeq - resetting the index for popseq_stats", time.ctime())
    popseq_stats = popseq_stats.reset_index(drop=False).set_index("scaffold_index", sorted=False)
    dask_logger.debug("%s Assigning PopSeq - computing Ncss", time.ctime())
    popseq_stats["popseq_Ncss"] = popseq_stats.groupby("scaffold_index")["N"].transform("sum", meta=int)
    # popseq_stats = popseq_stats.reset_index(drop=False)
    # So far we have congregated the different statistics about the centimorgans.
    # Now we have to consider .... ?
    popseq_stats = popseq_stats.reset_index(drop=False)
    _meta = dict(popseq_stats.dtypes)

    dask_logger.debug("%s Assigning PopSeq - starting to compute the best locations", time.ctime())
    best_locations = popseq_stats.groupby(
        "scaffold_index").apply(lambda df: df.nlargest(2, "N"), meta=_meta).reset_index(drop=True)

    best_locations = best_locations.groupby("scaffold_index").agg(
        {"popseq_chr": ["first", second_agg],
         "popseq_Ncss": ["first", second_agg]}
    ).astype({("popseq_chr", "first"): int, ("popseq_chr", "second"): float})
    best_locations.columns = ["popseq_chr", "popseq_chr2", "popseq_Ncss1", "popseq_Ncss2"]

    dask_logger.debug("%s Assigning PopSeq - merging popseq_stats with best_locations", time.ctime())
    popseq_stats = dd.merge(
        popseq_stats[
        ["scaffold_index", "popseq_chr",
         "popseq_Ncss", "popseq_cM", "popseq_cM_sd", "popseq_cM_mad"]],
        best_locations, on=["scaffold_index", "popseq_chr"], how="right", npartitions=best_locations.npartitions,
    )
    dask_logger.debug("%s Assigning PopSeq - calculating the popseq_pchr/12", time.ctime())
    popseq_stats["popseq_pchr"] = popseq_stats["popseq_Ncss1"].div(popseq_stats["popseq_Ncss"], fill_value=0)
    popseq_stats["popseq_p12"] = popseq_stats["popseq_Ncss2"].div(popseq_stats["popseq_Ncss1"], fill_value=0)

    dask_logger.debug("%s Assigning PopSeq - merging popseq_stats with wheatchr", time.ctime())
    pop1 = delayed(dd.merge)(
        popseq_stats,
        wheatchr.copy().rename(columns={"chr": "popseq_chr", "alphachr": "popseq_alphachr"}),
        on="popseq_chr", how="left")
    pop2 = delayed(dd.merge)(
        pop1,
        wheatchr.copy().rename(columns={"chr": "popseq_chr2", "alphachr": "popseq_alphachr2"}),
        how="left", on="popseq_chr2"
    )
    pop2 = client.compute(pop2).result().set_index("scaffold_index")
    del pop1
    dask_logger.debug("%s Assigning PopSeq - merging anchored_css with popseq_stats", time.ctime())
    anchored_css = dd.merge(anchored_css, pop2, on="scaffold_index", how="right", npartitions=anchored_css.npartitions)
    # anchored_css = anchored_css.set_index("scaffold_index")
    if anchored_css.index.name != "scaffold_index":
        print(anchored_css.columns)
        assert False
    anchored_css = anchored_css.reset_index(drop=False).set_index("scaffold_index")
    assert isinstance(anchored_css, dd.DataFrame)
    dask_logger.debug("%s Assigning PopSeq - finished", time.ctime())
    return anchored_css
