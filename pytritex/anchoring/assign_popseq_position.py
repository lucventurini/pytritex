import pandas as pd
import numpy as np
from scipy.stats import median_absolute_deviation
from ..utils import second_agg
import dask.dataframe as dd
import itertools as it
from dask.distributed import Client
from dask.delayed import delayed


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
    chunks = chunks.apply(lambda s: median_absolute_deviation(s, nan_policy="omit"))
    return chunks


mad = dd.Aggregation("mad", chunk=_mad1, agg=_mad2, finalize=None)


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
                           ["popseq_index", "popseq_alphachr", "popseq_cM"]
                       ][:].set_index("popseq_index").astype({"popseq_alphachr": str})
    assert isinstance(popseq_positions, dd.DataFrame)
    # This will be ["scaffold_index", "popseq_index"]
    right = client.scatter(cssaln[["popseq_index"]].compute().reset_index(drop=False))
    popseq_positions = client.scatter(popseq_positions)
    func = delayed(dd.merge)(popseq_positions, right, on="popseq_index", how="left")
    popseq_positions = client.compute(func).result()
    popseq_positions = popseq_positions.astype({"popseq_alphachr": object})
    # Group by "scaffold_index", "popseq_alphachr"
    popseq_stats = popseq_positions.groupby(["scaffold_index", "popseq_alphachr"])
    popseq_count = popseq_stats.size(
        split_out=popseq_positions.npartitions).to_dask_array()
    popseq_stats = popseq_stats.agg({"popseq_cM": [np.mean, np.std, mad]},
                                    split_out=popseq_positions.npartitions)
    popseq_stats.columns = ["popseq_cM", "popseq_cM_sd", "popseq_cM_mad"]
    popseq_stats["N"] = popseq_count
    popseq_stats = popseq_stats.reset_index(drop=False).set_index("scaffold_index", sorted=False)
    popseq_stats["popseq_Ncss"] = popseq_stats.groupby("scaffold_index")["N"].transform("sum", meta=int)
    # popseq_stats = popseq_stats.reset_index(drop=False)
    # So far we have congregated the different statistics about the centimorgans.
    # Now we have to consider .... ?
    popseq_stats = popseq_stats.reset_index(drop=False).persist()
    _meta = dict(popseq_stats.dtypes)

    best_locations = popseq_stats.groupby(
        "scaffold_index").apply(lambda df: df.nlargest(2, "N"), meta=_meta).reset_index(drop=True)

    best_locations = best_locations.groupby("scaffold_index").agg(
        {"popseq_alphachr": ["first", second_agg],
         "popseq_Ncss": ["first", second_agg]}
    ).astype({("popseq_alphachr", "first"): object, ("popseq_alphachr", "second"): object})
    best_locations.columns = ["popseq_alphachr", "popseq_alphachr2", "popseq_Ncss1", "popseq_Ncss2"]

    popseq_stats = delayed(dd.merge)(
        client.scatter(popseq_stats[
        ["scaffold_index", "popseq_alphachr", "popseq_Ncss", "popseq_cM", "popseq_cM_sd", "popseq_cM_mad"]]),
        client.scatter(best_locations), on=["scaffold_index", "popseq_alphachr"], how="right"
    )
    popseq_stats = client.compute(popseq_stats).result()
    popseq_stats["popseq_pchr"] = popseq_stats["popseq_Ncss1"].div(popseq_stats["popseq_Ncss"], fill_value=0)
    popseq_stats["popseq_p12"] = popseq_stats["popseq_Ncss2"].div(popseq_stats["popseq_Ncss1"], fill_value=0)

    pop1 = delayed(dd.merge)(
        client.scatter(popseq_stats),
        wheatchr.copy().rename(columns={"chr": "popseq_chr", "alphachr": "popseq_alphachr"}),
        on="popseq_alphachr", how="left")
    pop2 = delayed(dd.merge)(
        pop1,
        wheatchr.copy().rename(columns={"chr": "popseq_chr2", "alphachr": "popseq_alphachr2"}),
        how="left", on="popseq_alphachr2"
    )

    pop2 = client.compute(pop2).result().set_index("scaffold_index")
    func = delayed(dd.merge)(client.scatter(pop2),
                             client.scatter(anchored_css), on="scaffold_index", how="right")
    anchored_css = client.compute(func).result()
    # anchored_css = anchored_css.set_index("scaffold_index")
    anchored_css = anchored_css.persist()
    if anchored_css.index.name != "scaffold_index":
        print(anchored_css.columns)
        assert False
    # for column in ["popseq_Ncss", "popseq_Ncss1", "popseq_Ncss2"]:
    #     anchored_css.loc[:, column] = pd.to_numeric(anchored_css[column].fillna(0),
    #                                                 downcast="signed")
    return anchored_css
