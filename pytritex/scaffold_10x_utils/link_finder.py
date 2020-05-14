import pandas as pd
import numpy as np


def _initial_link_finder(info: pd.DataFrame, molecules: pd.DataFrame,
                         verbose=False, popseq_dist=1e4,
                         min_npairs=2, max_dist=1e3, min_nmol=2, min_nsample=2):

    if verbose:
        print("Finding links")

    z = molecules.loc[molecules["npairs"] >= min_npairs]
    z = info.loc[:, ["scaffold", "length"]].rename(columns={"length": "scaffold_length"}).merge(
        z, on="scaffold", how="right")
    z = z.loc[(z["end"] <= max_dist) | (z["scaffold_length"] - z["start"] <= max_dist), :]
    # Group by barcode and sample. Only keep those lines in the table where a barcode in a given sample
    # is linking two different scaffolds.
    z.loc[:, "nsc"] = z.groupby(["barcode", "sample"])["scaffold"].transform(
        lambda series: series.unique().shape[0])
    z = z.loc[z["nsc"] >= 2, :]
    x = z.loc[:, ["scaffold", "npairs", "start", "end", "sample", "barcode"]].rename(
        columns={"scaffold": "scaffold1", "npairs": "npairs1"}).assign(pos1=(z["start"] + z["end"]) // 2)[
        ["scaffold1", "npairs1", "pos1", "sample", "barcode"]]
    y = z.loc[:, ["scaffold", "npairs", "start", "end", "sample", "barcode"]].rename(
        columns={"scaffold": "scaffold2", "npairs": "npairs2"}).assign(pos2=(z["start"] + z["end"]) // 2)[
        ["scaffold2", "npairs2", "pos2", "sample", "barcode"]]
    xy = pd.merge(x, y, on=["sample", "barcode"], how="outer")
    link_pos = xy.copy()
    w = xy.groupby(["scaffold1", "scaffold2", "sample"]).size().to_frame("nmol").reset_index(drop=False)
    ww = w.loc[w["nmol"] >= min_nmol]
    ww2 = ww.groupby(["scaffold1", "scaffold2"]).agg(nsample=("sample", lambda series: series.unique().shape[0])).loc[
        lambda df: df["nsample"] >= min_nsample].reset_index(drop=False)
    basic = info.loc[:, ["scaffold", "popseq_chr", "length", "popseq_pchr", "popseq_cM"]]
    ww2 = basic.rename(columns=dict((col, col + "1") for col in basic.columns)).merge(
        ww2, on=["scaffold1"], how="right")
    ww2 = basic.rename(columns=dict((col, col + "2") for col in basic.columns)).merge(
        ww2, on=["scaffold2"], how="right")
    ww2.loc[:, "same_chr"] = ww2["popseq_chr2"] == ww2["popseq_chr1"]
    ww2.loc[:, "weight"] = -1 * np.log10((ww2["length1"] + ww2["length2"]) / 1e9)

    if popseq_dist > 0:
        links = ww2.loc[
            (ww2["popseq_chr2"] == ww2["popseq_chr1"]) &
            ((ww2["popseq_cM1"] - ww2["popseq_cM2"]).abs() <= popseq_dist)]
    else:
        links = ww2

    return ww2, links, link_pos
