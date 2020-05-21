import pandas as pd
from pytritex.utils import unique_count


def _initial_link_finder(info: pd.DataFrame, molecules: pd.DataFrame,
                         verbose=False, popseq_dist=1e4,
                         min_npairs=2, max_dist=1e3, min_nmol=2, min_nsample=2):

    if verbose:
        print("Finding links")

    z = molecules.loc[molecules["npairs"] >= min_npairs]
    try:
        z = info.loc[:, ["scaffold_index", "length"]].rename(columns={"length": "scaffold_length"}).merge(
            z, on="scaffold_index", how="right")
    except KeyError:
        raise KeyError((["scaffold_index", "length"], z.columns))
    z = z.loc[(z["end"] <= max_dist) | (z["scaffold_length"] - z["start"] <= max_dist), :]
    # Group by barcode and sample. Only keep those lines in the table where a barcode in a given sample
    # is linking two different scaffolds.
    grouped = z[["barcode_index", "sample", "scaffold_index"]
               ].drop_duplicates().groupby(["barcode_index", "sample"], observed=True)
    nsc = grouped["scaffold_index"].size().to_frame("nsc")
    z = nsc.merge(z, left_index=True, right_on=["barcode_index", "sample"]).loc[lambda df: df["nsc"] >= 2]
    x = z.loc[:, ["scaffold_index", "npairs", "start", "end", "sample", "barcode_index"]].rename(
        columns={"scaffold_index": "scaffold_index1", "npairs": "npairs1"}).assign(
        pos1=z.eval("floor((start + end) / 2)"))[["scaffold_index1", "npairs1", "pos1", "sample", "barcode_index"]]
    y = z.loc[:, ["scaffold_index", "npairs", "start", "end", "sample", "barcode_index"]].rename(
        columns={"scaffold_index": "scaffold_index2", "npairs": "npairs2"}).assign(
        pos2=z.eval("floor((start + end) / 2)"))[["scaffold_index2", "npairs2", "pos2", "sample", "barcode_index"]]
    xy = pd.merge(x, y, on=["sample", "barcode_index"], how="outer").apply(pd.to_numeric, downcast="signed")
    link_pos = xy.copy()
    w = xy.groupby(["scaffold_index1", "scaffold_index2", "sample"]).size().to_frame("nmol").reset_index(drop=False)
    ww = w.loc[w["nmol"] >= min_nmol]
    ww2 = ww.groupby(["scaffold_index1", "scaffold_index2"]).agg(nsample=("sample", unique_count)).loc[
        lambda df: df["nsample"] >= min_nsample].reset_index(drop=False)
    basic = info.loc[:, ["scaffold_index", "popseq_chr", "length", "popseq_pchr", "popseq_cM"]]
    ww2 = basic.rename(columns=dict((col, col + "1") for col in basic.columns)).merge(
        ww2, on=["scaffold_index1"], how="right")
    ww2 = basic.rename(columns=dict((col, col + "2") for col in basic.columns)).merge(
        ww2, on=["scaffold_index2"], how="right")
    ww2.loc[:, "same_chr"] = ww2["popseq_chr2"] == ww2["popseq_chr1"]
    ww2.loc[:, "weight"] = ww2.eval("-1 * log10((length1 + length2) / 1e9)")

    if popseq_dist > 0:
        links = ww2.query("(popseq_chr1 == popseq_chr2) & (abs(popseq_cM2 - popseq_cM1) <= @popseq_dist)")
    else:
        links = ww2

    return ww2, links, link_pos
