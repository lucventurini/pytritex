import pandas as pd
import dask.dataframe as dd


def _initial_link_finder(info: dd.DataFrame, molecules: dd.DataFrame,
                         verbose=False, popseq_dist=5,
                         min_npairs=2, max_dist=1e5, min_nmol=2, min_nsample=2):

    if verbose:
        print("Finding links")

    molecules_over_filter = molecules.query("npairs >= @min_npairs", local_dict=locals())
    movf = molecules_over_filter
    lengths = info[["length"]]
    lengths.columns = ["scaffold_length"]
    movf = movf.merge(lengths, left_index=True, right_index=True, how="left")
    movf = movf.query("(end <= @max_dist) | (scaffold_length - start <= @max_dist)",
                      local_dict=locals())
    # Group by barcode and sample. Only keep those lines in the table where a barcode in a given sample
    # is linking two different scaffolds.
    barcode_counts = movf[["barcode_index", "sample"]].groupby(
        ["barcode_index", "sample"]).size().to_frame("nsc").reset_index(drop=False)
    # Only keep those cases where a given barcode has been confirmed in at least two samples.
    movf = barcode_counts.merge(
        movf.reset_index(drop=False), on=["barcode_index", "sample"]).query("nsc >= 2")

    # These are the columns we are interested in
    side_columns = ["scaffold_index", "npairs", "start", "end", "sample", "barcode_index"]
    # Only rename the first two columns
    assert "start" in movf.columns
    scaffold1_side = movf.loc[:, side_columns].eval("pos1 = floor((start + end) / 2)")
    scaffold1_side = scaffold1_side.rename(
        columns={"scaffold_index": "scaffold_index1",
                 "npairs": "npairs1"}).drop(["start", "end"], axis=1)
    assert "start" in movf.columns

    scaffold1_side.index.name = "scaffold_index1"
    scaffold2_side = movf.loc[:, side_columns].eval("pos2 = floor((start + end) / 2)")
    scaffold2_side = scaffold2_side.rename(
        columns={"scaffold_index": "scaffold_index2",
                 "npairs": "npairs2"}).drop(["start", "end"], axis=1)
    # Now merge the two sides
    both_sides = scaffold1_side.merge(scaffold2_side, how="outer", on=["sample", "barcode_index"])
    link_pos = both_sides.copy()

    # First aggregate on the molecules ("barcode_index") and select only those
    # samples that have at least min_nmol
    mol_count = both_sides.groupby(
        ["scaffold_index1", "scaffold_index2", "sample"]
    )["barcode_index"].agg("size").to_frame(
        "nmol").query("nmol >= @min_nmol", local_dict=locals()).compute()
    # Then count how many samples pass the filter, and keep track of it.
    sample_count = mol_count[~mol_count.index.duplicated()].reset_index(
        level=2, drop=False).groupby(level=[0, 1]).agg(
        nsample=("sample", "size")).query("nsample >= @min_nsample")
    sample_count = sample_count.reset_index(level="scaffold_index2", drop=False)
    assert sample_count.index.name == "scaffold_index1", sample_count.head()

    basic = info.loc[:, ["popseq_chr", "length", "popseq_pchr", "popseq_cM"]]
    left = basic.rename(columns=dict((col, col + "1") for col in basic.columns))
    left.index = left.index.rename("scaffold_index1")
    assert left.index.name == sample_count.index.name, (left.head(), sample_count.head())
    sample_count = dd.merge(left, sample_count, on="scaffold_index1", how="right")
    sample_count = sample_count.reset_index(drop=False).set_index("scaffold_index2")
    left = basic.rename(columns=dict((col, col + "2") for col in basic.columns))
    left.index = left.index.rename("scaffold_index2")
    sample_count = left.merge(sample_count, on="scaffold_index2", how="right")
    # Check that chromosomes are not NAs
    print("Positions without an assigned chromosome:",
          sample_count.loc[sample_count["popseq_chr1"].isna()].shape[0])
    sample_count = sample_count.eval(
        "same_chr = ((popseq_chr1 == popseq_chr2) & (popseq_chr1 == popseq_chr1) & (popseq_chr2 == popseq_chr2))"
    )
    sample_count = sample_count.eval("weight = -1 * log10((length1 + length2) / 1e9)")
    sample_count = sample_count.reset_index(drop=False)
    if popseq_dist > 0:
        links = sample_count.copy().loc[(sample_count["same_chr"] == True) & (
            (sample_count["popseq_cM2"] - sample_count["popseq_cM1"]).abs() <= popseq_dist), :]
    else:
        links = sample_count[:]

    return sample_count, links, link_pos
