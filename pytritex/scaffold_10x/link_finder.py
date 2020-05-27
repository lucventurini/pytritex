import pandas as pd


def _initial_link_finder(info: pd.DataFrame, molecules: pd.DataFrame,
                         verbose=False, popseq_dist=5,
                         min_npairs=2, max_dist=1e3, min_nmol=2, min_nsample=2):

    if verbose:
        print("Finding links")

    molecules_over_filter = molecules.query("npairs >= @min_npairs").set_index("scaffold_index")
    movf = molecules_over_filter
    lengths = info[["scaffold_index", "length"]].reset_index(drop=True).set_index("scaffold_index")
    lengths.columns = ["scaffold_length"]
    movf = movf.merge(lengths, left_index=True, right_index=True, how="left")
    movf.query("(end <= @max_dist) | (scaffold_length - start <= @max_dist)", inplace=True)
    # Group by barcode and sample. Only keep those lines in the table where a barcode in a given sample
    # is linking two different scaffolds.
    barcode_counts = movf[["barcode_index", "sample"]].groupby(
        ["barcode_index", "sample"], observed=True).size().to_frame("nsc")
    movf = barcode_counts.merge(movf, left_index=True,
                                right_on=barcode_counts.index.names).query("nsc >= 2").reset_index(drop=False)
    side_columns = ["scaffold_index", "npairs", "start", "end", "sample", "barcode_index"]
    # Only rename the first two columns
    scaffold1_side = movf.loc[:, side_columns].eval("pos1 = floor((start + end) / 2)")
    scaffold1_side = scaffold1_side.rename(columns=dict((_, _ + "1") for _ in side_columns[:2]), errors="raise")
    scaffold1_side = scaffold1_side.drop(["start", "end"], axis=1).set_index(["sample", "barcode_index"])

    scaffold2_side = movf.loc[:, side_columns].eval("pos2 = floor((start + end) / 2)")
    scaffold2_side = scaffold2_side.rename(columns=dict((_, _ + "2") for _ in side_columns[:2]))
    scaffold2_side = scaffold2_side.drop(["start", "end"], axis=1).set_index(["sample", "barcode_index"])
    # Now merge the two sides
    both_sides = scaffold1_side.merge(scaffold2_side, how="outer", left_index=True, right_index=True).apply(
        pd.to_numeric, downcast="signed").reset_index(drop=False)

    link_pos = both_sides.copy()
    # First aggregate on the molecules ("barcode_index") and select only those samples that have at least min_nmol
    mol_count = both_sides.groupby(
        ["scaffold_index1", "scaffold_index2", "sample"]
    ).agg(nmol=("barcode_index", "size")).query("nmol >= @min_nmol")
    # Then count how many samples pass the filter, and keep track of it.
    sample_count = mol_count[~mol_count.index.duplicated()].reset_index(
        level=2, drop=False).groupby(level=[0, 1]).agg(nsample=("sample", "size")).query("nsample >= @min_nsample")
    sample_count = sample_count.reset_index(level="scaffold_index2", drop=False)
    basic = info.loc[:, ["scaffold_index", "popseq_chr", "length", "popseq_pchr", "popseq_cM"]].set_index(
        "scaffold_index")
    left = basic[:].add_suffix("1")
    left.index.names = ["scaffold_index1"]
    sample_count = left.merge(sample_count, on="scaffold_index1", how="right")
    sample_count = sample_count.reset_index(drop=False).set_index("scaffold_index2")
    left = basic[:].add_suffix("2")
    left.index.names = ["scaffold_index2"]
    sample_count = left.merge(sample_count, on="scaffold_index2", how="right")
    sample_count.eval("same_chr = (popseq_chr2 == popseq_chr1)", inplace=True)
    sample_count.eval("weight = -1 * log10((length1 + length2) / 1e9)", inplace=True)
    sample_count = sample_count.reset_index(drop=False)

    if popseq_dist > 0:
        links = sample_count.query("(popseq_chr1 == popseq_chr2) & (abs(popseq_cM2 - popseq_cM1) <= @popseq_dist)")
    else:
        links = sample_count[:]
    return sample_count, links, link_pos
