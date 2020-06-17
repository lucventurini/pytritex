import pandas as pd
import dask.dataframe as dd
import dask.array as da
import time
from ..utils import _rebalance_ddf
import os


def merger(left, right, left_on=None, right_on=None, left_index=False, right_index=False,
           how="inner"):
    return dd.merge(left, right, left_on=left_on, right_on=right_on,
                    left_index=left_index, right_index=right_index, how=how)


def _initial_link_finder(info: str, molecules: str, fai: str,
                         save_dir: str,
                         verbose=False, popseq_dist=5,
                         min_npairs=2, max_dist=1e5, min_nmol=2, min_nsample=2):

    if verbose:
        print("Finding links")

    os.makedirs(save_dir, exist_ok=True)
    info = dd.read_parquet(info, infer_divisions=True)
    molecules_over_filter = dd.read_parquet(molecules, infer_divisions=True,
                                filters=[
                                    ("npairs", ">=", min_npairs)
                                ])
    molecules_over_filter = molecules_over_filter[molecules_over_filter["npairs"] >= min_npairs]
    # rebalance
    molecules_over_filter = _rebalance_ddf(molecules_over_filter, target_memory=5 * 10**7)

    # We need to repartition. 100 partitions are just too few.
    movf = molecules_over_filter
    nparts = movf.npartitions
    fai = dd.read_parquet(fai, infer_divisions=True)
    lengths = fai[["length"]]
    lengths.columns = ["scaffold_length"]
    movf = dd.map_partitions(merger, movf, right=lengths,
                             left_index=True, right_index=True, how="left")

    # movf = dd.merge(movf, lengths, left_index=True, right_index=True, how="left")
    assert "scaffold_length" in movf.columns
    movf = movf.query("(end <= @max_dist) | (scaffold_length - start <= @max_dist)",
                      local_dict=locals()).reset_index(drop=False)

    # Group by barcode and sample. Only keep those lines in the table where a barcode in a given sample
    # is linking two different scaffolds.
    barcode_counts = movf[["barcode_index", "sample"]].groupby(
        ["barcode_index", "sample"]).size().to_frame("nsc").reset_index(drop=False)
    # let's write down this.
    print(time.ctime(), "Writing down the barcode counts")
    dd.to_parquet(barcode_counts, os.path.join(save_dir, "barcode_counts"),
                  compression="gzip", engine="pyarrow")
    del barcode_counts
    barcode_counts = dd.read_parquet(os.path.join(save_dir, "barcode_counts"),
                                     infer_divisions=True)
    print(time.ctime(), "Written down the barcode counts")
    # Only keep those cases where a given barcode has been confirmed in at least two samples.
    # movf = movf.map_partitions(lambda df: dd.merge(barcode_counts, df,
    #                                                how="right", on=["barcode_index", "sample"]))
    movf = dd.map_partitions(merger, movf, right=barcode_counts, how="right",
                             left_on=["barcode_index", "sample"],
                             right_on=["barcode_index", "sample"])
    assert "nsc" in movf.columns
    movf = movf[movf["nsc"] >= 2]
    assert isinstance(movf, dd.DataFrame)
    # These are the columns we are interested in
    side_columns = ["scaffold_index", "npairs", "start", "end", "sample", "barcode_index"]
    # Only rename the first two columns
    assert "start" in movf.columns
    scaffold1_side = movf.loc[:, side_columns]
    scaffold1_side["pos1"] = scaffold1_side.map_partitions(
        lambda df: df.eval("floor((start + end) / 2)")
    )
    scaffold1_side = scaffold1_side.rename(
        columns={"scaffold_index": "scaffold_index1",
                 "npairs": "npairs1"}).drop(["start", "end"], axis=1)
    assert "start" in movf.columns

    print(time.ctime(), "Prepared left (scaffold 1) side")
    scaffold2_side = movf.loc[:, side_columns]
    scaffold2_side["pos2"] = scaffold2_side.map_partitions(
        lambda df: df.eval("floor((start + end) / 2)")
    )
    scaffold2_side = scaffold2_side.rename(
        columns={"scaffold_index": "scaffold_index2",
                 "npairs": "npairs2"}).drop(["start", "end"], axis=1)
    # Now merge the two sides
    link_pos = dd.map_partitions(merger, scaffold1_side, right=scaffold2_side, how="outer",
                                 left_on=["sample", "barcode_index"],
                                 right_on=["sample", "barcode_index"])
    assert "scaffold_index2" in link_pos
    assert "npairs2" in link_pos
    # link_pos = scaffold1_side.merge(scaffold2_side, how="outer", on=["sample", "barcode_index"])
    # link_pos = both_sides.copy()
    link_pos = _rebalance_ddf(link_pos, target_memory= 5 * 10**7)
    link_pos_name = os.path.join(save_dir, "link_pos")
    dd.to_parquet(link_pos, link_pos_name, compression="gzip", engine="pyarrow")
    del link_pos
    # Reload from disk
    link_pos = dd.read_parquet(link_pos_name, infer_divisions=True)
    print(time.ctime(), "Prepared left (scaffold 2) side")

    # First aggregate on the molecules ("barcode_index") and select only those
    # samples that have at least min_nmol
    # computed = both_sides.compute()
    print(time.ctime(), "Arrived at merging both sides")
    mol_count = link_pos.groupby(
        ["scaffold_index1", "scaffold_index2", "sample"]
    )["barcode_index"].agg("size").to_frame("nmol").query("nmol >= @min_nmol", local_dict=locals())
    # Then count how many samples pass the filter, and keep track of it.
    sample_count = mol_count.reset_index(drop=False).drop_duplicates(
        subset=["scaffold_index1", "scaffold_index2", "sample"]).groupby(
        ["scaffold_index1", "scaffold_index2"])["sample"].size().rename("nsample")
    sample_count = sample_count[sample_count >= min_nsample]
    sample_count = sample_count.reset_index(drop=False).set_index("scaffold_index1", sorted=True)
    assert sample_count.index.name == "scaffold_index1", sample_count.head()
    basic = info.loc[:, ["popseq_chr", "length", "popseq_pchr", "popseq_cM"]]
    left = basic.rename(columns=dict((col, col + "1") for col in basic.columns))
    left.index = left.index.rename("scaffold_index1")
    assert left.index.name == sample_count.index.name, (left.head(), sample_count.head())
    sample_count = dd.map_partitions(merger, left, right=sample_count,
                                     left_index=True,
                                     right_index=True, how="right")
    sample_count = sample_count.reset_index(drop=False).set_index("scaffold_index2", sorted=False)
    left = basic.rename(columns=dict((col, col + "2") for col in basic.columns))
    left.index = left.index.rename("scaffold_index2")
    sample_count = dd.map_partitions(left, sample_count,
                                     left_index=True,
                                     right_index=True,
                                     how="right")
    # Check that chromosomes are not NAs
    print(time.ctime(),
          "Positions without an assigned chromosome:",
          sample_count.loc[sample_count["popseq_chr1"].isna()].shape[0].compute())

    sample_count["same_chr"] = sample_count.map_partitions(
        lambda df: df.eval(
            "((popseq_chr1 == popseq_chr2) & (popseq_chr1 == popseq_chr1) & (popseq_chr2 == popseq_chr2))"
        ))
    sample_count["weight"] = sample_count.map_partitions(
        lambda df: df.eval("-1 * log10((length1 + length2) / 1e9)")
    )

    sample_count = sample_count.reset_index(drop=False)
    sample_count = _rebalance_ddf(sample_count, target_memory=5 * 10**7)
    sample_count_name = os.path.join(save_dir, "sample_count")
    dd.to_parquet(sample_count, sample_count_name, compression="gzip", engine="pyarrow")
    del sample_count
    sample_count = dd.read_parquet(sample_count_name, infer_divisions=True)

    if popseq_dist > 0:
        links = sample_count.copy().loc[(sample_count["same_chr"] == True) & (
            (sample_count["popseq_cM2"] - sample_count["popseq_cM1"]).abs() <= popseq_dist), :]
        links_name = os.path.join(save_dir, "links")
        links = _rebalance_ddf(links, target_memory=5 * 10**7)
        dd.to_parquet(links, links_name, compression="gzip", engine="pyarrow")
    else:
        # links = sample_count[:]
        links_name = sample_count_name[:]

    return sample_count_name, links_name, link_pos_name
