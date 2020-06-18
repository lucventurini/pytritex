import pandas as pd
import dask.dataframe as dd
import dask.array as da
import time
from ..utils import _rebalance_ddf
import os
import numpy as np
from dask.delayed import delayed


def _initial_link_finder(info: str, molecules: str, fai: str,
                         save_dir: str,
                         client,
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
    molecules_over_filter = molecules_over_filter.astype(
        dict((col, np.int32) for col in molecules_over_filter.columns)
    )
    molecules_over_filter = molecules_over_filter[molecules_over_filter["npairs"] >= min_npairs]
    # rebalance
    # molecules_over_filter = _rebalance_ddf(molecules_over_filter, target_memory=5 * 10**7)

    # We need to repartition. 100 partitions are just too few.
    movf = molecules_over_filter
    nparts = movf.npartitions
    movf = client.scatter(movf)
    fai = dd.read_parquet(fai, infer_divisions=True)
    lengths = fai[["length"]]
    lengths.columns = ["scaffold_length"]
    merger = delayed(dd.merge)(movf, lengths, left_index=True, right_index=True, how="left")
    movf = client.compute(merger).result()
    print(movf.head(npartitions=-1))

    # movf = dd.merge(movf, lengths, left_index=True, right_index=True, how="left")
    # assert "scaffold_length" in movf.columns
    movf = movf[(movf["end"] <= max_dist) | (movf.eval("scaffold_length - start") <= max_dist)]

    # Group by barcode and sample. Only keep those lines in the table where a barcode in a given sample
    # is linking two different scaffolds.
    movf = client.scatter(movf)
    def chunk_computer(df):
        return df.map_partitions(len)
    chunks = delayed(chunk_computer)(movf)
    chunks = tuple(client.compute(chunks).result().compute().values.tolist())
    print("Chunks:", chunks, type(chunks))

    def pair_counter(df):
        return df[["barcode_index", "sample", "npairs"]].astype(np.int32).compute().groupby(
                      ["barcode_index", "sample"])["npairs"].transform("size")

    nsc = delayed(pair_counter)(movf)
    nsc = client.compute(nsc).result().values
    nsc = da.from_array(nsc, chunks=chunks)
    movf = client.gather(movf)
    movf = movf.assign(nsc=nsc).query("nsc >= 2")
    print(movf.head(npartitions=-1))
    # These are the columns we are interested in.
    side_columns = ["npairs", "start", "end", "sample", "barcode_index"]
    # Only rename the first two columns
    assert "start" in movf.columns
    base = movf.loc[:, side_columns].reset_index(drop=False)
    base["pos"] = base.map_partitions(lambda df: df.eval("floor((start + end) / 2)"))
    base = base.drop(["start", "end"], axis=1).persist()
    # Now merge the two sides
    link_pos = dd.merge(base, base, how="outer",
                        left_on=["sample", "barcode_index"],
                        right_on=["sample", "barcode_index"],
                        suffixes=("1", "2")).persist()
    link_pos = link_pos.query(
        "(scaffold_index1 != scaffold_index2) | ((scaffold_index1 == scaffold_index2) & (pos1 != pos2))")
    assert "scaffold_index2" in link_pos, link_pos.head()
    assert "npairs2" in link_pos
    # link_pos = scaffold1_side.merge(scaffold2_side, how="outer", on=["sample", "barcode_index"])
    # link_pos = both_sides.copy()
    # link_pos = _rebalance_ddf(link_pos, target_memory= 5 * 10**7)
    link_pos_name = os.path.join(save_dir, "link_pos")
    dd.to_parquet(link_pos, link_pos_name, compression="gzip", engine="pyarrow")
    del link_pos
    # Reload from disk
    link_pos = dd.read_parquet(link_pos_name, infer_divisions=True)
    # First aggregate on the molecules ("barcode_index") and select only those
    # samples that have at least min_nmol
    # computed = both_sides.compute()
    print(time.ctime(), "Arrived at merging both sides")
    mol_count = link_pos[
        ["scaffold_index1", "scaffold_index2", "sample", "barcode_index"]]
    mol_count = mol_count.groupby(["scaffold_index1", "scaffold_index2", "sample"]
    )["barcode_index"].size().to_frame("nmol")
    mol_count = mol_count[mol_count["nmol"] >= min_nmol]
    # Then count how many samples pass the filter, and keep track of it.
    sample_count = mol_count.reset_index(drop=False).drop_duplicates(
        subset=["scaffold_index1", "scaffold_index2", "sample"]).groupby(
        ["scaffold_index1", "scaffold_index2"])["sample"].size().rename("nsample")
    sample_count = sample_count[sample_count >= min_nsample].persist()
    print("Sample count")
    print(sample_count.head(npartitions=-1))
    sample_count = sample_count.reset_index(drop=False).set_index("scaffold_index1", sorted=True)
    assert sample_count.index.name == "scaffold_index1", sample_count.head()
    print("Setting index on SC")
    print(sample_count.head(npartitions=-1))
    basic = info.loc[:, ["popseq_chr", "length", "popseq_pchr", "popseq_cM"]]
    left = basic.rename(columns=dict((col, col + "1") for col in basic.columns))
    left.index = left.index.rename("scaffold_index1")
    assert left.index.name == sample_count.index.name, (left.head(), sample_count.head())
    left = client.scatter(left)
    merger = delayed(dd.merge)(left, right=sample_count, left_index=True,
                               right_index=True, how="right")
    sample_count = client.compute(merger).result()
    print("Partial sample count")
    print(sample_count.head(npartitions=-1))
    sample_count = sample_count.reset_index(drop=False).set_index("scaffold_index2", sorted=False)
    left = client.gather(left)
    left = basic.rename(columns=dict((col, col + "2") for col in basic.columns))
    left.index = left.index.rename("scaffold_index2")
    left = client.scatter(left)
    merger = delayed(dd.merge)(left, sample_count,
                               left_index=True, right_index=True, how="right")
    sample_count = client.compute(merger).result()
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
    # sample_count = _rebalance_ddf(sample_count, target_memory=5 * 10**7)
    sample_count_name = os.path.join(save_dir, "sample_count")
    dd.to_parquet(sample_count, sample_count_name, compression="gzip", engine="pyarrow")
    del sample_count
    sample_count = dd.read_parquet(sample_count_name, infer_divisions=True)

    if popseq_dist > 0:
        links = sample_count.loc[
                (sample_count["same_chr"] == True) & (
                (sample_count.eval("popseq_cM2 - popseq_cM1").astype(float)).abs() <= popseq_dist), :]
        print("Links: ")
        print(links.head(npartitions=-1))
        links_name = os.path.join(save_dir, "links")
        # links = _rebalance_ddf(links, target_memory=5 * 10**7)
        dd.to_parquet(links, links_name, compression="gzip", engine="auto")
    else:
        # links = sample_count[:]
        links_name = sample_count_name[:]

    return sample_count_name, links_name, link_pos_name
