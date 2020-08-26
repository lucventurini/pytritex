import dask.dataframe as dd
import dask.array as da
import time
import os
import numpy as np
from dask.delayed import delayed
from dask.distributed import Client
import logging
dask_logger = logging.getLogger("dask")


def _calculate_link_pos(molecules: str, fai: str, save_dir: str,
                        client, min_npairs: int, max_dist: float):
    molecules_over_filter = dd.read_parquet(molecules, infer_divisions=True,
                                            filters=[
                                                ("npairs", ">=", min_npairs)
                                            ], engine="auto")
    molecules_over_filter = molecules_over_filter.astype(
        dict((col, np.int32) for col in molecules_over_filter.columns)
    )
    molecules_over_filter = molecules_over_filter.query("npairs >= @min_npairs",
                                                        local_dict={"min_npairs": min_npairs})[:]
    # We need to repartition. 100 partitions are just too few.
    movf = molecules_over_filter
    assert isinstance(movf, dd.DataFrame)
    fai = dd.read_parquet(fai, infer_divisions=True, engine="auto")
    lengths = fai.query("to_use == True")[["length"]]
    lengths.columns = ["scaffold_length"]
    movf = dd.merge(movf, lengths, left_index=True, right_index=True, how="left", npartitions=movf.npartitions)
    assert isinstance(movf, dd.DataFrame)
    movf = movf[(movf["end"] <= max_dist) | (movf.eval("scaffold_length - start") <= max_dist)]

    # Group by barcode and sample. Only keep those lines in the table where a barcode in a given sample
    # is linking two different scaffolds.

    def chunk_computer(df):
        return df.map_partitions(len)
    chunks = tuple(chunk_computer(movf).compute().values.tolist())
    def pair_counter(df):
        return df[["barcode_index", "sample", "npairs"]].astype(np.int32).compute().groupby(
            ["barcode_index", "sample"])["npairs"].transform("size")

    nsc = pair_counter(movf).values
    nsc = da.from_array(nsc, chunks=chunks)
    movf = movf.assign(nsc=nsc).query("nsc >= 2")
    # These are the columns we are interested in.
    side_columns = ["npairs", "start", "end", "sample", "barcode_index"]
    # Only rename the first two columns
    assert "start" in movf.columns
    base = movf.loc[:, side_columns].reset_index(drop=False)
    base["pos"] = base.map_partitions(lambda df: df.eval("floor((start + end) / 2)"))
    base = base.drop(["start", "end"], axis=1)
    # Now merge the two sides
    link_pos = dd.merge(base, base, how="outer",
                        left_on=["sample", "barcode_index"],
                        right_on=["sample", "barcode_index"],
                        suffixes=("1", "2"))
    link_pos = link_pos.query(
        "(scaffold_index1 != scaffold_index2) | ((scaffold_index1 == scaffold_index2) & (pos1 != pos2))")
    assert "scaffold_index2" in link_pos, link_pos.head()
    assert "npairs2" in link_pos
    # link_pos = scaffold1_side.merge(scaffold2_side, how="outer", on=["sample", "barcode_index"])
    # link_pos = both_sides.copy()
    link_pos_name = os.path.join(save_dir, "link_pos")
    link_pos = link_pos.persist()
    dd.to_parquet(link_pos, link_pos_name, compression="gzip", engine="pyarrow", compute=True, schema="infer")
    link_pos = dd.read_parquet(link_pos_name)
    return link_pos, link_pos_name


def mol_counter(link_pos, min_nmol):
    mol_count = link_pos.reset_index(drop=False)[
        ["scaffold_index1", "scaffold_index2", "sample", "barcode_index"]].drop_duplicates()
    mol_count = mol_count.groupby(
        ["scaffold_index1", "scaffold_index2", "sample"])["barcode_index"].size().to_frame("nmol")
    mol_count = mol_count.query("nmol >= @min_nmol", local_dict={"min_nmol": min_nmol})
    return mol_count


def sample_counter(mol_count, min_nsample):
    sample_count = mol_count.reset_index(drop=False)[
        ["scaffold_index1", "scaffold_index2", "sample"]].drop_duplicates()
    sample_count = sample_count.groupby(
        ["scaffold_index1", "scaffold_index2"])["sample"].size().to_frame("nsample")
    sample_count = sample_count.query("nsample >= @min_nsample", local_dict={"min_nsample": min_nsample})
    sample_count = sample_count.reset_index(drop=False).set_index("scaffold_index1")
    return sample_count


def index_resetter(sample_count, index_name, sorted=False):
    return sample_count.reset_index(drop=False).set_index(index_name, sorted=sorted)


def add_cols(sample_count):
    sample_count["same_chr"] = sample_count.eval(
        "((popseq_chr1 == popseq_chr2) & (popseq_chr1 == popseq_chr1) & (popseq_chr2 == popseq_chr2))").astype(int)
    sample_count["weight"] = sample_count.eval("-1 * log10((length1 + length2) / 1e9)")
    return sample_count


def _initial_link_finder(info: str, molecules: str, fai: str,
                         save_dir: str,
                         client: Client,
                         verbose=False, popseq_dist=5,
                         min_npairs=2, max_dist=1e5, min_nmol=2, min_nsample=2):

    os.makedirs(save_dir, exist_ok=True)
    info = dd.read_parquet(info, infer_divisions=True, engine="auto")

    #  assembly$info -> info
    #  assembly$molecules[npairs >= min_npairs] -> z
    #  info[, .(scaffold, scaffold_length=length)][z, on="scaffold"]->z
    #  z[end <= max_dist | scaffold_length - start <= max_dist]->z
    #  z[, nsc := length(unique(scaffold)), key=.(barcode, sample)]
    #  z[nsc >= 2]->z
    #  z[, .(scaffold1=scaffold, npairs1=npairs, pos1=as.integer((start+end)/2), sample, barcode)]->x
    #  z[, .(scaffold2=scaffold, npairs2=npairs, pos2=as.integer((start+end)/2), sample, barcode)]->y
    #  y[x, on=.(sample, barcode), allow.cartesian=T][scaffold1 != scaffold2]->xy
    #  xy -> link_pos
    #  xy[, .(nmol=.N), key=.(scaffold1, scaffold2, sample)]->w
    #  w[nmol >= min_nmol]->ww
    #  ww[, .(nsample = length(unique(sample))), key=.(scaffold1, scaffold2)][nsample >= min_nsample]->ww2
    #  info[, .(scaffold1=scaffold, popseq_chr1=popseq_chr, length1=length, popseq_pchr1=popseq_pchr, popseq_cM1=popseq_cM)][ww2, on="scaffold1"]->ww2
    #  info[, .(scaffold2=scaffold, popseq_chr2=popseq_chr, length2=length, popseq_pchr2=popseq_pchr, popseq_cM2=popseq_cM)][ww2, on="scaffold2"]->ww2
    #  ww2[, same_chr := popseq_chr2 == popseq_chr1]
    #  ww2[, weight := -1 * log10((length1 + length2) / 1e9)]
    #
    #  if(popseq_dist > 0){
    #   ww2[(popseq_chr2 == popseq_chr1 & abs(popseq_cM1 - popseq_cM2) <= popseq_dist)] -> links
    #  } else {
    #   ww2 -> links
    #  }

    # Reload from disk
    link_pos, link_pos_name = _calculate_link_pos(molecules, fai, save_dir, client, min_npairs, max_dist)
    dask_logger.warning("{} Arrived at merging both sides".format(time.ctime()))

    # link_pos = client.scatter(link_pos)
    sample_counts = []
    link_pos = link_pos.set_index("scaffold_index1")

    dask_logger.warning("{} Computing the molecule and sample counts".format(time.ctime()))
    mol_count = link_pos.map_partitions(mol_counter, min_nmol).persist()
    dask_logger.warning("{} Computed the molecule counts".format(time.ctime()))
    sample_count = mol_count.map_partitions(sample_counter, min_nsample).persist()
    dask_logger.warning("{} Computed the sample counts".format(time.ctime()))

    basic = info.loc[:, ["popseq_chr", "length", "popseq_pchr", "popseq_cM"]]
    left = basic.rename(columns=dict((col, col + "1") for col in basic.columns))
    left.index = left.index.rename("scaffold_index1")
    sample_count = dd.merge(left, sample_count, left_index=True, right_index=True, how="right",
                            npartitions=sample_count.npartitions).persist()
    sample_count = index_resetter(sample_count, "scaffold_index2", sorted=False)
    left = basic.rename(
        columns=dict((col, col + "2") for col in basic.columns))
    left.index = left.index.rename("scaffold_index2")
    sample_count = dd.merge(left, sample_count, left_index=True, right_index=True,
                            how="right", npartitions=sample_count.npartitions).persist()
    sample_count = add_cols(sample_count.reset_index(drop=False))
    dask_logger.warning("{} Added the last columns to sample_count.".format(time.ctime()))
    sample_count_name = os.path.join(save_dir, "sample_count")
    dd.to_parquet(sample_count, sample_count_name, compression="gzip",
                  engine="pyarrow", compute=True, schema="infer")
    dask_logger.warning("Calculating the link positions.")
    if popseq_dist > 0:
        sample_count = dd.read_parquet(sample_count_name, infer_divisions=True, engine="auto")
        links = sample_count[
                (sample_count["same_chr"] == True) &
                ((sample_count.eval("popseq_cM2 - popseq_cM1").astype(float)).abs() <= popseq_dist)]
        links_name = os.path.join(save_dir, "links")
        dd.to_parquet(links, links_name, compression="gzip", engine="pyarrow", compute=True, schema="infer")
        client.cancel(links)
    else:
        # links = sample_count[:]
        links_name = sample_count_name[:]

    client.cancel(sample_count)
    client.cancel(link_pos)
    return sample_count_name, links_name, link_pos_name
