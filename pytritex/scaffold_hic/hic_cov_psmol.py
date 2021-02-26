import pandas as pd
import numpy as np
import dask.dataframe as dd
import functools
from ..sequencing_coverage import collapse_bins as cbn


def _group_analyser(group, binsize):
    """Count how many pairs cover a given bin; return into a column called "n"."""

    bins = group[["bin1", "bin2"]].to_numpy() + [0, 0]  # Need the sum to change the dtype
    counter = cbn(bins, binsize)
    counter = np.vstack([np.fromiter(counter.keys(), dtype=np.int), np.fromiter(counter.values(),
                                                                                dtype=np.int)]).T
    assigned = np.hstack([np.repeat(group.index.to_numpy()[0],
                                    counter.shape[0]).reshape(counter.shape[0], 1),
                          counter])
    return assigned


def hic_cov_psmol(hic_map: dict, binsize=int(1e3), binsize2=int(1e5), maxdist=int(1e6)):

    fpairs = hic_map["links"].copy()
    info = hic_map["chrlen"]
    fpairs = fpairs.rename(columns={"start1": "pos1", "start2": "pos2"})

    query = "chr1 == chr2 & pos1 < pos2"
    temp_values = fpairs.query(query)[["chr1", "pos1", "pos2"]].compute().to_numpy()
    temp_frame = pd.DataFrame({"chr": temp_values[:, 0],
                               "bin1": temp_values[:, 1] // binsize * binsize,
                               "bin2": temp_values[:, 2] // binsize * binsize}).astype(
        dtype={"chr": np.int32, "bin1": np.int32, "bin2": np.int32}).query(
        "(bin2 - bin1 > 2 * @binsize) & (bin2 - bin1 <= @maxdist)").set_index("chr")
    temp_frame.loc[:, "bin_group"] = temp_frame["bin1"] // binsize2
    temp_frame = dd.from_pandas(temp_frame, npartitions=fpairs.npartitions)
    _gr = functools.partial(_group_analyser, binsize=binsize)
    finalised = temp_frame.groupby(["chr", "bin_group"]).apply(_gr, meta=int).compute().values
    finalised = np.vstack(finalised)
    coverage_df = pd.DataFrame().assign(chr=finalised[:, 0],
                                        bin=finalised[:, 1],
                                        n=finalised[:, 2])
    coverage_df = coverage_df.groupby(
        ["chr", "bin"])["n"].sum().to_frame("n").reset_index(level=1)
    if info.index.name == "chr":
        left = info[["length"]]
    else:
        left = info[["chr", "length"]].set_index("chr")
    coverage_df = dd.merge(left, coverage_df, on="chr", how="right")
    arr = coverage_df[["bin", "length"]].to_dask_array(lengths=True)
    distance = dd.from_array(np.minimum(
        arr[:, 0],
        (arr[:, 1] - arr[:, 0]) // binsize * binsize)).to_frame("d")
    distance.index = coverage_df.index
    coverage_df = coverage_df.assign(d=distance["d"])
    assert coverage_df.query("d != d").shape[0].compute() == 0
    coverage_df["nbin"] = coverage_df.groupby(
        "chr").transform("size", meta=int).to_dask_array(lengths=True)
    mn = coverage_df.reset_index(drop=True)[["d", "n"]].groupby("d")["n"].mean().to_frame("mn")
    coverage_df = coverage_df.reset_index(drop=False)
    coverage_df = coverage_df.drop(
        "mn", axis=1, errors="ignore").set_index("d").merge(mn, how="left", on="d")
    coverage_df = coverage_df.reset_index(drop=False).set_index("chr")
    coverage_df = coverage_df.eval("r = log(n/mn) / log(2)")
    return coverage_df
