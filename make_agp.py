#!/usr/bin/env python

import sys
import os
import dask.dataframe as dd
import matplotlib
matplotlib.use("agg")  # Avoid pesky Gdk error messages
import pandas as pd
import argparse
import joblib
import re
import numpy as np
from pytritex.scaffold_10x.make_agp import make_agp
import multiprocessing as mp



def n50(array: np.array, p=0.5):
    array.sort()
    return array[np.where(array.cumsum() >= (1 - p) * array.sum())[0][0]]


def write_agp(mem_name, info_name, folder, gap_size=100):
    membership = dd.read_parquet(mem_name)
    info = dd.read_parquet(info_name)
    agp_dict = make_agp(membership, info, gap_size=100)
    with open(f"{folder}/genome.agp", "wt", newline="") as agp_out:
        agp_dict["agp"].to_csv(agp_out, sep="\t", header=False, index=False)    


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("results")
    parser.add_argument("--table-out", default="results.tsv", type=str)
    parser.add_argument("out", default="agp", type=str)
    args = parser.parse_args()

    results = joblib.load(args.results)

    pool = mp.Pool(10)
    
    rows = []
    for r_index, res in enumerate(results):
        row = res["row"]        
        folder = re.sub("(.*scaffold_10x/[^/]*)/.*", "\\1", res["result"])
        result = dd.read_parquet(res["result"])
        row["folder"] = folder
        row["n50"] = n50(result.length.values.compute())
        row["n50_sup"] = n50(result.query("n > 1").length.values.compute())
        row["sum_sup"] = result.query("n > 1").length.sum().compute()
        rows.append(row)        
        import time
        print(time.ctime(), "Started writing", r_index, file=sys.stderr)
        write_agp(res["membership"], os.path.join(os.path.dirname(res["info"]), "fai"), folder, 100)
        print(time.ctime(), "Finished writing", r_index, file=sys.stderr)        
        # pool.apply_async(write_agp, args=(res["membership"], os.path.join(os.path.dirname(res["info"]), "fai"), folder, 100))
        del result

    pool.close()
    pool.join()
    rtable = pd.concat(rows, axis=1).T.sort_values("n50", ascending=False)
    with open(args.table_out, "wt", newline="") as table_out:
        rtable.to_csv(table_out, sep="\t", header=True, index=False)

    # Now make AGP
    best = rtable.iloc[0]
    agp_src = os.path.relpath(os.path.join(best["folder"], "genome.agp"))
    os.symlink(agp_src, "genome.agp")
    return


if __name__ == "__main__":
    main()
