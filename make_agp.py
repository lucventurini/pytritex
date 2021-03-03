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
import glob
import json



def n50(array: np.array, p=0.5):
    array.sort()    
    return array[np.where(array.cumsum() >= (1 - p) * array.sum())[0][0]]


def create_agp(res, gap_size=100, do_write_agp=False):
    # row = res["row"]        
    row = dict()
    # npairs  nmol    nsample dist    folder  n50     n50_sup sum_sup
    # 2       3       1       90000   /home/lventuri/Tritex/wheat_uk_magic/hereward/scratch/dask-temp/../assembly/joblib/pytritex/scaffold_10x/239a82f057d0933746cde0298e68bff399a091da8c37be7ad04647f54676fb75       2689909 5186917 10347667022
    # ['min_npairs', 'max_dist', 'min_nmol', 'min_nsample', 'popseq_dist', 'max_dist_orientation', 'raw', 'info', 'molecules', 'fai', 'unanchored', 'result', 'hash', 'folder']
    for key in ["hash", "min_npairs", "min_nmol", "min_nsample", "max_dist", "popseq_dist", "max_dist_orientation"]:
        try:
            row[key] = res[key]
        except KeyError:
            raise KeyError(list(res.keys()))

    folder = res["folder"]
    result = dd.read_parquet(res["result"])
    row["folder"] = folder
    row["n50"] = n50(result.length.values.compute())
    row["num_sup"] = result.query("n > 1").shape[0].compute()    
    row["n50_sup"] = n50(result.query("n > 1").length.values.compute())
    row["sum_sup"] = result.query("n > 1").length.sum().compute()
    mem_name = os.path.join(row["folder"], "orientation", "membership")
    info_name = res["fai"]
    names = dd.read_parquet(os.path.join(os.path.dirname(res["fai"]), "cssaln"))[
        ["popseq_chr", "popseq_alphachr"]].drop_duplicates().compute().reset_index(drop=True)
    names = names.rename(columns={"popseq_chr": "chr", "popseq_alphachr": "alphachr"})
    if do_write_agp is True:
        write_agp(mem_name, info_name, names, folder, gap_size=gap_size)
    return list(row.items())
    
    
def write_agp(mem_name, info_name, names, folder, gap_size=100):
    membership = dd.read_parquet(mem_name)
    info = dd.read_parquet(info_name)
    agp_dict = make_agp(membership, info, names, gap_size=100)

    # Index(['super', 'super_start', 'super_end', 'super_index', 'gap',
    #       'orig_scaffold_index', 'orig_start', 'orig_end', 'orientation',
    #              'alphachr', 'cM', 'scaffold_index', 'length', 'super_length',
    #                     'super_size', 'scaffold', 'orig_scaffold_length', 'super_name'],
    #                           dtype='object')
    # 0       1       151729  1       False   16653   1       151729  1.0     7A      66.4445796324655     16653    151729  1643348 16      Triticum_aestivum_Cadenza_EIv1.1_scaffold_016653        151729  super_0
    # 0       151730  151829  2       True    -1      1       100                             -1      100  1643348  16      gap     100     super_0
    
    with open(f"{folder}/genome.agp", "wt", newline="") as agp_out:
        print(agp_out.name)
        agp = agp_dict["agp"].compute()
        _agp = agp[["super_name", "super_start", "super_end", "super_index",
                    "gap", "scaffold", "orig_start", "orig_end", "orientation",
                    "alphachr", "cM"]]
        is_gap = _agp["gap"] == True
        _agp.loc[is_gap, "orig_start"] = "scaffold"
        _agp.loc[is_gap, "orientation"] = "paired-ends"
        _agp.loc[is_gap, "gap"] = "U"
        _agp.loc[~is_gap, "gap"] = "W"
        assert isinstance(agp, (pd.DataFrame, dd.DataFrame)), type(agp)
        _agp.to_csv(agp_out, sep="\t", header=False, index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("results")
    parser.add_argument("-p", "--processes", default=10, type=int)
    parser.add_argument("--table-out", default="results.tsv", type=str)
    parser.add_argument("out", default="agp", type=str)
    args = parser.parse_args()

    dirs = [os.path.dirname(_) for _ in glob.glob(os.path.join(args.results, "*/orientation/res/_metadata"))]
    assert len(dirs) > 0
    results = []
    for key, val in enumerate(dirs, 1):
        result = json.load(open("/".join(val.split("/")[:6] + ["meta.json"] )))
        result["result"] = val
        result["hash"] = re.sub(".*scaffold_10x/([^/]*)/.*", "\\1", val)[:8]
        result["folder"] = re.sub("(.*scaffold_10x/[^/]*)/.*", "\\1", val)
        results.append(result)

    pool = mp.Pool(args.processes)
    rows = pool.map_async(create_agp, results)
    pool.close()
    pool.join()
    rows = rows.get()
    print(len(rows))
    assert len(rows) == len(dirs)
    rtable = pd.concat([pd.Series(dict(row)) for row in rows], axis=1).T.sort_values("n50", ascending=False)
    os.makedirs(os.path.dirname(os.path.abspath(args.table_out)), exist_ok=True)
    with open(args.table_out, "wt", newline="") as table_out:
        rtable.to_csv(table_out, sep="\t", header=True, index=False)

    # Now make AGP
    best = rtable.iloc[0]
    best_result = [result for result in results if result["folder"] == best["folder"]][0]
    create_agp(best_result, do_write_agp=True)
    agp_src = os.path.relpath(os.path.join(best["folder"], "genome.agp"), start=os.path.dirname(args.out))
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    if os.path.exists(args.out):
        os.remove(args.out)
    os.symlink(agp_src, args.out)
    return


if __name__ == "__main__":
    main()
