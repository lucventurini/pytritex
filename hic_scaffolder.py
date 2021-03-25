import argparse
import matplotlib

from pytritex.scaffold_hic.bin_hic_step import bin_hic_step

matplotlib.use("agg")  # Avoid pesky Gdk error messages
import dask
dask.config.set({'distributed.worker.multiprocessing-method': 'spawn'})
import dask.dataframe as dd
from pytritex.scaffold_hic.init_10x_assembly import init_10x_assembly
from pytritex.scaffold_hic.hic_map import hic_map, calculate_hic_link_weights
from pytritex.anchoring import anchor_scaffolds
from pytritex.utils import n50
from pytritex.scaffold_hic.read_frags import read_fragdata
import joblib
from joblib import Memory
import multiprocessing as mp
import os
import glob
import operator
import numexpr as ne
from pytritex.utils import return_size, parse_size
import logging
from dask.distributed import Client
from pytritex.sequencing_coverage.add_hic_cov import add_hic_cov
from pytritex.scaffold_10x import print_agp
import time
import pandas as pd
from pytritex.scaffold_hic.add_psmol_fpairs import add_psmol_fpairs
from pytritex.scaffold_hic.hic_cov_psmol import hic_cov_psmol
import json
logger = logging.getLogger("distributed.comm.tcp")
logger.setLevel(logging.ERROR)
import numpy as np
import hashlib
import itertools


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--processes", dest="procs", default=mp.cpu_count(), type=int)
    parser.add_argument("-dc", "--dask-cache", default="dask_data", type=str)
    parser.add_argument("-um", "--use-molecules", dest="use_mols", action="store_true", default=False)
    parser.add_argument("--min-length", nargs=3, type=int, default=[int(1e5), int(5e5), int(1e5)],
                        help="Minimum, maximum, and step for minimum length of scaffolds/super-scaffolds to try to combine using HiC.")
    parser.add_argument("--min-nlinks", nargs=3, type=int, default=[40, 100, 10],
                        help="Minimum, maximum, and step for minimum no. of links between scaffolds/super-scaffolds to try to combine using HiC.")
    parser.add_argument("--max-cM-dist", nargs=3, type=float, default=[5, 30, 5],
                        help="Minimum, maximum, and step for maximum distance in cM length of scaffolds/super-scaffolds to try to combine using HiC.")
    parser.add_argument("--min-nfrags", nargs=3, type=float, default=[20, 50, 10],
                        help="Minimum, maximum, and step for maximum distance in cM length of scaffolds/super-scaffolds to try to combine using HiC.")
    parser.add_argument("--mem", default="20GB", type=str)
    parser.add_argument("--save", default=False, action="store_true")
    parser.add_argument("--dir", default=None)
    parser.add_argument("--hash", default=None)
    parser.add_argument("--table", default=None)
    parser.add_argument("save_prefix", default="assembly",
                        help="Folder with the previous run (up to scaffold_10x, included)")
    parser.add_argument("fragments_bed")
    parser.add_argument("fragments")
    args = parser.parse_args()

    os.makedirs(args.dask_cache, exist_ok=True)
    dask.config.global_config.update({"temporary-directory": args.dask_cache})
    # dask.config.global_config.update({"distributed.comm.timeouts.tcp": "300s"})
    ne.set_num_threads(args.procs)
    worker_mem = return_size(parse_size(args.mem)[0] / 1, "GB")    

    # This is a STUPID hack, but I should have thought of dumping the results in a different way earlier ..
    assembly = sorted([(_, os.stat(_).st_mtime) for _ in
            glob.glob(f"{args.save_prefix}/joblib/pytritex/chimera_breaking/break_10x/break_10x/*/output.pkl")],
           key=operator.itemgetter(1), reverse=True)[0][0]
    assembly = joblib.load(assembly)["assembly"]
    if args.hash is not None:
        best_dir = [os.path.dirname(_) for _ in
            glob.glob(f"{args.save_prefix}/joblib/pytritex/scaffold_10x/{args.hash}*/orientation/res/_metadata")][0]
    elif args.table is not None:
        table = pd.read_csv(args.table, delimiter="\t")
        table = table.sort_values(["n50"], ascending=[False])
        best_row = table[table["n50"] == table["n50"].max()]
        mhash = best_row["hash"].values[0]
        best_dir = [os.path.dirname(_) for _ in
            glob.glob(f"{args.save_prefix}/joblib/pytritex/scaffold_10x/{mhash}*/orientation/res/_metadata")][0]
    elif args.dir is not None:
        best_dir = [os.path.dirname(_) for _ in
            glob.glob(f"{args.dir}")][0]
    else:
        best_dir = None

    client = Client(set_as_default=True, timeout=60, direct_to_workers=False, memory_limit=worker_mem,
                    nanny=True, address=None)
    client.cluster.scale(args.procs)
    # DO IT TWICE, sometimes the first time is not enough.
    client.cluster.scale(args.procs)
    memory = Memory(os.path.join(".", args.save_prefix), compress=("zlib", 6), verbose=1)
        
    if best_dir is not None and not os.path.exists(best_dir):
        best_dir = None
    if best_dir is not None:
        best_result = os.path.dirname(best_dir)
    else:
        dirs = [os.path.dirname(_) for _ in
        glob.glob(f"{args.save_prefix}/joblib/pytritex/scaffold_10x/*/orientation/res/_metadata")]
        results = dict()
        for directory in dirs:
            result = dd.read_parquet(directory)
            results[os.path.dirname(directory)] = n50(result.length.values.compute())
        best_result = sorted(results.items(), key=operator.itemgetter(1), reverse=True)[0][0]
        
    map_10x = {"membership": os.path.join(best_result, "membership"),
               "result": os.path.join(best_result, "res")}
    save_dir = os.path.join(args.save_prefix, "joblib", "pytritex", "hic_map")
    assembly_10x = memory.cache(init_10x_assembly)(
        assembly=assembly, map_10x=map_10x, gap_size=200, molecules=args.use_mols, save=save_dir)

    # Print the AGP
    out_agp = os.path.join(save_dir, "scaffold10x.agp")
    print_agp(assembly_10x["agp"], out_agp)
    
    # Now anchor the scaffolds
    assembly_10x = memory.cache(anchor_scaffolds, ignore=["client"])(
        assembly_10x, save_dir, species="wheat", client=client)
    assembly_10x = memory.cache(add_hic_cov)(assembly_10x, save_dir=save_dir,
                                             binsize=1e4, binsize2=1e6, minNbin=100, innerDist=3e5)
    # fai, fragfile, map_10x, savedir=None
    fragment_data = memory.cache(read_fragdata)(fai=assembly["fai"],
                                                fragfile=args.fragments_bed,
                                                map_10x=assembly_10x,
                                                savedir=save_dir)
    # frag_data$info[!is.na(hic_chr) & length >= 3e5, .(scaffold, nfrag, chr=hic_chr, cM=popseq_cM)]->hic_info
    # # make Hi-C map
    # hic_map(info=hic_info, assembly=assembly_v2, frags=frag_data$bed, species="wheat", ncores=21,
    # 	min_nfrag_scaffold=50, max_cM_dist = 50,
    # 	binsize=1e5, min_nfrag_bin=20, gap_size=100)->hic_map_v1

    assembly_10x = memory.cache(calculate_hic_link_weights)(assembly_10x, save_dir)
    for cM, min_nlinks, min_length, min_nfrags in itertools.product(
        np.arange(args.max_cM_dist[0], args.max_cM_dist[1] + args.max_cM_dist[2], args.max_cM_dist[2], dtype=float),
        np.arange(args.min_nlinks[0], args.min_nlinks[1] + args.min_nlinks[2], args.min_nlinks[2], dtype=int),
        np.arange(args.min_length[0], args.min_length[1] + args.min_length[2], args.min_length[2], dtype=int),
        np.arange(args.min_nfrags[0], args.min_nfrags[1] + args.min_nfrags[2], args.min_nfrags[2], dtype=int),
    ):
        sha = hashlib.sha256()
        params = (cM, min_nlinks, min_length)
        sha.update(str(params).encode())
        hash_string = sha.hexdigest()
        hic_save_dir = os.path.join(save_dir, hash_string)
        os.makedirs(hic_save_dir, exist_ok=True)
        params_dict = {"max_cM_dist": int(cM), "min_length": int(min_length), "min_nlinks": int(min_nlinks)}
        with open(os.path.join(hic_save_dir, "meta.json"), "wt") as meta_out:
            json.dump(params_dict, meta_out)
        try:
            memory.cache(hic_map, ignore=["client"])(
                assembly=assembly_10x, client=client, fragment_data=fragment_data, species="wheat",
                ncores=args.procs, min_length=min_length, save_dir=hic_save_dir, max_cM_dist=cM,
                min_nfrag_scaffold=min_nfrags, min_nlinks=min_nlinks)
        except KeyboardInterrupt:
            raise
        except Exception as exc:
            logger.error("Run {} with following parameters FAILED: {}".format(hash_string[:10], params_dict))
            logger.error(exc)
            raise exc

    # Now select the best
    best_hic = None
    # something = add_psmol_fpairs(assembly=assembly_10x, hic_map=best_hic, nucfile=args.fragments)

    # best_hic["hic_1Mb"] = bin_hic_step(hic=best_hic["links"], frags=best_hic["frags"],
    #                                    binsize=10**6, chrlen=best_hic["chrlen"], chrs=np.arange(1, 22, 1, dtype=int))
    # best_hic["hic_1Mb"]["norm"] = normalize_cis()

    # bin_hic_step(hic=hic_map_v1$links, frags=hic_map_v1$frags, binsize=1e6,
    # 	    chrlen=hic_map_v1$chrlen, chrs=1:21, cores=21)->hic_map_v1$hic_1Mb
    # normalize_cis(hic_map_v1$hic_1Mb, ncores=21, percentile=0, omit_smallest=1)->hic_map_v1$hic_1Mb$norm
    #
    # find_inversions(hic_map=hic_map_v1, links=hic_map_v1$hic_1Mb$norm, species="wheat", cores=30)->inv
    #
    # hic_cov_psmol(hic_map=hic_map_v1, binsize=1e3, binsize2=1e6, maxdist=1e6, cores=40) -> cov


    # f <- 'Triticum_aestivum_Claire_EIv1.1_DpnII_fragments_30bp_split.nuc.txt'
    # add_psmol_fpairs(assembly=assembly_v1, hic_map=hic_map_v1, map_10x=assembly_v1_10x,
    # 		 assembly_10x=assembly_v2, nucfile=f)->hic_map_v1


    # bin_hic_step(hic=hic_map_v1$links, frags = hic_map_v1$frags, binsize = 1e6,
    #                                                                        chrlen = hic_map_v1$chrlen, chrs = 1:21, cores = 21)->hic_map_v1$hic_1Mb

    # normalize_cis(hic_map_v1$hic_1Mb, ncores = 21, percentile = 0, omit_smallest = 1)->hic_map_v1$hic_1Mb$norm

    return
    # # exclude scaffolds <= 300 kb from Hi-C map construction
    # frag_data$info[!is.na(hic_chr) & length >= 3e5, .(scaffold, nfrag, chr=hic_chr, cM=popseq_cM)]->hic_info
    #
    # # make Hi-C map
    # hic_map(info=hic_info, assembly=assembly_v2, frags=frag_data$bed, species="wheat", ncores=21,
    # 	min_nfrag_scaffold=50, max_cM_dist = 50,
    # 	binsize=1e5, min_nfrag_bin=20, gap_size=100)->hic_map_v1
    #
    # # add links data and compute contact matrices
    # f <- 'Triticum_aestivum_Claire_EIv1.1_DpnII_fragments_30bp_split.nuc.txt'
    # add_psmol_fpairs(assembly=assembly_v1, hic_map=hic_map_v1, map_10x=assembly_v1_10x,
    # 		 assembly_10x=assembly_v2, nucfile=f)->hic_map_v1

if __name__ == "__main__":
    mp.freeze_support()
    main()
