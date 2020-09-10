import argparse
import matplotlib
matplotlib.use("agg")  # Avoid pesky Gdk error messages
import dask
dask.config.set({'distributed.worker.multiprocessing-method': 'spawn'})
import dask.dataframe as dd
from pytritex.scaffold_hic.init_10x_assembly import init_10x_assembly
from pytritex.anchoring import anchor_scaffolds
from pytritex.utils import n50
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
import time
logger = logging.getLogger("distributed.comm.tcp")
logger.setLevel(logging.ERROR)


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--processes", dest="procs", default=mp.cpu_count(), type=int)
    parser.add_argument("-dc", "--dask-cache", default="dask_data", type=str)
    parser.add_argument("-um", "--use-molecules", dest="use_mols", action="store_true", default=False)
    parser.add_argument("--mem", default="20GB", type=str)
    parser.add_argument("--save", default=False, action="store_true")
    parser.add_argument("save-prefix", dest="save_prefix", default="assembly",
                        help="Folder with the previous run (up to scaffold_10x, included)")
    parser.add_argument("fragments_bed")
    parser.add_argument("fragments")
    args = parser.parse_args()

    os.makedirs(args.dask_cache, exist_ok=True)
    dask.config.global_config.update({"temporary-directory": args.dask_cache})
    # dask.config.global_config.update({"distributed.comm.timeouts.tcp": "300s"})
    ne.set_num_threads(args.procs)
    worker_mem = return_size(parse_size(args.mem)[0] / 1, "GB")
    client = Client(set_as_default=True, timeout=60, direct_to_workers=False, memory_limit=worker_mem,
                    nanny=True, address=None)
    client.cluster.scale(args.procs)
    # DO IT TWICE, sometimes the first time is not enough.
    client.cluster.scale(args.procs)
    memory = Memory(os.path.join(".", args.save_prefix), compress=("zlib", 6), verbose=1)

    # This is a STUPID hack, but I should have thought of dumping the results in a different way earlier ..
    assembly = sorted([(_, os.stat(_).st_mtime) for _ in
            glob.glob(f"{args.save_prefix}/joblib/pytritex/chimera_breaking/break_10x/break_10x/*/output.pkl")],
           key=operator.itemgetter(1), reverse=True)
    assembly = joblib.load(assembly)["assembly"]
    dirs = [os.path.dirname(_) for _ in
            glob.glob(f"{args.save_prefix}/joblib/pytritex/scaffold_10x/*/orientation/res/_metadata")]
    results = dict()
    for directory in dirs:
        result = dd.read_parquet(directory)
        results[os.path.dirname(directory)] = n50(result.length.values.compute())
    best_result = sorted(results.items(), key=operator.itemgetter(1), reverse=True)[0][0]
    map_10x = {"membership": os.path.join(best_result, "membership"),
               "result": os.path.join(best_result, "result")}
    save_dir = os.path.join(args.save_prefix, "joblib", "pytritex", "hic_map")
    assembly_10x = memory.cache(init_10x_assembly)(
        assembly=assembly, map_10x=map_10x, gap_size=200, molecules=args.use_mols, save=save_dir)

    # Now anchor the scaffolds
    assembly_10x = memory.cache(anchor_scaffolds, ignore=["client"])(
        assembly_10x, save_dir, species="wheat", client=client)
    assembly_10x = memory.cache(add_hic_cov)(assembly, save_dir=save_dir,
                                             binsize=1e4, binsize2=1e6, minNbin=100, innerDist=3e5)
    