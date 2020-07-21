import pandas as pd
import argparse
import multiprocessing as mp
import os
import dask
import matplotlib
matplotlib.use("agg")  # Avoid pesky Gdk error messages
dask.config.set({'distributed.worker.multiprocessing-method': 'spawn'})
from pytritex.initialisation import initial
from pytritex.anchoring.anchor_scaffolds import anchor_scaffolds
from pytritex.sequencing_coverage.add_molecule_cov import add_molecule_cov
from pytritex.sequencing_coverage.add_hic_cov import add_hic_cov
# from pytritex.chimera_breaking.find_10x_breaks import find_10x_breaks
from pytritex.chimera_breaking.break_10x import break_10x
import itertools
from pytritex.scaffold_10x.__init__ import scaffold_10x
from pytritex.utils import n50
from joblib import Memory, dump, load
import numexpr as ne
import dask
import dask.dataframe as dd
import logging
from dask.distributed import Client, SpecCluster
# from dask.distributed import Scheduler, Worker, Nanny
from pytritex.utils import return_size, parse_size
import logging
import joblib
import time
logger = logging.getLogger("distributed.comm.tcp")
logger.setLevel(logging.ERROR)


def dispatcher(assembly, save_dir, memory, client, row, ncores):
    membership, result = scaffold_10x(assembly,
                          memory=memory,
                          save_dir=save_dir,
                          client=client,
                          min_npairs=row.npairs,
                          max_dist=row.dist, popseq_dist=5, max_dist_orientation=5,
                          min_nsample=row.nsample,
                          min_nmol=row.nmol, unanchored=False, ncores=ncores)

    logger.warning("""{ctime} Parameters: {row}\n
Result: {res}\n""".format(ctime=time.ctime(),
                          row=row, res=n50(dd.read_parquet(result)["length"].values.compute())
                          ))
    return {"membership": membership, "result": result, "info": assembly["info"],
            "row": row}


def grid_evaluation(assembly, args, client, memory):

    npairs = list(range(args.npairs[0], args.npairs[1] + 1))
    nmol = list(range(args.nmol[0], args.nmol[1] + 1))
    nsamples = list(range(args.nsamples[0], args.nsamples[1] + 1))
    dist = list(range(args.dist[0], args.dist[1] + args.dist[2],
                      args.dist[2]))
    grid = pd.DataFrame(dict(zip(("npairs", "nmol", "nsample", "dist"),
                                 list(zip(*itertools.product(npairs,
                                                             nmol,
                                                             nsamples,
                                                             dist))))))
    print("Starting grid evaluation")
    results = []
    for _index, row in grid.iterrows():
        worker_mem = return_size(parse_size(args.mem)[0] / 1, "GB")
        # DO IT TWICE, sometimes the first time is not enough.
        client.cluster.scale(args.procs)
        result = dispatcher(
            assembly, save_dir=args.save_prefix,
            memory=memory, row=row, client=client, ncores=args.procs)
        results.append(result)
        logger.warning("%s Finished row %s (%s)", time.ctime(), _index, row)
    return results


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--processes", dest="procs", default=mp.cpu_count(), type=int)
    parser.add_argument("-s", "--save-prefix", default="assembly")
    parser.add_argument("-dc", "--dask-cache", default="dask_data", type=str)
    parser.add_argument("--nsamples", nargs=2, type=int, default=[1, 2],
                        help="Inclusive range to investigate for the min_nsample parameter. Default: [2, 3]")
    parser.add_argument("--nmol", nargs=2, type=int, default=[2, 3],
                        help="Inclusive range to investigate for the min_nmol parameter. Default: [2, 3]")
    parser.add_argument("--npairs", nargs=2, type=int, default=[2, 3],
                        help="Inclusive range to investigate for the min_npairs parameter. Default: [2, 3]")
    parser.add_argument("--dist", nargs=3, type=int, default=[6 * 10**4, 9 * 10**4, 10**4],
                        help="Inclusive range to investigate for the dist parameter, with step size.\
                        Default: [60,000, 90,000, 10,000]")
    # parser.add_argument("-umfs", "--use-memory_fs", default=False, action="store_true")
    parser.add_argument("--mem", default="20GB", type=str)
    parser.add_argument("--save", default=False, action="store_true")
    parser.add_argument("--10x", dest="tenx")
    parser.add_argument("popseq")
    parser.add_argument("css")
    parser.add_argument("hic")
    parser.add_argument("fasta")
    mp.set_start_method("spawn", force=True)
    args = parser.parse_args()

    # Initial set-up
    os.makedirs(args.dask_cache, exist_ok=True)
    dask.config.global_config.update({"temporary-directory": args.dask_cache})
    # dask.config.global_config.update({"distributed.comm.timeouts.tcp": "300s"})
    ne.set_num_threads(args.procs)
    worker_mem = return_size(parse_size(args.mem)[0] / 1, "GB")
    client = Client(set_as_default=True, timeout=60, direct_to_workers=False, memory_limit=worker_mem,
                    nanny=False, dashboard_address=None)
    client.cluster.scale(args.procs)
    # DO IT TWICE, sometimes the first time is not enough.
    client.cluster.scale(args.procs)

    # Initial assembly
    memory = Memory(os.path.join(".", args.save_prefix), compress=("zlib", 6), verbose=1)
    assembly = memory.cache(initial, ignore=["cores", "client", "memory", "ram"])(
        args.popseq, args.fasta, args.css, args.tenx, args.hic, args.save_prefix, client=client, memory=memory,
        ram=args.mem)
    assembly = memory.cache(anchor_scaffolds, ignore=["client"])(
        assembly, os.path.join(args.save_prefix, "joblib", "pytritex", "anchoring"), species="wheat", client=client)
    cov_base = os.path.join(args.save_prefix, "joblib", "pytritex", "sequencing_coverage")
    assembly = memory.cache(add_molecule_cov)(assembly, binsize=200, save_dir=cov_base)
    assembly = memory.cache(add_hic_cov)(assembly, save_dir=cov_base,
                                         binsize=5e3, binsize2=5e4, minNbin=50, innerDist=3e5)
    assembly_v1 = memory.cache(break_10x, ignore=["cores", "client"])(
        assembly, client=client,
        ratio=-3, interval=5e4, minNbin=20, dist=2e3, save_dir=args.save_prefix,
        slop=2e2, species="wheat", intermediate=False, cores=args.procs)["assembly"]
    results = grid_evaluation(assembly_v1, args, client=client, memory=memory)
    results_name = os.path.join(args.save_prefix, "joblib", "pytritex", "scaffold_10x", "results.pkl")
    joblib.dump(results, results_name, compress=("gzip", 6))
    client.close()
    print("Broken chimeras")
    return


if __name__ == "__main__":
    mp.freeze_support()
    main()
