import pandas as pd
import argparse
import multiprocessing as mp
import os
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
# import dask.dataframe as dd
import logging
from dask.distributed import Client, LocalCluster



def dispatcher(assembly, row):
    result = scaffold_10x(assembly,
                          prefix="scaffold_10x", min_npairs=row.npairs,
                          max_dist=row.dist, popseq_dist=5, max_dist_orientation=5,
                          min_nsample=row.nsample,
                          min_nmol=row.nmol, unanchored=False, ncores=1)
    print("""Parameters: {row}\n
Result: {res}\n""".format(row=row, res=n50(result["info"]["length"])))
    return result


def grid_evaluation(assembly, args):

    grid = pd.DataFrame(dict(zip(("npairs", "nmol", "nsample", "dist"),
                                 list(zip(*itertools.product((2, 3),
                                                             (2, 3),
                                                             (1, 2),
                                                             range(6 * 10**4, 10**5, 10**4)))))))
    print("Starting grid evaluation")
    pool = mp.Pool(processes=args.procs)
    # results = pool.starmap(dispatcher, [(assembly, row) for index, row in grid.iterrows()])
    _index, row = next(grid.iterrows())
    result = dispatcher(assembly, row)
    return result


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--processes", dest="procs", default=mp.cpu_count(), type=int)
    parser.add_argument("-s", "--save-prefix", default="assembly")
    parser.add_argument("-dc", "--dask-cache", default="dask_data", type=str)
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
    dask.config.global_config.update({"distributed.comm.timeouts.tcp": "300s"})
    ne.set_num_threads(args.procs)
    cluster = LocalCluster(n_workers=args.procs, threads_per_worker=1, processes=True,
                           scheduler_port=0, protocol="tcp://",
                           silence_logs=logging.WARNING)
    cluster.adapt(maximum_memory=args.mem)
    client = Client(cluster, set_as_default=True)

    # Initial assembly
    memory = Memory(os.path.join(".", args.save_prefix), compress=("zlib", 6), verbose=10)
    assembly = memory.cache(initial, ignore=["cores", "client"])(
        args.popseq, args.fasta, args.css, args.tenx, args.hic, args.save_prefix, client=client)
    assembly = memory.cache(anchor_scaffolds, ignore=["client"])(
        assembly, args.save_prefix, client=client, species="wheat")
    assembly = memory.cache(add_molecule_cov, ignore=["cores", "client"])(
        assembly, cores=args.procs, client=client, binsize=200, save_dir=args.save_prefix)
    assembly = memory.cache(add_hic_cov, ignore=["cores", "client"])(
        assembly, save_dir=args.save_prefix, client=client,
        cores=args.procs, binsize=5e3, binsize2=5e4, minNbin=50, innerDist=3e5)
    return
    assembly_v1 = memory.cache(break_10x, ignore=["cores"])(
        assembly,
        ratio=-3, interval=5e4, minNbin=20, dist=2e3,
        slop=2e2, species="wheat", intermediate=False, cores=args.procs)["assembly"]
    return
    print("Broken chimeras")
    grid_evaluation(assembly_v1, args)
    return


if __name__ == "__main__":
    mp.freeze_support()
    main()
