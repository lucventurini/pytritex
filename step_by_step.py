import pandas as pd
import argparse
import multiprocessing as mp
import pickle
from time import ctime
from pytritex.anchoring.anchor_scaffolds import anchor_scaffolds
from pytritex.add_molecule_cov import add_molecule_cov
from pytritex.add_hic_cov import add_hic_cov
from pytritex.find_10x_breaks import find_10x_breaks
from pytritex.break_10x import break_10x
import io
import itertools
from pytritex.scaffold_10x import scaffold_10x
from pytritex.utils import n50
from pytritex.initialisation import initial
import subprocess as sp


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
    parser.add_argument("--save", default=False, action="store_true")
    parser.add_argument("--10x", dest="tenx")
    parser.add_argument("popseq")
    parser.add_argument("css")
    parser.add_argument("hic")
    parser.add_argument("fasta")
    args = parser.parse_args()
    popseq = pd.read_pickle(args.popseq)
    popseq.columns = popseq.columns.str.replace("morex", "css")
    assembly = initial(args, popseq)
    initial_fname = args.save_prefix + ".initial_assembly.pickle"
    with io.BufferedWriter(open(initial_fname, "wb")) as dump:
        pickle.dump(assembly, dump)
    sp.call(["gzip", "-f", initial_fname])
    assembly = anchor_scaffolds(assembly, popseq=popseq, species="wheat")
    assembly = add_molecule_cov(assembly, cores=args.procs, binsize=200)
    assembly = add_hic_cov(assembly, cores=args.procs, binsize=5e3, binsize2=5e4, minNbin=50, innerDist=3e5)
    breaks = find_10x_breaks(assembly["molecule_cov"])
    b = breaks[breaks["d"] >= 1e4].sort_values("d", ascending=False).head(100)
    if args.save is True:
        print(ctime(), "Started to save the data")
        fname = args.save_prefix + ".anchored_assembly.pickle"
        with io.BufferedWriter(open(fname, "wb")) as dump:
            pickle.dump(assembly, dump)
        sp.call(["gzip", "-f", fname])
        fname = args.save_prefix + ".initial_breaks.pickle"
        with io.BufferedWriter(open(fname, "wb")) as dump:
            pickle.dump(breaks, dump)
        sp.call(["gzip", "-f", fname])
        print(ctime(), "Finished saving the data")
    return
    a = break_10x(
        assembly, prefix="scaffold_corrected", ratio=-3,
        interval=5e4, minNbin=20, dist=2e3, slop=2e2, species="wheat", intermediate=False, ncores=8)
    print("Broken chimeras")
    assembly_v1 = a["assembly"]
    fname = args.save_prefix + ".assembly_v1.pickle",
    with io.BufferedWriter(open(fname, "wb")) as dump:
        pickle.dump(assembly_v1, dump)
    sp.call(["gzip", "-f", fname])
    fname = args.save_prefix + ".breaks.pickle.gz"
    with io.BufferedWriter(open(fname, "wb")) as dump:
        pickle.dump(a["breaks"], dump)
    sp.call(["gzip", "-f", fname])
    # grid_evaluation(assembly, args)

    return


main()
