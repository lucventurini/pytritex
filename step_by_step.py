import pandas as pd
import argparse
import os
import pysam
import pytritex.utils.read_cssaln
import multiprocessing as mp
from pytritex.utils.read_10x import read_10x_molecules
from pytritex.init_assembly import init_assembly
import subprocess as sp
import pickle
import gzip
from time import ctime
from pytritex.anchor_scaffolds import anchor_scaffolds
from pytritex.add_molecule_cov import add_molecule_cov
from pytritex.add_hic_cov import add_hic_cov
from pytritex.find_10x_breaks import find_10x_breaks
from pytritex.break_10x import break_10x
import io
import itertools
from pytritex.scaffold_10x import scaffold_10x


def initial(args, popseq):

    if not os.path.exists(args.fasta + ".fai"):
        pysam.Fastafile(args.fasta)
    fai = pd.read_csv(args.fasta + ".fai", sep="\t",
                      usecols=[0, 1], names=["scaffold", "length"])
    # Alignment of genetic markers used for the genetic map. In this example, the Morex WGS assembly by IBSC (2012).
    cssaln = pytritex.utils.read_cssaln.read_morexaln_minimap(
        paf=args.css, popseq=popseq, minqual=30, minlen=500, ref=True
    )

    # Read the list of Hi-C links.
    fpairs_command = 'find {} -type f | grep "_fragment_pairs.tsv.gz$"'.format(args.hic)
    fpairs = pd.concat(
        [
            pd.read_csv(
                fname.decode().rstrip(), sep="\t", header=None, names=["scaffold1", "pos1", "scaffold2", "pos2"]
            )
            for fname in sp.Popen(fpairs_command, shell=True, stdout=sp.PIPE).stdout
        ]
    )

    tenx_command = 'find {} -type f | grep "molecules.tsv.gz$"'.format(args.tenx)
    tenx_files = [line.rstrip().decode() for line in
                  sp.Popen(tenx_command, shell=True, stdout=sp.PIPE).stdout]
    # names(f) <- c("S1", "S2")
    if tenx_files:
        tenx_dict = []
        for fname in tenx_files:
            sample = os.path.basename(os.path.dirname(fname))
            tenx_dict.append((sample, fname))
        molecules = read_10x_molecules(tenx_dict, ncores=args.procs)
    else:
        print("No 10X molecules! Error")
        print("Command:", tenx_command)
        import sys
        sys.exit(1)

    assembly = init_assembly(fai=fai, cssaln=cssaln, fpairs=fpairs, molecules=molecules, rename=None)
    return assembly


def dispatcher(assembly, row):
    result = scaffold_10x(assembly, **{"prefix": "scaffold_10x", "min_npairs": row.npairs, "max_dist": row.dist,
              "popseq_dist": 5, "max_dist_orientation": 5, "min_nsample": row.nsample, "min_mol": row.nmol,
              "unanchored": False, "ncores": 1})
    print("""Parameters: {row}\n
Result: {res}""".format(row=row, res=n50(result["info"]["length"])))
    print("N50")
    return result


def grid_evaluation(assembly, args):

    grid = pd.DataFrame(dict(zip(("npairs", "nmol", "nsample", "dist"),
                                 list(zip(*itertools.product((2, 3),
                                                             (2, 3),
                                                             (1, 2),
                                                             range(6 * 10**4, 10**5, 10**4)))))))
    print("Starting grid evaluation")
    pool = mp.Pool(processes=args.procs)
    results = pool.starmap(dispatcher, [(assembly, row) for index, row in grid.iterrows()])

    return

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--processes", dest="procs", default=mp.cpu_count(), type=int)
    parser.add_argument("-s", "--save-prefix", default="assembly")
    parser.add_argument("--10x", dest="tenx")
    parser.add_argument("popseq")
    parser.add_argument("css")
    parser.add_argument("hic")
    parser.add_argument("fasta")
    args = parser.parse_args()
    popseq = pd.read_pickle(args.popseq)
    popseq.columns = popseq.columns.str.replace("morex", "css")
    assembly = initial(args, popseq)
    with io.BufferedWriter(gzip.open(args.save_prefix + ".initial_assembly.pickle.gz", "wb")) as dump:
        pickle.dump(assembly, dump)
    assembly = anchor_scaffolds(assembly, popseq=popseq, species="wheat")
    assembly = add_molecule_cov(assembly, cores=args.procs, binsize=200)
    assembly = add_hic_cov(assembly, cores=args.procs, binsize=5e3, binsize2=5e4, minNbin=50, innerDist=3e5)
    breaks = find_10x_breaks(assembly)
    b = breaks[breaks["d"] >= 1e4].sort_values("d", ascending=False).head(100)
    print(ctime(), "Started to save the data")
    with io.BufferedWriter(gzip.open(args.save_prefix + ".anchored_assembly.pickle.gz",
                                     "wb", compresslevel=1)) as dump:
        pickle.dump(assembly, dump)
    with io.BufferedWriter(gzip.open(args.save_prefix + ".initial_breaks.pickle.gz",
                                     "wb", compresslevel=1)) as dump:
        pickle.dump(breaks, dump)
    print(ctime(), "Finished saving the data")
    a = break_10x(
        assembly, prefix="scaffold_corrected", ratio=-3,
        interval=5e4, minNbin=20, dist=2e3, slop=2e2, species="wheat", intermediate=False, ncores=8)
    print("Broken chimeras")
    assembly_v1 = a["assembly"]
    with io.BufferedWriter(gzip.open(args.save_prefix + ".assembly_v1.pickle.gz",
                                     "wb", compresslevel=1)) as dump:
        pickle.dump(assembly_v1, dump)
    with io.BufferedWriter(gzip.open(args.save_prefix + ".breaks.pickle.gz",
                                     "wb", compresslevel=1)) as dump:
        pickle.dump(a["breaks"], dump)
    print("Starting grid evaluation")



    return


main()
