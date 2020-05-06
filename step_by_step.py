import pandas as pd
import argparse
import os
import pysam
import pytritex.utils.read_cssaln
import multiprocessing as mp
from pytritex.utils.read_10x import read_10x_molecules
import subprocess as sp

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--processes", dest="procs", default=mp.cpu_count(), type=int)
    parser.add_argument("--10x", dest="tenx")
    parser.add_argument("popseq")
    parser.add_argument("css")
    parser.add_argument("hic")
    parser.add_argument("fasta")
    args = parser.parse_args()

    popseq = pd.read_pickle(args.popseq)
    popseq.columns = popseq.columns.str.replace("morex", "css")
    if not os.path.exists(args.fasta + ".fai"):
        pysam.Fastafile(args.fasta)
    fai = pd.read_csv(args.fasta + ".fai", sep="\t",
                      usecols=[0, 1], names=["scaffold", "length"])
    # Alignment of genetic markers used for the genetic map. In this example, the Morex WGS assembly by IBSC (2012).
    morexaln = pytritex.utils.read_cssaln.read_morexaln_minimap(
        paf=args.css, popseq=popseq, minqual=30, minlen=500
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
    return


main()