import pysam
import pandas as pd
from ..utils.read_cssaln import read_morexaln_minimap
import os
import subprocess as sp
from .read_10x import read_10x_molecules
from .init_assembly import init_assembly


def initial(args, popseq):

    if not os.path.exists(args.fasta + ".fai"):
        pysam.Fastafile(args.fasta)
    fai = pd.read_csv(args.fasta + ".fai", sep="\t",
                      usecols=[0, 1], names=["scaffold", "length"]).reset_index(drop=False).rename(
        columns={"index": "scaffold_index"}
    )

    # Alignment of genetic markers used for the genetic map. In this example, the Morex WGS assembly by IBSC (2012).
    cssaln = read_morexaln_minimap(paf=args.css, popseq=popseq, fai=fai, minqual=30, minlen=500, ref=True)

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
    # Now let us change the scaffold1 and scaffold2
    fpairs = fai.loc[["scaffold", "scaffold_index"]].rename(columns={"scaffold": "scaffold1",
                                                                     "scaffold_index": "scaffold_index1"}).merge(
        fpairs, on="scaffold1").drop("scaffold1", axis=1)
    fpairs = fai.loc[["scaffold", "scaffold_index"]].rename(columns={"scaffold": "scaffold2",
                                                                     "scaffold_index": "scaffold_index2"}).merge(
        fpairs, on="scaffold2").drop("scaffold2", axis=1)

    tenx_command = 'find {} -type f | grep "molecules.tsv.gz$"'.format(args.tenx)
    tenx_files = [line.rstrip().decode() for line in
                  sp.Popen(tenx_command, shell=True, stdout=sp.PIPE).stdout]
    # names(f) <- c("S1", "S2")
    if tenx_files:
        tenx_dict = []
        for fname in tenx_files:
            sample = os.path.basename(os.path.dirname(fname))
            tenx_dict.append((sample, fname))
        molecules = read_10x_molecules(tenx_dict, fai, ncores=args.procs)
    else:
        print("No 10X molecules! Error")
        print("Command:", tenx_command)
        import sys
        sys.exit(1)

    assembly = init_assembly(fai=fai, cssaln=cssaln, fpairs=fpairs, molecules=molecules, rename=None)
    return assembly
