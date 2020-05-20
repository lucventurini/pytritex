import pysam
import pandas as pd
from pytritex.initialisation.read_cssaln import read_morexaln_minimap
import os
import subprocess as sp
from .read_10x import read_10x_molecules
import numpy as np


def initial(args, popseq):

    if not os.path.exists(args.fasta + ".fai"):
        pysam.Fastafile(args.fasta)
    fai = pd.read_csv(args.fasta + ".fai", sep="\t",
                      usecols=[0, 1], names=["scaffold", "length"]).reset_index(drop=False).rename(
        columns={"index": "scaffold_index"}
    )

    fai.loc[:, "scaffold_index"] = pd.to_numeric(fai["scaffold_index"] + int(1), downcast="signed")
    fai.loc[:, "length"] = pd.to_numeric(fai["length"], downcast="signed")
    fai.loc[:, "orig_scaffold_index"] = fai["scaffold_index"]
    fai.loc[:, "start"] = fai.loc[:, "orig_start"] = np.array([1], dtype=fai["length"].dtype)
    fai.loc[:, "end"] = fai.loc[:, "orig_end"] = fai["length"]

    # fai = fai.convert_dtypes()
    # Alignment of genetic markers used for the genetic map. In this example, the Morex WGS assembly by IBSC (2012).
    cssaln = read_morexaln_minimap(paf=args.css, popseq=popseq, fai=fai, minqual=30, minlen=500, ref=True)
    # cssaln = cssaln.convert_dtypes()

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
    fpairs = fai[["scaffold", "scaffold_index"]].rename(columns={"scaffold": "scaffold1",
                                                                     "scaffold_index": "scaffold_index1"}).merge(
        fpairs, on="scaffold1").drop("scaffold1", axis=1)
    fpairs = fai[["scaffold", "scaffold_index"]].rename(columns={"scaffold": "scaffold2",
                                                                     "scaffold_index": "scaffold_index2"}).merge(
        fpairs, on="scaffold2").drop("scaffold2", axis=1)
    fpairs.loc[:, "orig_scaffold_index1"] = fpairs["scaffold_index1"]
    fpairs.loc[:, "orig_scaffold_index2"] = fpairs["scaffold_index2"]
    fpairs.loc[:, "pos1"] = pd.to_numeric(fpairs["pos1"], downcast="signed")
    fpairs.loc[:, "pos2"] = pd.to_numeric(fpairs["pos2"], downcast="signed")
    fpairs.loc[:, "orig_pos1"] = fpairs["pos1"]
    fpairs.loc[:, "orig_pos2"] = fpairs["pos2"]

    # fpairs = fpairs.convert_dtypes()
    tenx_command = 'find {} -type f | grep "molecules.tsv.gz$"'.format(args.tenx)
    tenx_files = [line.rstrip().decode() for line in
                  sp.Popen(tenx_command, shell=True, stdout=sp.PIPE).stdout]
    if tenx_files:
        tenx_dict = []
        for index, fname in enumerate(tenx_files):
            sample = os.path.basename(os.path.dirname(fname))
            tenx_dict.append((index, sample, fname))
        samples = list(zip(*tenx_dict))
        samples = pd.DataFrame({"index": samples[0], "sample": samples[1], "fname": samples[2]})
        molecules, barcodes = read_10x_molecules(samples, fai, ncores=args.procs)
    else:
        print("No 10X molecules! Error")
        print("Command:", tenx_command)
        import sys
        sys.exit(1)

    # molecules = molecules.convert_dtypes()
    assembly = {"popseq": popseq,
                "fai": fai,
                "cssaln": cssaln,
                "10x_samples": samples,
                "fpairs": fpairs,
                "molecules": molecules,
                "barcodes": barcodes}

    return assembly
