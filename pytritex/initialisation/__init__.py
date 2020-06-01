import pysam
import pandas as pd
from pytritex.initialisation.read_cssaln import read_morexaln_minimap
import os
import subprocess as sp
from .read_10x import read_10x_molecules
import numpy as np
import dask.dataframe as dd


def initial(args, popseq):

    if not os.path.exists(args.fasta + ".fai"):
        pysam.Fastafile(args.fasta)
    fai = pd.read_csv(args.fasta + ".fai", sep="\t",
                      usecols=[0, 1], names=["scaffold", "length"],
                      dtype={"scaffold": "str",
                             "length": np.int32}).reset_index(drop=False).rename(
        columns={"index": "scaffold_index"}).astype({"scaffold_index": np.int32})

    # fai.loc[:, "scaffold_index"] = pd.to_numeric(fai["scaffold_index"] + int(1), downcast="signed")
    # fai.loc[:, "length"] = pd.to_numeric(fai["length"], downcast="signed")
    fai["scaffold_index"] += 1
    fai["orig_scaffold_index"] = fai["scaffold_index"]
    fai["start"] = fai["orig_start"] = 1
    fai["end"] = fai["length"]
    fai["orig_end"] = fai["length"]
    fai = fai.set_index("scaffold_index")
    fai = dd.from_pandas(fai, npartitions=(fai.shape[0] // 1e5) + 1)

    # Alignment of genetic markers used for the genetic map. In this example, the Morex WGS assembly by IBSC (2012).
    cssaln = read_morexaln_minimap(paf=args.css, popseq=popseq, fai=fai, minqual=30, minlen=500, ref=True)
    # cssaln = cssaln.convert_dtypes()

    # Read the list of Hi-C links.
    fpairs_command = 'find {} -type f | grep "_fragment_pairs.tsv.gz$"'.format(args.hic)

    fpairs = [
            dd.read_csv(
                fname.decode().rstrip(), sep="\t", header=None, names=["scaffold1", "pos1", "scaffold2", "pos2"],
                compression="gzip", blocksize=None)
            for fname in sp.Popen(fpairs_command, shell=True, stdout=sp.PIPE).stdout
        ]

    fpairs = dd.concat(fpairs).reset_index(drop=True)
    # Now let us change the scaffold1 and scaffold2
    left = fai[["scaffold"]].reset_index(drop=False)
    left1 = left.rename(columns={"scaffold": "scaffold1", "scaffold_index": "scaffold_index1"})
    fpairs = dd.merge(left1, fpairs, how="right", on="scaffold1").drop("scaffold1", axis=1)
    left2 = left.copy()
    left2 = left.rename(columns={"scaffold": "scaffold2", "scaffold_index": "scaffold_index2"})
    fpairs = dd.merge(left2, fpairs, how="right", on="scaffold2").drop("scaffold2", axis=1)
    fpairs["orig_scaffold_index1"] = fpairs["scaffold_index1"]
    fpairs["orig_scaffold_index2"] = fpairs["scaffold_index2"]
    fpairs["orig_pos1"] = fpairs["pos1"]
    fpairs["orig_pos2"] = fpairs["pos2"]
    fpairs = fpairs.persist()

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
    cssaln = cssaln.persist()
    fpairs = fpairs.persist()
    molecules = molecules.persist()
    print(molecules.shape[0])
    fai = fai.persist()
    assert "scaffold" in fai.columns
    assembly = {"popseq": popseq,
                "fai": fai,
                "cssaln": cssaln,
                "10x_samples": samples,
                "fpairs": fpairs,
                "molecules": molecules,
                "barcodes": barcodes}

    return assembly
