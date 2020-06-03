import pysam
import pandas as pd
from pytritex.initialisation.read_cssaln import read_morexaln_minimap
import os
import subprocess as sp
from .read_10x import read_10x_molecules
import numpy as np
import dask.dataframe as dd
import time
from joblib import load


def fai_reader(fasta):
    fai = pd.read_csv(fasta + ".fai", sep="\t",
                      usecols=[0, 1], names=["scaffold", "length"],
                      dtype={"scaffold": "str",
                             "length": np.int32}).reset_index(drop=False).rename(
        columns={"index": "scaffold_index"}).astype({"scaffold_index": np.int32})
    fai["scaffold_index"] += 1
    fai["orig_scaffold_index"] = fai["scaffold_index"]
    fai["start"] = fai["orig_start"] = 1
    fai["end"] = fai["length"]
    fai["orig_end"] = fai["length"]
    fai = fai.set_index("scaffold_index")
    fai = dd.from_pandas(fai, npartitions=(fai.shape[0] // 1e5) + 1)
    return fai


def read_fpairs(hic, fai):
    fpairs_command = 'find {} -type f | grep "_fragment_pairs.tsv.gz$"'.format(hic)
    fpairs = [
        dd.read_csv(
            fname.decode().rstrip(), sep="\t", header=None, names=["scaffold1", "pos1", "scaffold2", "pos2"],
            compression="gzip", blocksize=None)
        for fname in sp.Popen(fpairs_command, shell=True, stdout=sp.PIPE).stdout
    ]
    if len(fpairs) > 0:
        fpairs = dd.concat(fpairs).reset_index(drop=True)
        # Now let us change the scaffold1 and scaffold2
        left = fai[["scaffold"]].reset_index(drop=False)
        left1 = left.rename(columns={"scaffold": "scaffold1", "scaffold_index": "scaffold_index1"})
        fpairs = dd.merge(left1, fpairs, how="right", on="scaffold1").drop("scaffold1", axis=1)
        left2 = left.copy()
        left2 = left2.rename(columns={"scaffold": "scaffold2", "scaffold_index": "scaffold_index2"})
        fpairs = dd.merge(left2, fpairs, how="right", on="scaffold2").drop("scaffold2", axis=1)
        fpairs["orig_scaffold_index1"] = fpairs["scaffold_index1"]
        fpairs["orig_scaffold_index2"] = fpairs["scaffold_index2"]
        fpairs["orig_pos1"] = fpairs["pos1"]
        fpairs["orig_pos2"] = fpairs["pos2"]
    else:
        fpairs = pd.DataFrame().assign(scaffold_index1=[], scaffold_index2=[],
                                       pos1=[], pos2=[], orig_pos1=[], orig_pos2=[],
                                       orig_scaffold_index1=[], orig_scaffold_index2=[])
        fpairs = dd.from_pandas(fpairs)
    fpairs = fpairs.persist()
    return fpairs


def initial(popseq, fasta, css, tenx, hic, save, cores=1):

    popseq = dd.from_pandas(load(popseq), npartitions=10)
    popseq.columns = popseq.columns.str.replace("morex", "css")

    if not os.path.exists(fasta + ".fai"):
        pysam.Fastafile(fasta)

    print(time.ctime(), "Reading the FAIDX index")
    fai = fai_reader(fasta)
    print(time.ctime(), "Read the FAIDX index")
    # Alignment of genetic markers used for the genetic map. In this example, the Morex WGS assembly by IBSC (2012).
    print(time.ctime(), "Reading the CSS alignment")
    cssaln = read_morexaln_minimap(
        paf=css, popseq=popseq, fai=fai, minqual=30, minlen=500, ref=True)
    print(time.ctime(), "Read the CSS alignment")
    # cssaln = cssaln.convert_dtypes()

    # Read the list of Hi-C links.
    print(time.ctime(), "Reading the HiC links")
    fpairs = read_fpairs(hic, fai)
    print(time.ctime(), "Read the HiC links")

    tenx_command = 'find {} -type f | grep "molecules.tsv.gz$"'.format(tenx)
    tenx_files = [line.rstrip().decode() for line in
                  sp.Popen(tenx_command, shell=True, stdout=sp.PIPE).stdout]
    if tenx_files:
        print(time.ctime(), "Reading the 10X links")
        tenx_dict = []
        for index, fname in enumerate(tenx_files):
            sample = os.path.basename(os.path.dirname(fname))
            tenx_dict.append((index, sample, fname))
        samples = list(zip(*tenx_dict))
        samples = pd.DataFrame({"index": samples[0], "sample": samples[1], "fname": samples[2]})
        molecules, barcodes = read_10x_molecules(samples, fai, ncores=cores)
        print(time.ctime(), "Read the 10X links")
    else:
        print("No 10X molecules! Error")
        print("Command:", tenx_command)
        import sys
        sys.exit(1)

    save_dir = os.path.join(save, "joblib", "pytritex", "initialisation")
    os.makedirs(save_dir, exist_ok=True)

    dd.to_parquet(cssaln, os.path.join(save_dir, "cssaln"), compression="gzip")
    cssaln = os.path.join(save_dir, "cssaln")
    dd.to_parquet(fpairs, os.path.join(save_dir, "fpairs"), compression="gzip")
    fpairs = os.path.join(save_dir, "fpairs")
    dd.to_parquet(molecules, os.path.join(save_dir, "molecules"), compression="gzip")
    molecules = os.path.join(save_dir, "molecules")
    dd.to_parquet(fai, os.path.join(save_dir, "fai"), compression="gzip")
    barcodes = dd.to_parquet(dd.from_pandas(barcodes, npartitions=100),
                             os.path.join(save_dir, "barcodes"), compression="gzip")
    popseq = dd.to_parquet(popseq, os.path.join(save_dir, "popseq"), compression="gzip")

    fai = os.path.join(save_dir, "fai")
    assembly = {"popseq": popseq,
                "fai": fai,
                "cssaln": cssaln,
                "10x_samples": samples,
                "fpairs": fpairs,
                "molecules": molecules,
                "barcodes": barcodes}

    return assembly
