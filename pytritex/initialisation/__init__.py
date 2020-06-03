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


def fai_reader(fasta, save_dir):
    fai = dd.read_csv(fasta + ".fai", sep="\t",
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
    fname = os.path.join(save_dir, "fai")
    dd.to_parquet(fai, os.path.join(save_dir, "fai"), compression="gzip", engine="pyarrow")
    fai = dd.read_parquet(fname, engine="pyarrow")
    return fai, fname


def read_fpairs(hic, fai, save_dir):
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
    fname = os.path.join(save_dir, "fpairs")
    dd.to_parquet(fpairs, fname, compression="gzip", compute=True, engine="pyarrow")
    return fname


def initial(popseq, fasta, css, tenx, hic, save, client, cores=1):

    save_dir = os.path.join(save, "joblib", "pytritex", "initialisation")
    os.makedirs(save_dir, exist_ok=True)
    popseq = dd.from_pandas(load(popseq), npartitions=10)
    popseq.columns = popseq.columns.str.replace("morex", "css")
    pop_name = os.path.join(save_dir, "popseq")
    dd.to_parquet(popseq, pop_name, compute=True)
    popseq = dd.read_parquet(pop_name)

    if not os.path.exists(fasta + ".fai"):
        pysam.Fastafile(fasta)

    print(time.ctime(), "Reading the FAIDX index")
    fai, fai_name = fai_reader(fasta, save_dir)
    print(time.ctime(), "Read the FAIDX index")
    # Alignment of genetic markers used for the genetic map. In this example, the Morex WGS assembly by IBSC (2012).
    print(time.ctime(), "Reading the CSS alignment")
    cssaln = read_morexaln_minimap(
        paf=css, popseq=popseq, save_dir=save_dir, fai=fai, minqual=30, minlen=500, ref=True)
    print(time.ctime(), "Read the CSS alignment")

    # Read the list of Hi-C links.
    print(time.ctime(), "Reading the HiC links")
    fpairs = read_fpairs(hic, fai, save_dir)
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
        molecules, barcodes = read_10x_molecules(samples, fai, save_dir, client=client)
        print(time.ctime(), "Read the 10X links")
    else:
        print("No 10X molecules! Error")
        print("Command:", tenx_command)
        import sys
        sys.exit(1)

    assembly = {"popseq": pop_name,
                "fai": fai_name,
                "cssaln": cssaln,
                "10x_samples": samples,
                "fpairs": fpairs,
                "molecules": molecules,
                "barcodes": barcodes}

    # Test it out
    print(dd.read_parquet(assembly["cssaln"]).head(5))
    return assembly
