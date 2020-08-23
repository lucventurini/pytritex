import pysam
import pandas as pd
from pytritex.initialisation.read_cssaln import read_morexaln_minimap
import os
import subprocess as sp
from .read_10x import read_10x_molecules
import dask.dataframe as dd
from dask import delayed
import time
from joblib import load
import numpy as np
import dask.array as da
import logging
dask_logger = logging.getLogger("dask")


def fai_reader(fasta, save_dir):
    fasta = pysam.Fastafile(fasta)
    fai = pd.DataFrame().assign(
        scaffold=fasta.references,
        length=fasta.lengths).astype({"scaffold": "str", "length": np.int32})

    fai["scaffold_index"] = pd.Series(range(1, fai["length"].shape[0] + 1))
    fai["orig_scaffold_index"] = fai["scaffold_index"]
    fai["start"] = fai["orig_start"] = 1
    fai["end"] = fai["length"]
    fai["orig_end"] = fai["length"]
    fai = fai.set_index("scaffold_index")
    fname = os.path.join(save_dir, "fai")
    dd.to_parquet(dd.from_pandas(fai,
                                 chunksize=int(1e5)), os.path.join(save_dir, "fai"),
                  compression="gzip", engine="pyarrow")
    fai = dd.read_parquet(fname, engine="pyarrow")
    return fai, fname


def column_switcher(fpairs: dd.DataFrame):
    query = "scaffold1 == scaffold2 & pos1 > pos2"
    dask_logger.warning("%s Extracting the positions for column switching", time.ctime())
    values = fpairs[["pos1", "pos2"]].to_dask_array(lengths=True)
    dask_logger.warning("%s Extracted the positions for column switching", time.ctime())    
    if values.shape[0] == 0:
        return fpairs
    dask_logger.warning("%s Creating the mask", time.ctime())
    mask = np.repeat(fpairs.eval(query).compute().to_numpy(), 2).reshape(values.shape)
    dask_logger.warning("%s Created the mask, creating the values", time.ctime())
    values = np.where(mask, values[:, [1, 0]], values)
    dask_logger.warning("%s Created the values, switching columns", time.ctime())    
    fpairs["pos1"] = values[:, 0]
    fpairs["pos2"] = values[:, 1]
    dask_logger.warning("%s Switched columns, dropping duplicates", time.ctime())
    fpairs = fpairs.drop_duplicates()
    dask_logger.warning("%s Dropped duplicates, finished", time.ctime())
    return fpairs


def read_fpairs(hic, fai, save_dir):
    fpairs_command = 'find {} -type f | grep "_fragment_pairs.tsv.gz$"'.format(hic)
    fpairs = [fname.decode().rstrip() for fname in sp.Popen(fpairs_command, shell=True, stdout=sp.PIPE).stdout]
    fnames = []
    
    if len(fpairs) > 0:
        dask_logger.warning("%s Decompressing the FPAIRS file(s)")
        for num, fname in enumerate(fpairs, 1):
            buf = open(os.path.join(save_dir, "{}.tsv".format(num)), "wb")
            fnames.append(buf.name)
            command = "zcat {fname}".format(**locals())
            reader = sp.Popen(command, shell=True, stdout=buf)
            reader.communicate()
            buf.close()
        
        dask_logger.warning("%s Starting to read the CSV(s): %s.", time.ctime(), ",".join(fnames))
        fpairs = dd.read_csv(fnames, sep="\t", header=None, names=["scaffold1", "pos1", "scaffold2", "pos2"])
        fpairs["pair_index"] = da.from_array(np.arange(1, fpairs.shape[0].compute() + 1, dtype=np.int),
                                             chunks=tuple(fpairs.map_partitions(len).compute().values.tolist()))
        dask_logger.warning("%s Indexing.", time.ctime())
        fpairs = fpairs.set_index("pair_index", sorted=True)
        dask_logger.warning("%s Switching columns.", time.ctime())
        fpairs = column_switcher(fpairs)
        dask_logger.warning("%s Merging on scaffold1.", time.ctime())
        # Now let us change the scaffold1 and scaffold2
        left = fai[["scaffold"]].reset_index(drop=False)
        left1 = left[:].rename(columns={"scaffold_index": "scaffold_index1", "scaffold": "scaffold1"})
        # assert "scaffold1" in left1.columns and "scaffold_index1" in left1.columns and left1.columns.shape[0] == 2
        fpairs = dd.merge(left1, fpairs, how="right", on="scaffold1").drop("scaffold1", axis=1).persist()
        dask_logger.warning("%s Merging on scaffold2.", time.ctime())        
        left2 = left[:].rename(columns={"scaffold_index": "scaffold_index2", "scaffold": "scaffold2"})
        # assert "scaffold2" in left2.columns and "scaffold_index2" in left2.columns and left2.columns.shape[0] == 2
        fpairs = dd.merge(left2, fpairs, how="right", on="scaffold2").drop("scaffold2", axis=1).persist()
        dask_logger.warning("%s Merged on scaffold2.", time.ctime())
        fpairs["orig_scaffold_index1"] = fpairs["scaffold_index1"]
        fpairs["orig_scaffold_index2"] = fpairs["scaffold_index2"]
        fpairs["orig_pos1"] = fpairs["pos1"]
        fpairs["orig_pos2"] = fpairs["pos2"]
        # Now create an index.
        dask_logger.warning("%s Creating the pair_index column.", time.ctime())
        fpairs = fpairs.drop("pair_index", axis=1, errors="ignore")
        fpairs["pair_index"] = da.from_array(np.arange(1, fpairs.shape[0].compute() + 1, dtype=np.int),
                                             chunks=tuple(fpairs.map_partitions(len).compute().values.tolist()))
        dask_logger.warning("%s Setting the index on pair_index.", time.ctime())        
        fpairs = fpairs.set_index("pair_index")
        dask_logger.warning("%s Created and indexed with the pair_index column.", time.ctime())        
        # Remove double lines.
    else:
        fpairs = pd.DataFrame().assign(pair_index=[], scaffold_index1=[], scaffold_index2=[],
                                       pos1=[], pos2=[], orig_pos1=[], orig_pos2=[],
                                       orig_scaffold_index1=[], orig_scaffold_index2=[]).set_index("pair_index")
        fpairs = dd.from_pandas(fpairs, chunksize=10 ** 5)
    fname = os.path.join(save_dir, "fpairs")
    dd.to_parquet(fpairs, fname, compression="gzip", compute=True, engine="pyarrow")
    [os.remove(fname) for fname in fnames]    
    return fname


def load_popseq(popseq, save_dir):
    popseq = dd.from_pandas(load(popseq), npartitions=10)
    popseq.columns = popseq.columns.str.replace("morex", "css")
    pop_name = os.path.join(save_dir, "popseq")
    dd.to_parquet(popseq, pop_name, compute=True)
    popseq = dd.read_parquet(pop_name)
    return popseq, pop_name


def initial(popseq, fasta, css, tenx, hic, save, client, memory, ram="20GB", ref=None,
            cores=1):

    save_dir = os.path.join(save, "joblib", "pytritex", "initialisation")
    os.makedirs(save_dir, exist_ok=True)
    def submit_popseq(popseq, save_dir):
        pop_function = delayed(load_popseq)(popseq, save_dir)
        popped = client.compute(pop_function)
        # result = client.compute(popped)
        popseq, pop_name = popped.result()
        return popseq, pop_name

    popseq, pop_name = memory.cache(submit_popseq)(popseq, save_dir)

    def fai_submitter(fasta, save_dir):
        func = delayed(fai_reader)(fasta, save_dir)
        fairead = client.compute(func)
        fai, fai_name = fairead.result()
        return fai, fai_name

    print(time.ctime(), "Reading the FAIDX index")
    fai, fai_name = memory.cache(fai_submitter)(fasta, save_dir)
    print(time.ctime(), "Read the FAIDX index")
    # Alignment of genetic markers used for the genetic map. In this example, the Morex WGS assembly by IBSC (2012).
    print(time.ctime(), "Reading the CSS alignment")

    def cssaln_submitter(fai, css, popseq, save_dir, minqual, minlen, ref):
        cssaln = client.submit(read_morexaln_minimap,
                               paf=css, popseq=popseq, save_dir=save_dir,
                               fai=fai, minqual=minqual, minlen=minlen, ref=ref)
        cssaln = client.compute(cssaln).result()
        return cssaln

    cssaln = memory.cache(cssaln_submitter, ignore=["fai", "css", "popseq"])(
        fai, css, popseq, save_dir, minqual=30, minlen=500, ref=ref)
    print(time.ctime(), "Read the CSS alignment")
    # Read the list of Hi-C links.
    print(time.ctime(), "Reading the HiC links")

    def fpairs_reader_submitter(hic, fai, save_dir):
        fpairs = client.submit(read_fpairs, hic, fai, save_dir)
        fpairs = client.compute(fpairs).result()
        return fpairs

    fpairs = memory.cache(fpairs_reader_submitter, ignore=["hic", "fai"])(hic, fai, save_dir)
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
        molecules, barcodes = memory.cache(read_10x_molecules,
                                           ignore=["fai", "client",
                                                   "cores", "memory"])(samples, fai, save_dir, client=client,
                                                                       memory=ram, cores=cores)
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
