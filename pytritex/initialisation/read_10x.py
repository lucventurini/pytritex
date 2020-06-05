import pandas as pd
import dask.dataframe as dd
from dask.distributed import Client
import numpy as np
import os
import time
from ..utils import parse_size, return_size
from dask import delayed


def _10xreader(item):
    sample, fname = item
    df = dd.read_csv(fname, header=None, names=("scaffold", "start", "end", "barcode", "npairs"), sep="\t",
                     compression="gzip", blocksize=None,
                     dtype={"scaffold": str,
                            "barcode": str,
                            "start": np.uint32, "end": np.uint32, "npairs": np.uint32},
                     ).set_index("scaffold")
    df["sample"] = sample
    # df = df.reset_index(drop=False).set_index("index")
    # for key in ["start", "end", "npairs"]:
    #     df[key] = pd.to_numeric(df[key].compute(), downcast="signed")
    return df


def read_10x_molecules(samples: pd.DataFrame, fai: pd.DataFrame, save_dir, client: Client,
                       cores=1, memory="20GB"):
    """Read the files as produced by run_10x_mapping.zsh"""
    # pool = mp.Pool(ncores)

    proc_memory, unit = parse_size(memory)
    proc_memory = return_size(proc_memory / cores, unit)

    molecules = [client.submit(
        _10xreader, row) for row in samples[["index", "fname"]].itertuples(index=False, name=None)]
    molecules = client.gather(molecules)
    print(time.ctime(), "Concatenating 10X files")
    mol = dd.concat(molecules)
    print(time.ctime(), "Concatenated 10X files")
    mol = client.scatter(mol)
    print(time.ctime(), "Scattered 10X files")
    f = fai[["scaffold"]].reset_index(drop=False).set_index("scaffold")
    try:
        print(time.ctime(), "Merging FAI into molecules")
        molfunc = delayed(dd.merge)(f, mol, on="scaffold", how="right")
        mol = client.compute(molfunc).result()
        mol = mol.reset_index(drop=True)
        print(time.ctime(), "Merged FAI into molecules")
    except KeyError:
        print(mol.head(10))
        print()
        print(f.head(10))
        print()
        raise
    mol["length"] = mol.eval("end - start")
    mol["start"] += 1
    mol["orig_start"] = mol["start"]
    mol["orig_end"] = mol["end"]
    mol["orig_scaffold_index"] = mol["scaffold_index"]
    barcode = mol["barcode"].compute()
    print(time.ctime(), "Changing barcodes with indices")
    barcodes = pd.DataFrame({"barcode_index": np.arange(barcode.shape[0], dtype=np.int32),
                             "barcode": barcode})
    molfunc = delayed(dd.merge)(barcodes, mol, how="right", on="barcode",
                                npartitions=100)
    mol = client.compute(molfunc)
    mol = dd.from_delayed([mol])
    mol = mol.drop("barcode", axis=1).set_index("scaffold_index")
    print(time.ctime(), "Storing data to disk")
    fname = os.path.join(save_dir, "molecules")
    dd.to_parquet(mol, fname, compression="gzip", engine="pyarrow", compute=True)
    return fname, barcodes
