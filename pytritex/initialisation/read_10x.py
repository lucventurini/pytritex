import pandas as pd
import numpy as np
import multiprocessing as mp
from time import ctime


def _10xreader(item):
    sample, fname = item
    df = pd.read_csv(fname, header=None, names=("scaffold", "start", "end", "barcode", "npairs"), sep="\t")
    df.loc[:, "sample"] = sample
    for key in ["start", "end", "npairs"]:
        df.loc[:, key] = pd.to_numeric(df[key], downcast="unsigned")
    return df


def read_10x_molecules(samples: pd.DataFrame, fai: pd.DataFrame, ncores=1):
    """Read the files as produced by run_10x_mapping.zsh"""
    pool = mp.Pool(ncores)
    mol = pd.concat(pool.map(_10xreader, samples[["index", "fname"]].itertuples(index=False, name=None)))
    pool.close()
    mol = fai[["scaffold", "scaffold_index"]].merge(mol, on=["scaffold"], how="right")
    assert (~mol["scaffold_index"].isna()).all()
    mol.drop("scaffold", axis=1, inplace=True)
    mol["length"] = mol.end - mol.start
    mol.start += 1
    mol.loc[:, "orig_start"] = mol["start"]
    mol.loc[:, "orig_end"] = mol["end"]
    mol.loc[:, "orig_scaffold_index"] = mol["scaffold_index"]
    barcodes = pd.DataFrame({"barcode_index": np.arange(mol.shape[0], dtype=np.int),
                             "barcode": mol["barcode"]})
    barcodes.loc[:, "barcode_index"] = pd.to_numeric(barcodes["barcode_index"], downcast="unsigned")
    mol = barcodes.merge(mol, how="right", on="barcode").drop("barcode", axis=1)
    return mol, barcodes
