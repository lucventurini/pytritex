import pandas as pd
import numpy as np
import multiprocessing as mp
from time import ctime


def _10xreader(item):
    sample, fname = item
    df = pd.read_csv(fname, header=None, names=("scaffold", "start", "end", "barcode", "npairs"), sep="\t")
    df["sample"] = sample
    return df


def read_10x_molecules(files: list, fai: pd.DataFrame, ncores=1):
    """Read the files as produced by run_10x_mapping.zsh"""
    pool = mp.Pool(ncores)
    mol = pd.concat(pool.map(_10xreader, files))
    pool.close()
    mol = fai.loc[["scaffold", "scaffold_index"]].merge(mol, on=["scaffold"], how="right")
    assert ~mol["scaffold_index"].isna()
    mol.drop("scaffold", axis=1, inplace=True)
    mol["length"] = mol.end - mol.start
    mol.start += 1
    return mol
