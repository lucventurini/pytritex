import pandas as pd
import numpy as np
import multiprocessing as mp


def _10xreader(item):
    sample, fname = item
    df = pd.read_csv(fname, header=None, names=("scaffold", "start", "end", "barcode", "npairs"))
    df["sample"] = sample
    return df


def read_10x_molecules(files: dict, ncores=1):
    """Read the files as produced by run_10x_mapping.zsh"""
    pool = mp.Pool(ncores)
    mol = pd.concat(pool.map(_10xreader, files.items()))
    mol["length"] = mol.end - mol.start
    mol.start += 1
    return mol
