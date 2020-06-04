import pandas as pd
import dask.dataframe as dd
import numpy as np
import os


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


def read_10x_molecules(samples: pd.DataFrame, fai: pd.DataFrame, save_dir, client):
    """Read the files as produced by run_10x_mapping.zsh"""
    # pool = mp.Pool(ncores)

    molecules = [client.submit(
        _10xreader, row) for row in samples[["index", "fname"]].itertuples(index=False, name=None)]
    molecules = client.gather(molecules)
    mol = dd.concat(molecules)
    f = fai[["scaffold"]].reset_index(drop=False).set_index("scaffold")
    try:
        mol = dd.merge(f, mol, on="scaffold", how="right").reset_index(drop=True)
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
    mol["orig_scaffold_index"] = mol.index.values
    barcode = mol["barcode"].compute()
    barcodes = pd.DataFrame({"barcode_index": np.arange(barcode.shape[0], dtype=np.int32),
                             "barcode": barcode})
    mol = dd.merge(barcodes, mol, how="right", on="barcode").drop("barcode", axis=1).set_index("scaffold_index")
    fname = os.path.join(save_dir, "molecules")
    dd.to_parquet(mol, fname, compression="gzip", engine="pyarrow", compute=True)
    return fname, barcodes
