import pandas as pd
from ...utils.rolling_join import rolling_join
from ...sequencing_coverage import add_molecule_cov
import dask.dataframe as dd
import os
from dask.distributed import Client
from joblib import Memory
import numpy as np


def _transpose_molecule_cov(new_assembly, fai, assembly, save_dir: str, memory: Memory,
                            client: Client, cores=1):

    molecules = assembly.get("molecule_cov", None)
    if molecules is not None:
        molecules = dd.read_parquet(molecules)
    fai = dd.read_parquet(fai)
    info = dd.read_parquet(new_assembly["info"])

    if molecules is not None and molecules.shape[0].compute() > 0:
        scaffolds = fai[fai["derived_from_split"] == True].index.values.compute()
        old_to_keep = fai[fai["derived_from_split"] == False].index.values.compute()
        assert "mr_10x" not in info.columns
        coverage = memory.cache(add_molecule_cov,
                                ignore=["cores", "client"])(
            new_assembly, save_dir=save_dir, client=client, scaffolds=scaffolds,
            binsize=assembly["mol_binsize"], cores=cores)
        old_info = assembly["info"].loc[assembly["info"].scaffold_index.isin(old_to_keep)]
        new_assembly["info"] = pd.concat([old_info, coverage["info"]]).reset_index(drop=True)
        new_assembly["mol_binsize"] = assembly["mol_binsize"]
        old_coverage = assembly["molecule_cov"].loc[assembly["molecule_cov"].scaffold_index.isin(old_to_keep)]
        if coverage["molecule_cov"].shape[0] > 0:
            new_assembly["molecule_cov"] = pd.concat([old_coverage, coverage["molecule_cov"]]).reset_index(drop=True)
        else:
            new_assembly["molecule_cov"] = old_coverage
        for key in ["cov", "info"]:
            fname = os.path.join(save_dir, key + "_hic")
            # BugFix for pyarrow not handling float16 and int16
            if (new_assembly[key].dtypes == np.int16).any():
                for index, col in enumerate(new_assembly[key].dtypes == np.int16):
                    if col is False:
                        continue
                    col = new_assembly[key].dtypes.index[index]
                    new_assembly[key][col] = new_assembly[key][col].astype(np.int32)
            if (new_assembly[key].dtypes == np.float16).any():
                for index, col in enumerate(new_assembly[key].dtypes == np.float16):
                    if col is False:
                        continue
                    col = new_assembly[key].dtypes.index[index]
                    new_assembly[key][col] = new_assembly[key][col].astype(np.float32)

            dd.to_parquet(new_assembly[key], fname, compression="gzip", engine="pyarrow", compute=True)
            new_assembly[key] = fname

    else:
        new_assembly["molecule_cov"] = pd.DataFrame()
    return new_assembly


def _transpose_molecules(molecules: str, fai: str, save_dir: str):
    # if("molecules" %in% names(assembly) && nrow(molecules) > 0){
    #   cat("Transpose molecules\n")
    #   copy(molecules) -> z
    #   z[, scaffold := NULL]
    #   fai[, .(scaffold, orig_scaffold, orig_start, s_length=orig_end - orig_start + 1, orig_pos=orig_start)][z, on=c("orig_scaffold", "orig_start"), roll=T]->z
    #   z[, start := orig_start - orig_pos + 1]
    #   z[, end := orig_end - orig_pos + 1]
    #   z[end <= s_length]->z
    #   z[, orig_pos := NULL]
    #   z[, s_length  := NULL]
    #   assembly_new$molecules <- z
    #  } else {
    #   assembly_new$molecules <- data.table()
    #  }

    fai = dd.read_parquet(fai)
    if molecules is not None:
        molecules = dd.read_parquet(molecules)

    if molecules is not None and molecules.shape[0].compute() > 0:
        derived = fai.loc[fai.derived_from_split == True]
        to_keep = fai.loc[fai.derived_from_split == False]
        derived = derived[["orig_scaffold_index", "orig_start", "length"]].rename(
            columns={"length": "s_length"}).copy()
        derived = derived.reset_index(drop=False).assign(orig_pos=lambda df: df["orig_start"])
        molecules_up = molecules[molecules["orig_scaffold_index"].isin(
            to_keep["orig_scaffold_index"].compute())]
        assert molecules.index.name == "scaffold_index"
        molecules_down = molecules.loc[molecules["orig_scaffold_index"].isin(
            derived["orig_scaffold_index"].compute())].reset_index(drop=True)
        assert molecules_up.shape[0].compute() + molecules_down.shape[0].compute() == molecules.shape[0].compute()
        molecules_down = rolling_join(derived, molecules_down, on="orig_scaffold_index", by="orig_start")
        try:
            molecules_down["start"] = molecules_down.eval("orig_start - orig_pos + 1")
            molecules_down["end"] = molecules_down.eval("orig_end - orig_pos + 1")
        except pd.core.computation.ops.UndefinedVariableError:
            print(molecules_up.compute().head())
            print("\n###\n")
            print(molecules_down.head())
            print("\n###\n")
            print(derived.compute().head())
            print("\n###")
            raise

        molecules_down = molecules_down.set_index("scaffold_index")
        molecules = dd.concat([molecules_up, molecules_down]).reset_index(drop=True)
        molecules = molecules.loc[molecules["end"] <= molecules["s_length"]].drop("orig_pos", axis=1).drop(
            "s_length", axis=1)
    elif molecules is None:
        molecules = pd.DataFrame().assign(
            scaffold_index=[], barcode_index=[], start=[], end=[],
            npairs=[], sample=[], length=[], orig_start=[], orig_end=[], orig_scaffold_index=[])
        molecules = dd.from_pandas(molecules, npartitions=1)

    fname = os.path.join(save_dir, "molecules")
    dd.to_parquet(molecules, fname, engine="pyarrow", compression="gzip", compute=True)
    return fname
