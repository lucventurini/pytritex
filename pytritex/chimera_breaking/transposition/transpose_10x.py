import pandas as pd
from ...utils.rolling_join import rolling_join
from ...sequencing_coverage import add_molecule_cov
import dask.dataframe as dd
import os
from dask.distributed import Client
from joblib import Memory
import numpy as np
from ...utils import _rebalance_ddf


def _transpose_molecule_cov(new_assembly, fai: dd.DataFrame, save_dir: str,
                            assembly, client: Client, cores=1):

    molecules = assembly.get("molecule_cov", None)
    info = new_assembly["info"]
    assert isinstance(info, dd.DataFrame)
    if molecules is not None and molecules.shape[0].compute() > 0:
        scaffolds = fai[fai["derived_from_split"] == True].index.values.compute()
        old_to_keep = fai[fai["derived_from_split"] == False].index.values.compute()
        assert "mr_10x" not in info.columns
        coverage = add_molecule_cov(
            new_assembly, save_dir=save_dir, client=client, scaffolds=scaffolds,
            binsize=assembly["mol_binsize"], cores=cores)
        old_info = assembly["info"]
        assert old_info.index.name == "scaffold_index"
        _index = old_info.index.compute()
        present = np.unique(_index.values[_index.isin(old_to_keep)])
        old_info = old_info.loc[present]
        # We have to reset the index to trigger the sorting.
        assert old_info.index.name == coverage["info"].index.name == "scaffold_index"
        new_assembly["info"] = dd.concat([old_info.reset_index(drop=False),
                                          coverage["info"].reset_index(drop=False)
                                          ]).set_index("scaffold_index")

        new_assembly["mol_binsize"] = assembly["mol_binsize"]
        old_coverage = dd.read_parquet(assembly["molecule_cov"], infer_divisions=True)
        _index = old_coverage.index.compute()
        present = np.unique(_index.values[_index.isin(old_to_keep)])
        old_coverage = old_coverage.loc[present]
        if coverage["molecule_cov"].shape[0].compute() > 0:
            assert "scaffold_index" == old_coverage.index.name == coverage["molecule_cov"].index.name
            new_assembly["molecule_cov"] = dd.concat(
                [old_coverage.reset_index(drop=False),
                 coverage["molecule_cov"].reset_index(drop=False)
                 ]).set_index("scaffold_index")
        else:
            new_assembly["molecule_cov"] = old_coverage
        #
    else:
        new_assembly["molecule_cov"] = pd.DataFrame()
    return new_assembly


def _transpose_molecules(molecules: dd.DataFrame, fai: dd.DataFrame, save_dir: str) -> str:
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

    if isinstance(molecules, str):
        molecules = dd.read_parquet(molecules, infer_divisions=True)
    else:
        assert isinstance(molecules, dd.DataFrame)

    if isinstance(fai, str):
        fai = dd.read_parquet(fai, infer_divisions=True)

    if molecules is not None:
        orig_shape = molecules.shape[0].compute()
    else:
        orig_shape = 0

    if molecules is not None and orig_shape > 0:
        derived = fai.loc[fai.derived_from_split == True]
        to_keep = fai.loc[fai.derived_from_split == False]
        derived = derived[["orig_scaffold_index", "orig_start", "length"]].rename(
            columns={"length": "s_length"}).copy()
        derived = derived.reset_index(drop=False).assign(orig_pos=lambda df: df["orig_start"])
        molecules_up = molecules[molecules["orig_scaffold_index"].isin(
            to_keep["orig_scaffold_index"].compute())]
        assert molecules_up.index.name == "scaffold_index"
        molecules_down = molecules.loc[molecules["orig_scaffold_index"].isin(
            derived["orig_scaffold_index"].compute())].reset_index(drop=True)
        # assert molecules_up.shape[0].compute() + molecules_down.shape[0].compute() == orig_shape
        molecules_down = rolling_join(derived, molecules_down, on="orig_scaffold_index", by="orig_start")
        molecules_down["start"] = molecules_down.eval("orig_start - orig_pos + 1")
        molecules_down["end"] = molecules_down.eval("orig_end - orig_pos + 1")
        molecules_down = molecules_down.query("end <= s_length")
        # molecules_down = molecules_down.set_index("scaffold_index")
        assert "scaffold_index" in molecules_down
        # assert molecules_up.index.name == "scaffold_index", molecules_up.compute().head()
        molecules = dd.concat([molecules_up.reset_index(drop=False),
                               molecules_down])
        # assert molecules.index.values.compute().min() == min_idx
        molecules = molecules.drop("orig_pos", axis=1).drop("s_length", axis=1)
        # assert molecules.index.name == "scaffold_index", molecules.compute().head()
        # assert molecules.shape[0].compute() >= molecules_up.shape[0].compute()
        molecules = molecules.astype(dict((key, np.int) for key in
                                          ["start", "end", "length", "scaffold_index", "orig_scaffold_index",
                                           "orig_start", "orig_end"]))
        molecules = molecules.set_index("scaffold_index").persist()

    elif molecules is None:
        molecules = pd.DataFrame().assign(
            scaffold_index=[], barcode_index=[], start=[], end=[],
            npairs=[], sample=[], length=[], orig_start=[], orig_end=[], orig_scaffold_index=[])
        molecules = dd.from_pandas(molecules, npartitions=1)

    molecules = _rebalance_ddf(molecules, target_memory=5*10**7)
    fname = os.path.join(save_dir, "molecules")
    dd.to_parquet(molecules, fname, engine="pyarrow", compression="gzip", compute=True)
    return fname
