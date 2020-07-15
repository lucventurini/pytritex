import pandas as pd
from ...utils.rolling_join import rolling_join
from ...sequencing_coverage import add_molecule_cov
import dask.dataframe as dd
import os
from dask.distributed import Client
from joblib import Memory
import numpy as np
from ...utils import _rebalance_ddf


def _transpose_molecule_cov(new_assembly,
                            fai: dd.DataFrame,
                            save_dir: str,
                            scaffolds: np.ndarray,
                            assembly, client: Client, cores=1):

    molecules = assembly.get("molecule_cov", None)
    info = new_assembly["info"]
    assert isinstance(info, dd.DataFrame)
    if isinstance(molecules, str):
        molecules = dd.read_parquet(molecules, infer_divisions=True)
    if molecules is not None and molecules.shape[0].compute() > 0:

        assert "mr_10x" not in info.columns
        coverage = add_molecule_cov(
            new_assembly, save_dir=save_dir, client=client, scaffolds=scaffolds,
            binsize=assembly["mol_binsize"], cores=cores)
        old_info = assembly["info"]
        if isinstance(old_info, str):
            old_info = dd.read_parquet(old_info, infer_divisions=True)
        assert old_info.index.name == "scaffold_index"
        _index = np.unique(old_info.index.compute())
        present = _index[np.in1d(_index, fai[fai.to_use == True].index.values.compute()) &
                          ~np.in1d(_index, scaffolds)]
        old_info = old_info.loc[present]
        # We have to reset the index to trigger the sorting.
        assert old_info.index.name == coverage["info"].index.name == "scaffold_index"
        new_assembly["info"] = dd.concat([old_info.reset_index(drop=False),
                                          coverage["info"].reset_index(drop=False)
                                          ]).set_index("scaffold_index")

        new_assembly["mol_binsize"] = assembly["mol_binsize"]
        old_coverage = dd.read_parquet(assembly["molecule_cov"], infer_divisions=True)
        _index = np.unique(old_coverage.index.compute())
        present = _index[np.in1d(_index, fai[fai.to_use == True].index.values.compute()) &
                         ~np.in1d(_index, scaffolds)]
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

        molecules_index = molecules.index.unique().compute()
        mol_up_index = np.unique(
            molecules_index.intersection(fai.loc[fai["to_use"] == True].index.compute().values).values)
        mol_down_index = molecules_index.difference(mol_up_index).values
        molecules_up = molecules.loc[mol_up_index]
        assert molecules_up.index.name == "scaffold_index"
        molecules_down = molecules.loc[mol_down_index].reset_index(drop=True).set_index("orig_scaffold_index")
        derived = fai.loc[
            (fai.orig_scaffold_index.isin(fai.loc[mol_down_index, "orig_scaffold_index"].values.compute())) & (
                        fai["to_use"] == True), ["orig_scaffold_index", "orig_start", "length"]][:]
        derived = derived.rename(columns={"length": "s_length"}).reset_index(drop=False)
        derived = derived.assign(orig_pos=derived["orig_start"])
        molecules_down = rolling_join(derived, molecules_down, on="orig_scaffold_index", by="orig_start")

        molecules_down["start"] = molecules_down.eval("orig_start - orig_pos + 1")
        molecules_down["end"] = molecules_down.eval("orig_end - orig_pos + 1")
        # Remove 10X links that straddle the broken chimerism.
        molecules_down = molecules_down.query("end <= s_length").set_index("scaffold_index")
        # molecules_down = molecules_down.set_index("scaffold_index")
        molecules_down = molecules_down.drop("orig_pos", axis=1).drop("s_length", axis=1)
        # assert molecules_up.index.name == "scaffold_index", molecules_up.compute().head()
        molecules = dd.concat([molecules_up,
                               molecules_down]).reset_index(drop=False)
        # assert molecules.index.values.compute().min() == min_idx
        # assert molecules.index.name == "scaffold_index", molecules.compute().head()
        # assert molecules.shape[0].compute() >= molecules_up.shape[0].compute()
        molecules = molecules.astype(dict((key, np.int) for key in
                                          ["start", "end", "length", "scaffold_index", "orig_scaffold_index",
                                           "orig_start", "orig_end"]))
        molecules = molecules.set_index("scaffold_index")
        assert isinstance(molecules, dd.DataFrame)

    elif molecules is None:
        molecules = pd.DataFrame().assign(
            scaffold_index=[], barcode_index=[], start=[], end=[],
            npairs=[], sample=[], length=[], orig_start=[], orig_end=[], orig_scaffold_index=[])
        molecules = dd.from_pandas(molecules, npartitions=1)

    assert isinstance(molecules, dd.DataFrame)
    molecules = _rebalance_ddf(molecules, target_memory=5*10**7)
    fname = os.path.join(save_dir, "molecules")
    dd.to_parquet(molecules, fname, engine="pyarrow", compression="gzip", compute=True)
    return fname
