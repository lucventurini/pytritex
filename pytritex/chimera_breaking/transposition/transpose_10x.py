import pandas as pd
from ...utils.rolling_join import rolling_join
from ...sequencing_coverage import add_molecule_cov
import dask.dataframe as dd
import os
from dask.distributed import Client
# from joblib import Memory
import numpy as np
import logging
dask_logger = logging.getLogger("dask")
import time


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
        dask_logger.debug("%s Starting transposing the molecules; npartitions: %s",
                            time.ctime(), molecules.npartitions)
        molecules_index = molecules.index.unique().compute()
        mol_up_index = np.unique(
            molecules_index.intersection(fai.loc[fai["to_use"] == True].index.compute().values).values)
        mol_down_index = molecules_index.difference(mol_up_index).values
        molecules_up = molecules.loc[mol_up_index].persist()
        assert molecules_up.index.name == "scaffold_index"
        molecules_down = molecules.loc[mol_down_index].reset_index(drop=True).set_index("orig_scaffold_index").persist()
        derived = fai.loc[
            (fai.orig_scaffold_index.isin(fai.loc[mol_down_index, "orig_scaffold_index"].values.compute())) & (
                        fai["to_use"] == True), ["orig_scaffold_index", "orig_start", "length"]][:]
        derived = derived.rename(columns={"length": "s_length"}).reset_index(drop=False)
        derived = derived.assign(orig_pos=derived["orig_start"])
        molecules_down = rolling_join(derived, molecules_down,
                                      on="orig_scaffold_index", by="orig_start")
        assert "orig_scaffold_index" in molecules_down.columns
        assert "scaffold_index" in molecules_down.columns
        molecules_down = molecules_down.set_index("scaffold_index")
        molecules_down["start"] = molecules_down.eval("orig_start - orig_pos + 1")
        molecules_down["end"] = molecules_down.eval("orig_end - orig_pos + 1")
        # Remove 10X links that straddle the broken chimerism.
        molecules_down = molecules_down.query("end <= s_length")
        # molecules_down = molecules_down.set_index("scaffold_index")
        molecules_down = molecules_down.drop("orig_pos", axis=1).drop("s_length", axis=1)
        # assert molecules_up.index.name == "scaffold_index", molecules_up.compute().head()
        nparts = molecules.npartitions
        orig_columns = molecules.columns[:]
        molecules = dd.concat([molecules_up,
                               molecules_down]).reset_index(drop=False)
        assert molecules.scaffold_index.isna().any().compute() == False
        molecules = molecules.astype({"scaffold_index": np.int})
        molecules = molecules.set_index("scaffold_index")
        molecules = molecules.astype(dict((key, np.int) for key in
                                          ["start", "end", "length", "orig_scaffold_index", "orig_start", "orig_end"]))
        molecules = molecules.repartition(npartitions=nparts)
        molecules = molecules[orig_columns]
        assert isinstance(molecules, dd.DataFrame)
        dask_logger.debug("%s Finished transposing the molecules; npartitions: %s", time.ctime(), molecules.npartitions)

    elif molecules is None:
        molecules = pd.DataFrame().assign(
            scaffold_index=[], barcode_index=[], start=[], end=[],
            npairs=[], sample=[], length=[], orig_start=[], orig_end=[], orig_scaffold_index=[])
        molecules = dd.from_pandas(molecules, npartitions=1)

    assert isinstance(molecules, dd.DataFrame)
    fname = os.path.join(save_dir, "molecules")
    dd.to_parquet(molecules, fname, engine="pyarrow", compression="gzip", compute=True)
    return fname
