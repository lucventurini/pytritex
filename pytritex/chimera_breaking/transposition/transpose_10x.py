import pandas as pd
from ...utils.rolling_join import rolling_join
from ...sequencing_coverage import add_molecule_cov
import dask.dataframe as dd
import os


def _transpose_molecule_cov(new_assembly, fai, assembly, cores=1):
    #  if("molecule_cov" %in% names(assembly) & nrow(molecules) > 0){
    #   cat("10X molecule coverage\n")
    #   add_molecule_cov(assembly_new, scaffolds=fai[split == T]$scaffold, binsize=assembly$mol_binsize, cores=cores)->cov
    #   info[!breaks$scaffold, on="scaffold"]->x
    #   x[, scaffold := paste0(prefix, sub(regex2, "\\2", scaffold))]
    #   x[, split := F]
    #   rbind(x[, names(cov$info), with=F], cov$info)->assembly_new$info
    #   assembly_new$mol_binsize <- assembly$mol_binsize
    #
    #   assembly$molecule_cov[!breaks$scaffold, on="scaffold"]->x
    #   x[, scaffold := paste0(prefix, sub(regex2, "\\2", scaffold))]
    #   if(nrow(cov$molecule_cov) > 0){
    #    rbind(x, cov$molecule_cov)->assembly_new$molecule_cov
    #   } else {
    #    x -> assembly_new$molecule_cov
    #   }
    #  } else {
    #   assembly_new$molecule_cov <- data.table()
    #  }

    if assembly.get("molecule_cov", None) is not None and assembly["molecule_cov"].shape[0] > 0:
        scaffolds = fai.loc[fai["derived_from_split"], "scaffold_index"]
        old_to_keep = fai.loc[~fai["derived_from_split"], "scaffold_index"]
        assert "mr_10x" not in new_assembly["info"].columns
        coverage = add_molecule_cov(new_assembly,
                                    scaffolds=scaffolds, binsize=assembly["mol_binsize"], cores=cores)
        old_info = assembly["info"].loc[assembly["info"].scaffold_index.isin(old_to_keep)]
        new_assembly["info"] = pd.concat([old_info, coverage["info"]]).reset_index(drop=True)
        new_assembly["mol_binsize"] = assembly["mol_binsize"]
        old_coverage = assembly["molecule_cov"].loc[assembly["molecule_cov"].scaffold_index.isin(old_to_keep)]
        if coverage["molecule_cov"].shape[0] > 0:
            new_assembly["molecule_cov"] = pd.concat([old_coverage, coverage["molecule_cov"]]).reset_index(drop=True)
        else:
            new_assembly["molecule_cov"] = old_coverage
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
