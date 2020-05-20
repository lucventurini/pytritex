import pandas as pd
from ...utils.rolling_join import rolling_join
from ...sequencing_coverage import add_molecule_cov


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


def _transpose_molecules(molecules, fai):
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

    if molecules is not None and molecules.shape[0] > 0:
        molecules = molecules.copy().drop("scaffold_index", axis=1)
        fai = fai[["scaffold_index", "orig_scaffold_index", "orig_start", "length"]].rename(
            columns={"length": "s_length"}).copy().assign(orig_pos=lambda df: df["orig_start"])
        molecules = rolling_join(fai, molecules, on="orig_scaffold_index", by="orig_start")
        # TODO I think this is wrong, we should switch start and pos
        molecules.loc[:, "start"] = molecules["orig_start"] - molecules["orig_pos"] + 1
        molecules.loc[:, "end"] = molecules["orig_end"] - molecules["orig_pos"] + 1
        molecules = molecules.loc[molecules["end"] <= molecules["s_length"]].drop("orig_pos", axis=1).drop(
            "s_length", axis=1)
        return molecules
    else:
        return pd.DataFrame()
