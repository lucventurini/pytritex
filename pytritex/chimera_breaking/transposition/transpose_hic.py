import pandas as pd
from ...utils.rolling_join import rolling_join
from ...sequencing_coverage import add_hic_cov


def _transpose_hic_cov(new_assembly, old_info, fai, coverage, fpairs, cores=1):
    if coverage is not None and fpairs.shape[0] > 0:
        binsize, minNbin, innerDist = new_assembly["binsize"], new_assembly["minNbin"], new_assembly["innerDist"]
        scaffolds = fai.loc[fai["derived_from_split"], "scaffold_index"]
        old_to_keep = fai.loc[~fai["derived_from_split"], "scaffold_index"]
        new_coverage = add_hic_cov(new_assembly, scaffolds=scaffolds,
                    binsize=binsize, minNbin=minNbin, innerDist=innerDist, cores=cores)
        previous_to_keep = coverage[coverage["scaffold_index"].isin(old_to_keep)]
        if new_coverage["cov"].shape[0] > 0:
            new_assembly["cov"] = pd.concat([previous_to_keep, new_coverage["cov"]])
        else:
            new_assembly["cov"] = previous_to_keep
        old_info = old_info.loc[~new_assembly["info"]["scaffold_index"].isin(old_to_keep)].drop(
            "mr_10x", axis=1, errors="ignore")
        new_assembly["info"] = pd.concat([
            old_info, new_coverage["info"]
        ])
    else:
        new_assembly["fpairs"] = pd.DataFrame()
    return new_assembly


def _transpose_fpairs(fpairs, fai):
    # if("fpairs" %in% names(assembly) && nrow(fpairs) > 0){
    #   cat("Transpose fpairs\n")
    #   assembly$fpairs[, .(orig_scaffold1, orig_scaffold2, orig_pos1, orig_pos2)]->z
    #   fai[, .(scaffold1=scaffold, orig_scaffold1=orig_scaffold, orig_start1=orig_start, orig_pos1=orig_start)
    #   ][z, on=c("orig_scaffold1", "orig_pos1"), roll=T]->z
    #   fai[, .(scaffold2=scaffold, orig_scaffold2=orig_scaffold, orig_start2=orig_start, orig_pos2=orig_start)][
    #   z, on=c("orig_scaffold2", "orig_pos2"), roll=T]->z
    #   z[, pos1 := orig_pos1 - orig_start1 + 1]
    #   z[, pos2 := orig_pos2 - orig_start2 + 1]
    #   z[, orig_start1 := NULL]
    #   z[, orig_start2 := NULL]
    #   assembly_new$fpairs <- z
    #  } else {
    #   assembly_new$fpairs <- data.table()
    #  }
    print("Transpose fpairs")
    if fpairs is not None and fpairs.shape[0] > 0:
        try:
            pairs = fpairs[["orig_scaffold_index1", "orig_scaffold_index2", "orig_pos1", "orig_pos2"]]
        except KeyError:
            raise KeyError(fpairs.columns)
        left = fai[["scaffold_index", "orig_scaffold_index", "orig_start"]]
        left = left.assign(orig_pos=left["orig_start"])
        pairs = rolling_join(left.rename(columns=dict((_, _ + "1") for _ in left.columns)),
                             pairs, on="orig_scaffold_index1", by="orig_pos1")
        pairs = rolling_join(left.rename(columns=dict((_, _ + "2") for _ in left.columns)),
                             pairs, on="orig_scaffold_index2", by="orig_pos2")
        pairs.loc[:, "pos1"] = pairs["orig_pos1"] - pairs["orig_start1"] + 1
        pairs.loc[:, "pos2"] = pairs["orig_pos2"] - pairs["orig_start2"] + 1
        pairs = pairs.drop("orig_start1", axis=1).drop("orig_start2", axis=1)
        return pairs
    else:
        return pd.DataFrame()
