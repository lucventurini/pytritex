import pandas as pd
from ...utils.rolling_join import rolling_join
from ...sequencing_coverage import add_hic_cov
import dask.dataframe as dd
from dask.distributed import Client
import numpy as np
from ...utils import _rebalance_ddf
import os


def _transpose_hic_cov(new_assembly: dict,
                       old_info: dd.DataFrame,
                       fai: dd.DataFrame,
                       save_dir: str,
                       coverage: dd.DataFrame,
                       fpairs: dd.DataFrame,
                       client: Client, cores=1):

    if fpairs is not None and isinstance(fpairs, str):
        fpairs = dd.read_parquet(fpairs)
    else:
        assert isinstance(fpairs, dd.DataFrame) or fpairs is None

    if coverage is not None and isinstance(coverage, str):
        coverage = dd.read_parquet(coverage)
    else:
        assert isinstance(coverage, dd.DataFrame) or coverage is None

    if isinstance(fai, str):
        fai = dd.read_parquet(fai)
    else:
        assert isinstance(fai, dd.DataFrame)

    if isinstance(old_info, str):
        old_info = dd.read_parquet(old_info)
    else:
        assert isinstance(old_info, dd.DataFrame)

    if coverage is not None and fpairs.shape[0].compute() > 0:
        binsize, minNbin, innerDist = new_assembly["binsize"], new_assembly["minNbin"], new_assembly["innerDist"]
        scaffolds = fai.loc[fai["derived_from_split"] == True].index.values.compute()
        old_to_keep = fai.loc[fai["derived_from_split"] == False].index.values.compute()
        # This returns a dictionary with "info" and "cov"
        new_coverage = add_hic_cov(
            new_assembly,
            scaffolds=scaffolds,
            save_dir=save_dir,
            client=client, binsize=binsize, minNbin=minNbin, innerDist=innerDist, cores=cores)

        # First let's get the new coverage
        present = coverage.index.compute()
        inters = present.isin(old_to_keep)
        present = np.unique(present.values[inters])
        previous_to_keep = coverage.loc[present]
        if new_coverage["cov"].shape[0].compute() > 0:
            assert "scaffold_index" == previous_to_keep.index.name == new_coverage["cov"].index.name
            new_assembly["cov"] = dd.concat([
                previous_to_keep.reset_index(drop=False),
                # TODO: for some reason the program crashes here unless I explicitly compute the result
                new_coverage["cov"].compute().reset_index(drop=False)]).set_index(
                "scaffold_index")
        else:
            new_assembly["cov"] = previous_to_keep

        # Now extract the info which is still valid
        present = old_info.index.compute()
        present = np.unique(present.values[present.isin(old_to_keep)])
        old_info = old_info.loc[present].drop("mr_10x", axis=1, errors="ignore")
        if isinstance(new_coverage["info"], str):
            new_coverage["info"] = dd.read_parquet(new_coverage["info"])
        new_assembly["info"] = dd.concat([
            old_info.reset_index(drop=False),
            new_coverage["info"].reset_index(drop=False)
        ]).set_index("scaffold_index")

    return new_assembly


def _transpose_fpairs(fpairs: dd.DataFrame, fai: dd.DataFrame, save_dir: str):
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
    final_columns = ['scaffold_index2', 'orig_scaffold_index2', 'scaffold_index1',
                     'orig_scaffold_index1', 'orig_pos1', 'orig_pos2', 'pos1', 'pos2']

    if fpairs is not None and isinstance(fpairs, str):
        fpairs = dd.read_parquet(fpairs, infer_divisions=True)

    if fpairs is not None and fpairs.shape[0].compute() > 0:
        print("Starting transposition")
        derived = fai[fai["derived_from_split"] == True]
        # "Scaffold_index" is the index name
        left = derived[["orig_scaffold_index", "orig_start"]]
        left = left.assign(orig_pos=left["orig_start"]).reset_index(drop=False)

        bait1 = fpairs.orig_scaffold_index1.isin(derived.orig_scaffold_index.compute())
        bait2 = fpairs.orig_scaffold_index2.isin(derived.orig_scaffold_index.compute())
        stable_pairs = fpairs[~(bait1 | bait2)][final_columns]
        unstable_pairs = fpairs[bait1 & bait2]
        unstable_1_stable_2 = fpairs[bait1 & ~bait2]
        stable_1_unstable_2 = fpairs[~bait1 & bait2]

        # This is the easy one.
        unstable_pairs = unstable_pairs[
            ["orig_scaffold_index1", "orig_scaffold_index2", "orig_pos1", "orig_pos2"]]
        unstable_pairs = rolling_join(left.rename(columns=dict((_, _ + "1") for _ in left.columns)),
                                      unstable_pairs, on="orig_scaffold_index1", by="orig_pos1")
        unstable_pairs = unstable_pairs.reset_index(drop=False)
        assert "scaffold_index1" in unstable_pairs.columns, unstable_pairs.compute().head()
        unstable_pairs = rolling_join(left.rename(columns=dict((_, _ + "2") for _ in left.columns)),
                                      unstable_pairs, on="orig_scaffold_index2", by="orig_pos2")
        assert "scaffold_index2" in unstable_pairs.columns, unstable_pairs.head()
        unstable_pairs["pos1"] = unstable_pairs.eval("orig_pos1 - orig_start1 + 1")
        unstable_pairs["pos2"] = unstable_pairs.eval("orig_pos2 - orig_start2 + 1")
        unstable_pairs = unstable_pairs[final_columns]

        # Unstable 1, stable2: only roll join the scaffold1
        unstable_1_stable_2 = unstable_1_stable_2[
            ["orig_scaffold_index1", "orig_scaffold_index2", "orig_pos1", "orig_pos2",
             'scaffold_index2', 'pos2']]
        unstable_1_stable_2 = rolling_join(left.rename(columns=dict((_, _ + "1") for _ in left.columns)),
                                           unstable_1_stable_2, on="orig_scaffold_index1", by="orig_pos1")
        unstable_1_stable_2["pos1"] = unstable_1_stable_2.eval("orig_pos1 - orig_start1 + 1")
        unstable_1_stable_2 = unstable_1_stable_2.reset_index(drop=False)
        unstable_1_stable_2 = unstable_1_stable_2[final_columns]

        # Stable 1, unstable2: only roll join the scaffold1
        stable_1_unstable_2 = stable_1_unstable_2[
            ["orig_scaffold_index1", "orig_scaffold_index2", "orig_pos1", "orig_pos2",
             'scaffold_index1', 'pos1']]
        stable_1_unstable_2 = rolling_join(left.rename(columns=dict((_, _ + "2") for _ in left.columns)),
                                           stable_1_unstable_2, on="orig_scaffold_index2", by="orig_pos2")
        stable_1_unstable_2["pos2"] = stable_1_unstable_2.eval("orig_pos2 - orig_start2 + 2")
        stable_1_unstable_2 = stable_1_unstable_2[final_columns]

        fpairs = dd.concat([
            stable_pairs,
            unstable_pairs,
            unstable_1_stable_2,
            stable_1_unstable_2]).reset_index(drop=True)
        fpairs = fpairs.astype(dict((_, np.int) for _ in
                                    ["orig_scaffold_index1", "orig_scaffold_index2", "orig_pos1", "orig_pos2",
                                     'scaffold_index1', 'pos1']))
        fpairs = fpairs.persist()
    else:
        fpairs = pd.DataFrame().assign(**dict((column, list()) for column in final_columns))
        fpairs = dd.from_pandas(fpairs, chunksize=1)

    fpairs = _rebalance_ddf(fpairs, target_memory=5 * 10**7).persist()
    fpairs_name = os.path.join(save_dir, "anchored_hic_links")
    dd.to_parquet(fpairs, fpairs_name, compute=True, engine="pyarrow", compression="gzip")
    return fpairs_name
