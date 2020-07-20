import pandas as pd
from ...utils.rolling_join import rolling_join
from ...sequencing_coverage import add_hic_cov
import dask.dataframe as dd
from dask.distributed import Client
import numpy as np
import os
import time
import logging
dask_logger = logging.getLogger("dask")


def _transpose_hic_cov(new_assembly: dict,
                       fai: dd.DataFrame,
                       scaffolds: np.ndarray,
                       old_info: dd.DataFrame,
                       save_dir: str,
                       coverage: dd.DataFrame,
                       fpairs: dd.DataFrame):

    if fpairs is not None and isinstance(fpairs, str):
        fpairs = dd.read_parquet(fpairs)
    else:
        assert isinstance(fpairs, dd.DataFrame) or fpairs is None

    if coverage is not None and isinstance(coverage, str):
        coverage = dd.read_parquet(coverage)
    else:
        assert isinstance(coverage, dd.DataFrame) or coverage is None

    if isinstance(old_info, str):
        old_info = dd.read_parquet(old_info)
    else:
        assert isinstance(old_info, dd.DataFrame)

    if coverage is not None and fpairs.shape[0].compute() > 0:
        binsize, minNbin, innerDist = new_assembly["binsize"], new_assembly["minNbin"], new_assembly["innerDist"]
        # scaffolds = fai.loc[fai["derived_from_split"] == True].index.values.compute()
        # old_to_keep = fai.loc[fai["derived_from_split"] == False].index.values.compute()
        # This returns a dictionary with "info" and "cov"
        new_coverage = add_hic_cov(
            new_assembly,
            scaffolds=None,
            save_dir=save_dir,
            binsize=binsize, minNbin=minNbin, innerDist=innerDist)
        return new_coverage

    return new_assembly


def _transpose_fpairs(fpairs: dd.DataFrame, fai: dd.DataFrame, save_dir: str) -> str:
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
    final_columns = ['scaffold_index1', 'orig_scaffold_index1', 'orig_pos1', 'pos1',
                     'scaffold_index2', 'orig_scaffold_index2', 'orig_pos2', 'pos2']

    if fpairs is not None and isinstance(fpairs, str):
        dask_logger.debug("%s Loading Fpairs from %s", time.ctime(), fpairs)
        fpairs = dd.read_parquet(fpairs, infer_divisions=True)
        assert fpairs.known_divisions is True
    else:
        assert fpairs.known_divisions is True

    if fpairs is not None and fpairs.shape[0].compute() > 0:
        nparts = fpairs.npartitions
        fai_index = fai[fai.to_use == True].index.compute()
        # derived = fai[fai["derived_from_split"] == True]
        # "Scaffold_index" is the index name
        dask_logger.debug("%s Getting the left frame", time.ctime())
        left = fai.loc[(fai.orig_scaffold_index.isin(
            fai.loc[fai["to_use"] == False, "orig_scaffold_index"].values.compute())) & (fai.to_use == True),
                       ["orig_scaffold_index", "orig_start"]]
        left = left.assign(orig_pos=left["orig_start"]).reset_index(drop=False)
        # left = derived[["orig_scaffold_index", "orig_start"]]
        # left = left.assign(orig_pos=left["orig_start"]).reset_index(drop=False)
        dask_logger.debug("%s Getting the baits", time.ctime())
        bait1 = ~fpairs.scaffold_index1.isin(fai_index.values).compute()
        bait2 = ~fpairs.scaffold_index2.isin(fai_index.values).compute()
        
        stable_index = np.intersect1d(bait1[bait1 == False].index.values, bait2[bait2 == False].index.values)
        unstable_index = np.intersect1d(bait1[bait1 == True].index.values, bait2[bait2 == True].index.values)
        unstable1_stable2_index = np.intersect1d(bait1[bait1 == True].index.values, bait2[bait2 == False].index.values)
        stable1_unstable2_index = np.intersect1d(bait1[bait1 == False].index.values, bait2[bait2 == True].index.values)

        stable_pairs = fpairs.loc[stable_index, final_columns].reset_index(drop=False)
        unstable_pairs = fpairs.loc[unstable_index, final_columns].reset_index(drop=False).persist()
        unstable_1_stable_2 = fpairs.loc[unstable1_stable2_index, final_columns].reset_index(drop=False).persist()
        stable_1_unstable_2 = fpairs.loc[stable1_unstable2_index, final_columns].reset_index(drop=False).persist()

        final_columns = [fpairs.index.name] + final_columns

        if unstable_pairs.shape[0].compute() > 0:
            dask_logger.debug("%s Roll join USS", time.ctime())
            partial_columns = ["hic_index", "orig_scaffold_index1", "orig_scaffold_index2",
                               "orig_pos1", "orig_pos2"]
            # This is the easy one.
            dask_logger.debug("%s Roll join USS: 1st. USS: %s; left: %s",
                                time.ctime(), unstable_pairs.shape[0].compute(), left.shape[0].compute())
            unstable_pairs = unstable_pairs[partial_columns]
            unstable_pairs = rolling_join(
                left.rename(columns=dict((_, _ + "1") for _ in left.columns)),
                unstable_pairs, on="orig_scaffold_index1", by="orig_pos1")
            assert "scaffold_index1" in unstable_pairs.columns, unstable_pairs.compute().head()
            dask_logger.debug("%s Roll join USS: 2nd. USS: %s; left: %s",
                                time.ctime(), unstable_pairs.shape[0].compute(), left.shape[0].compute())
            unstable_pairs = rolling_join(
                left.rename(columns=dict((_, _ + "2") for _ in left.columns)),
                unstable_pairs, on="orig_scaffold_index2", by="orig_pos2")
            assert "scaffold_index2" in unstable_pairs.columns, unstable_pairs.head()
            unstable_pairs["pos1"] = unstable_pairs.eval("orig_pos1 - orig_start1 + 1")
            unstable_pairs["pos2"] = unstable_pairs.eval("orig_pos2 - orig_start2 + 1")

        unstable_pairs = unstable_pairs[final_columns]

        # Unstable 1, stable2: only roll join the scaffold1
        if unstable_1_stable_2.shape[0].compute() > 0:
            partial_columns = ["hic_index", "orig_scaffold_index1", "orig_scaffold_index2", "orig_pos1", "orig_pos2",
                               "scaffold_index2", "pos2"]
            dask_logger.debug("%s Roll join US1S2", time.ctime())
            dask_logger.debug("%s Roll join US1S2. US1S2: %s; left: %s",
                                time.ctime(), unstable_1_stable_2.shape[0].compute(), left.shape[0].compute())
            unstable_1_stable_2 = unstable_1_stable_2[partial_columns]
            unstable_1_stable_2 = rolling_join(
                left.rename(columns=dict((_, _ + "1") for _ in left.columns)),
                unstable_1_stable_2, on="orig_scaffold_index1", by="orig_pos1")
            unstable_1_stable_2["pos1"] = unstable_1_stable_2.eval("orig_pos1 - orig_start1 + 1")
        unstable_1_stable_2 = unstable_1_stable_2[final_columns]

        # Stable 1, unstable2: only roll join the scaffold1
        if stable_1_unstable_2.shape[0].compute() > 0:
            partial_columns = ["hic_index", "orig_scaffold_index1", "orig_scaffold_index2", "orig_pos1", "orig_pos2",
                               "scaffold_index1", "pos1"]
            dask_logger.debug("%s Roll join S1US2", time.ctime())
            dask_logger.debug("%s Roll join S1US2. S1US2: %s; left: %s",
                                time.ctime(), stable_1_unstable_2.shape[0].compute(), left.shape[0].compute())
            stable_1_unstable_2 = stable_1_unstable_2[partial_columns]
            stable_1_unstable_2 = rolling_join(
                left.rename(columns=dict((_, _ + "2") for _ in left.columns)),
                stable_1_unstable_2, on="orig_scaffold_index2", by="orig_pos2")
            stable_1_unstable_2["pos2"] = stable_1_unstable_2.eval("orig_pos2 - orig_start2 + 2")
        stable_1_unstable_2 = stable_1_unstable_2[final_columns]

        dask_logger.debug("%s Concatenating fpairs", time.ctime())
        fpairs = dd.concat([
            stable_pairs,
            unstable_pairs.astype(dict(stable_pairs.dtypes)),
            unstable_1_stable_2.astype(dict(stable_pairs.dtypes)),
            stable_1_unstable_2.astype(dict(stable_pairs.dtypes))
        ]).reset_index(drop=True).set_index("hic_index")
        assert fpairs.known_divisions is True
        if fpairs.index.isna().any() is True:
            dd.to_parquet(fpairs, "/tmp/fpairs", compute=True, engine="pyarrow", compression="gzip")
            raise AssertionError
        fpairs = fpairs.astype(dict((_, np.int) for _ in
                                    ["orig_scaffold_index1", "orig_scaffold_index2", "orig_pos1", "orig_pos2",
                                     'scaffold_index1', 'pos1']))
        fpairs = fpairs.repartition(npartitions=nparts)
        dask_logger.debug("%s Finished fpairs", time.ctime())
    else:
        fpairs = pd.DataFrame().assign(**dict((column, list()) for column in final_columns))
        fpairs = dd.from_pandas(fpairs, chunksize=1)

    fpairs_name = os.path.join(save_dir, "fpairs")
    dd.to_parquet(fpairs, fpairs_name, compute=True, engine="pyarrow", compression="gzip")
    return fpairs_name
