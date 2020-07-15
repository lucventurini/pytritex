from ...utils import rolling_join
import dask.dataframe as dd
import numpy as np
from ...utils import _rebalance_ddf
import os
from time import ctime
import logging
dask_logger = logging.getLogger("dask")


def _transpose_cssaln(cssaln: str, fai: dd.DataFrame, save_dir: str) -> str:
    """Transpose the CSS alignment given the newly split scaffolds."""

    #  cat("Transpose cssaln\n")
    #  copy(cssaln) -> z
    #  z[, scaffold_length := NULL]
    #  z[, scaffold := NULL]
    #  fai[, .(scaffold, scaffold_length=length, orig_scaffold, orig_start, orig_pos=orig_start)][
    #  z, on=c("orig_scaffold", "orig_pos"), roll=T]->z
    #  z[, pos := orig_pos - orig_start + 1]
    #  z[, orig_start := NULL]
    #  assembly_new$cssaln <- z

    # fai = dd.read_parquet(fai)
    if isinstance(cssaln, str):
        _cssaln_str = cssaln[:]
        cssaln = dd.read_parquet(cssaln)
    else:
        _cssaln_str = None
        assert isinstance(cssaln, dd.DataFrame)

    css_index = cssaln.index.compute()
    css_up_index = np.unique(
        css_index.intersection(fai.loc[fai["to_use"] == True].index.compute().values).values)
    cssaln_up = cssaln.loc[css_up_index].copy()
    css_down_index = css_index.difference(css_up_index).values
    _cssaln_down = cssaln.loc[css_down_index].copy()

    derived = fai.loc[
                  (fai.orig_scaffold_index.isin(fai.loc[css_down_index, "orig_scaffold_index"].values.compute())) & (
                          fai["to_use"] == True), ["orig_scaffold_index", "orig_start", "length"]][:]

    derived = derived.copy().assign(orig_pos=derived["orig_start"])

    dask_logger.warning("To keep: %s; derived: %s", css_up_index.shape[0], derived.shape[0].compute())
    # assert np.isnan(cssaln.orig_scaffold_index.values.compute()).any() == False, _cssaln_str

    _cssaln_down = _cssaln_down.drop("length", axis=1).reset_index(drop=True).set_index("orig_scaffold_index")
    # assert np.isnan(_cssaln_down.index.values.compute()).any() == False
    assert derived.index.isna().any().compute() == False and derived.index.name == "scaffold_index"
    cssaln_down = rolling_join(derived.reset_index(drop=False),
                               _cssaln_down, on="orig_scaffold_index", by="orig_pos")
    assert np.isnan(cssaln_down.scaffold_index.values.compute()).any() == False, \
        cssaln_down[cssaln_down["scaffold_index"].isna()].head(npartitions=-1)
    assert np.isnan(cssaln_down.index.values.compute()).any() == False
    assert "scaffold_index" in cssaln_down.columns or cssaln_down.index.name == "scaffold_index", (
        (derived.columns, cssaln_down.columns)
    )
    assert cssaln_down[cssaln_down["scaffold_index"].isna()].shape[0].compute() == 0

    cssaln_down["pos"] = cssaln_down.eval("orig_pos - orig_start + 1")
    cssaln_down = cssaln_down.drop("orig_start", axis=1)
    cssaln_down = cssaln_down.astype(
        {"scaffold_index": cssaln_up.index.dtype}).set_index("scaffold_index")
    cssaln = dd.concat(
        [cssaln_up,
         cssaln_down[cssaln_up.columns]])
    cssaln = cssaln.reset_index(drop=False).set_index("scaffold_index")
    # assert set(cols.values.tolist()) == set(cssaln.columns.values.tolist())
    cssaln = cssaln.persist()
    dask_logger.debug("%s Finished calculating CSS-ALN, rebalancing", ctime())

    cssaln = _rebalance_ddf(cssaln, target_memory=5 * 10**7)
    dask_logger.debug("%s Rebalanced CSS-ALN", ctime())
    cssaln_name = os.path.join(save_dir, "cssaln")
    dd.to_parquet(cssaln, cssaln_name, engine="pyarrow", compression="gzip", compute=True)
    dask_logger.warning("%s Saved CSS-ALN in %s", ctime(), cssaln_name)
    return cssaln_name
