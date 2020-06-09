from ...utils import rolling_join
import dask.dataframe as dd
import pandas as pd
import os


def _transpose_cssaln(cssaln: str, fai: str, save_dir: str):
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

    fai = dd.read_parquet(fai)
    cssaln = dd.read_parquet(cssaln)
    derived = fai[fai["derived_from_split"] == True]
    to_keep = fai[fai["derived_from_split"] == False]
    to_keep_index = to_keep.index.values.compute()
    assert to_keep_index.shape[0] > 0
    derived = derived[["length", "orig_scaffold_index", "orig_start"]]
    derived = derived.copy().assign(orig_pos=derived["orig_start"])
    assert "orig_start" in derived.columns

    # These are not split
    cols = cssaln.columns[:]
    cssaln_up = cssaln[cssaln["orig_scaffold_index"].isin(to_keep_index)].copy()

    # Now change all broken scaffolds
    _cssaln_down = cssaln[~cssaln["orig_scaffold_index"].isin(to_keep_index)].copy()
    original_length = _cssaln_down.shape[0].compute()
    _cssaln_down = _cssaln_down.drop("length", axis=1).reset_index(
        drop=True).set_index("orig_scaffold_index")
    assert "scaffold_index" in derived.columns or derived.index.name == "scaffold_index"
    cssaln_down = rolling_join(derived, _cssaln_down, on="orig_scaffold_index", by="orig_pos")
    assert cssaln_down.shape[0].compute() >= max(0, original_length), (cssaln_down.shape[0].compute(),
                                                                      original_length)
    assert "scaffold_index" in cssaln_down.columns or cssaln_down.index.name == "scaffold_index", (
        (derived.columns, cssaln_down.columns)
    )
    cssaln_down["pos"] = cssaln_down["orig_pos"] - cssaln_down["orig_start"] + 1
    if original_length > 0:
        assert cssaln_down.shape[0].compute() > 0
    cssaln_down["pos"] = cssaln_down.eval("orig_pos - orig_start + 1")
    cssaln_down = cssaln_down.drop("orig_start", axis=1)
    cssaln_down = cssaln_down.set_index("scaffold_index")
    assert sorted(cssaln_down.columns) == sorted(cssaln_up.columns)
    assert cssaln_down.index.name == cssaln_up.index.name
    cssaln_up = cssaln_up.categorize()
    cssaln_down = cssaln_down.categorize()
    npartitions = cssaln.npartitions
    cssaln = dd.from_pandas(pd.concat([
            cssaln_up.compute(), cssaln_down[cssaln_up.columns].compute()]),
            npartitions=npartitions)
    assert set(cols.values.tolist()) == set(cssaln.columns.values.tolist())

    fname = os.path.join(save_dir, "cssaln")
    dd.to_parquet(cssaln, fname, compression="gzip", engine="pyarrow")
    return fname
