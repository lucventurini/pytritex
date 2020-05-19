import pandas as pd
from ...utils import rolling_join


def _transpose_cssaln(cssaln: pd.DataFrame, fai: pd.DataFrame):
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

    cssaln = cssaln.copy().drop("length", axis=1).drop("scaffold_index", axis=1)
    fai = fai[["scaffold_index", "length", "orig_scaffold_index", "orig_start"]].copy().rename(
        columns={"length": "scaffold_length"}).assign(orig_pos=lambda df: df["orig_start"])
    cssaln = rolling_join(fai, cssaln, on="orig_scaffold_index", by="orig_pos")
    cssaln.loc[:, "pos"] = cssaln.eval("orig_pos- orig_start + 1")
    cssaln.drop("orig_start", axis=1, inplace=True)
    return cssaln
