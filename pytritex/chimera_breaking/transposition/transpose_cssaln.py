import pandas as pd
from ...utils import rolling_join
import numexpr as ne


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

    derived = fai.loc[fai.derived_from_split]
    to_keep = fai.loc[~fai.orig_scaffold_index.isin(derived.orig_scaffold_index)]
    derived = derived[["scaffold_index", "length", "orig_scaffold_index", "orig_start"]].copy().assign(
        orig_pos=lambda df: df["orig_start"])

    # These are not split
    cols = cssaln.columns[:]
    cssaln_up = cssaln.copy().loc[cssaln["orig_scaffold_index"].isin(to_keep["orig_scaffold_index"])]
    cssaln_down= cssaln.copy().loc[~cssaln["orig_scaffold_index"].isin(to_keep["orig_scaffold_index"])]
    cssaln_down = cssaln_down.copy().drop("length", axis=1).drop("scaffold_index", axis=1)
    cssaln_down = rolling_join(derived, cssaln_down, on="orig_scaffold_index", by="orig_pos")
    cssaln = pd.concat([cssaln_up, cssaln_down]).sort_values(["scaffold_index", "orig_pos"]).reset_index(drop=True)
    cssaln.loc[:, "pos"] = cssaln.eval("orig_pos - orig_start + 1")
    cssaln.drop("orig_start", axis=1, inplace=True)
    print("Original css columns:", cols)
    print("New css columns:", cssaln.columns)
    return cssaln
