import pandas as pd
import numpy as np
import dask.dataframe as dd
from pytritex.utils.chrnames import chrNames


def make_agp(membership: dd.DataFrame, info: dd.DataFrame, names=None, gap_size=100):

    # names = chrNames()
    if "scaffold_index" not in membership.columns:
        assert "scaffold_index" == membership.index.name
        initial = membership[
            ["length", "super", "bin", "orientation", "rank", "chr", "cM"]].reset_index(drop=False).compute()
    else:
        assert membership.index.name is None or membership.index.name == "index", membership.index.name
        initial = membership[
            ["scaffold_index", "length", "super", "bin", "orientation", "rank", "chr", "cM"]
        ].reset_index(drop=True).compute()
    if names is None:
        names = chrNames()
    initial = initial.sort_values(["super", "bin", "length"], ascending=[True, True, False]).reset_index(drop=True)
    initial = initial.merge(names, on="chr", how="left")
    initial.loc[:, "index"] = pd.Series(np.arange(1, initial.shape[0] + 1) * 2 - 1, index=initial.index)
    initial.loc[:, "gap"] = False
    initial.loc[:, "start"] = 1

    # Now create the gap rows
    gaps = pd.DataFrame().assign(super=initial["super"].values,
                                 bin=np.nan,
                                 start=1,  # This is the "gap_type" column
                                 length=gap_size,
                                 orientation=np.nan,
                                 index=initial["index"].values + 1,
                                 alphachr=np.nan,
                                 cM=np.nan,
                                 scaffold_index=-1,
                                 gap=True,
                                 with_gaps=(initial.groupby("super")["scaffold_index"].transform("size") > 1).values).query(
        "with_gaps == True")[:].drop("with_gaps", axis=1)
    # Probably "gap" is a reserved keyword
    assert (gaps["scaffold_index"].unique() == [-1]).all()
    data = pd.concat([initial, gaps]).sort_values("index").reset_index(drop=True)
    # Assign the super_start position
    cusum = data.groupby("super")["length"].transform("cumsum")
    data.loc[:, "super_start"] = (cusum - data.groupby("super")["length"].shift(0) + 1).values
    data.loc[:, "super_end"] = (cusum).values
    # assign the final index
    data["super_index"] = data.groupby("super").apply(lambda g: range(1, g.shape[0] + 1)).explode().values
    original = info[info.derived_from_split == False][["scaffold", "length"]].rename(
        columns={"scaffold": "original_scaffold", "length": "orig_length"}).compute().reset_index(drop=False)
    scaffolds = info[["orig_scaffold_index", "scaffold", "orig_start", "orig_end"]].compute()
    # We need to isolate the scaffold_index column or it will be rewritten over during the merge.
    scaffolds["scaffold_key"] = scaffolds.index.values
    scaffolds = original.merge(scaffolds, left_on="scaffold_index", right_on="orig_scaffold_index", how="right")
    scaffolds["scaffold_index"] = scaffolds["scaffold_key"].values
    scaffolds = scaffolds.drop("scaffold_key", axis=1)
    assert scaffolds.query("orig_scaffold_index == scaffold_index & scaffold != original_scaffold").shape[0] == 0
    scaffolds = scaffolds.append({"scaffold_index": -1, "scaffold": "gap", "orig_start": 1,
                                  "orig_end": gap_size, "original_scaffold": "gap", "orig_length": gap_size},
                                 ignore_index=True)
    data = data.merge(scaffolds, how="left", on="scaffold_index")
    data.loc[:, "orig_start"] = data["orig_start"].fillna(0).astype(int)
    agp = data[["super", "super_start", "super_end", "super_index", "gap",
                "original_scaffold", "orig_start", "orig_end", "orientation", "alphachr", "cM",
                "scaffold_index", "length", "orig_length"]][:]  # We'll remove the length later
    nulls = agp.loc[data["orig_end"].isna()]
    if nulls.shape[0] > 0:
        print(nulls.head())
        import sys
        sys.exit(1)
    
    for column in ["super_start", "super_end", "orig_start", "orig_end", "length", "orig_length"]:
        try:
            agp.loc[:, column] = agp[column].astype(int)
        except ValueError:
            raise ValueError(column)
    
    # agp["super"] = agp["super"].astype(str)
    # agp.loc[:, "super"] = "super_" + agp["super"].astype(str)
    # assert agp["super"].dtype == "object", agp["super"].head()
    # agp.loc[data["gap"] == True, "orig_start"] = "scaffold"  # This is the "gap_type" column
    # agp.loc[data["gap"] == True, "orig_end"] = "yes"  # This is the "linkage" column
    agp.loc[:, "gap"] = agp["gap"].map({True: "U", False: "W"})
    # agp.loc[:, "orientation"] = agp["orientation"].map({1: "+", -1: "-", "paired-ends;map": "paired-ends;map"})
    # agp.loc[agp["ssize"] == 1, "super"] = agp.loc[agp["ssize"] == 1, "original_scaffold"]
    # bait = (agp["ssize"] == 1) & (((agp["orig_start"] != "scaffold") & (agp["orig_start"] != 1)) |
    #                               ((agp["orig_end"] != agp["orig_length"]) & (agp["orig_end"] != "yes")))
    # agp.loc[bait, "super"] = agp.loc[bait, ["original_scaffold", "orig_start", "orig_end"]].astype(str).apply(
    #     lambda row: row["original_scaffold"] + ":" + row["orig_start"] + "-" + row["orig_end"], axis=1)
    agp = agp.drop("orig_length", axis=1)
    # z[, .(scaffold=scaffold, bed_start=0, bed_end=scaffold_length,
    # name=scaffold, score=1, strand=ifelse(is.na(orientation) | orientation == 1, "+", "-"), agp_chr=agp_chr)]->agp_bed
    # agp = agp.rename(columns={"original_scaffold": "scaffold"})

    agp_bed = agp[:].query("original_scaffold != 'gap'")[["super", "original_scaffold",
                                                          "orig_start", "orig_end", "orientation", "alphachr"]]
    agp_bed = agp_bed.rename(columns={"orig_start": "bed_start", "orig_end": "bed_end", "orientation": "strand",
                                      "original_scaffold": "scaffold"})
    agp_bed["score"] = 1
    agp_bed["bed_start"] -= 1
    # agp_bed["name"] = agp_bed["original_scaffold"].astype(str) + ":" + agp_bed["bed_start"].astype(str) + "-" + agp_bed["bed_end"].astype(str)
    agp_bed["strand"] = agp_bed["strand"].mask(agp_bed["strand"] == 1, "+")
    agp_bed["strand"] = agp_bed["strand"].mask(agp_bed["strand"] == -1, "-")

    return {"agp": agp, "agp_bed": agp_bed}
