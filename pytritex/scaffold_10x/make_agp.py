import pandas as pd
import numpy as np
import dask.dataframe as dd
from pytritex.utils.chrnames import chrNames
import dask.array as da


def make_agp(membership: dd.DataFrame, info: dd.DataFrame, names=None, gap_size=100):

    # names = chrNames()
    if "scaffold_index" not in membership.columns:
        assert "scaffold_index" == membership.index.name
        initial = membership[
            ["length", "super", "bin", "orientation", "rank", "chr", "cM"]].reset_index(drop=False)
    else:
        assert membership.index.name is None or membership.index.name == "index", membership.index.name
        initial = membership[
            ["scaffold_index", "length", "super", "bin", "orientation", "rank", "chr", "cM"]
        ].reset_index(drop=True)
    if names is None:
        names = chrNames()
    initial = dd.from_pandas(initial.compute().sort_values(["super", "bin", "length"], ascending=[True, True, False]).reset_index(drop=True),
                             chunksize=int(1e6))
    initial = initial.merge(names, on="chr", how="left")
    initial["index"] = da.from_array(np.arange(1, initial.shape[0].compute() + 1) * 2 - 1,
                                            chunks=initial.map_partitions(len).compute().values.tolist())
    initial = initial.set_index("index")
    initial["gap"] = False
    initial["start"] = 1

    # Now create the gap rows
    gaps = pd.DataFrame().assign(super=initial["super"].values.compute(),
                                 bin=np.nan,
                                 start=1,  # This is the "gap_type" column
                                 length=gap_size,
                                 orientation=np.nan,
                                 index=initial.index.values.compute() + 1,
                                 alphachr=np.nan,
                                 cM=np.nan,
                                 scaffold_index=-1,
                                 gap=True,
                                 with_gaps=(initial.groupby("super")["scaffold_index"].transform("size", meta=int) > 1).values.compute()).query(
        "with_gaps == True")[:].drop("with_gaps", axis=1).set_index("index")
    # Probably "gap" is a reserved keyword
    assert (gaps["scaffold_index"].unique() == [-1]).all()
    data = dd.concat([initial, gaps]).reset_index(drop=False).set_index("index")
    # Assign the super_start position
    cusum = data.groupby("super")["length"].transform("cumsum", meta=int).compute()
    data["super_start"] = da.from_array((cusum - data.groupby("super")["length"].apply(lambda x: x.shift(0), meta=int).values.compute() + 1).values,
                                        chunks=data.map_partitions(len).compute().values.tolist())
    data["super_end"] = da.from_array((cusum).values, chunks=data.map_partitions(len).compute().values.tolist())
    # assign the final index
    data["super_index"] = da.from_array(data.groupby("super").apply(lambda g: list(range(1, g.shape[0] + 1)), meta=list).compute().explode().values,
                                        chunks=data.map_partitions(len).compute().values.tolist())
    original = info[info.derived_from_split == False][["scaffold", "length"]].rename(
        columns={"scaffold": "original_scaffold", "length": "orig_length"}).reset_index(drop=False)
    scaffolds = info[["orig_scaffold_index", "scaffold", "orig_start", "orig_end"]]
    # We need to isolate the scaffold_index column or it will be rewritten over during the merge.
    scaffolds["scaffold_key"] = scaffolds.index.values
    scaffolds = original.merge(scaffolds, left_on="scaffold_index", right_on="orig_scaffold_index", how="right")
    scaffolds["scaffold_index"] = scaffolds["scaffold_key"].values
    scaffolds = scaffolds.drop("scaffold_key", axis=1)
    assert scaffolds.query("scaffold_index != scaffold_index").shape[0].compute() == 0
    scaffolds = scaffolds.set_index("scaffold_index")
    assert scaffolds.query("orig_scaffold_index == scaffold_index & scaffold != original_scaffold").shape[0].compute() == 0
    gap = {"scaffold_index": [-1], "scaffold": ["gap"], "orig_start": [1], "orig_end": [gap_size],
           "original_scaffold": ["gap"], "orig_length": [gap_size], "orig_scaffold_index": [-1]}
    gap = pd.DataFrame(data=gap).set_index("scaffold_index")[scaffolds.columns]
    assert gap.loc[-1, "orig_end"] == gap_size
    assert gap.loc[-1, "orig_length"] == gap_size
    assert gap.loc[-1, "orig_start"] == 1
    gap = dd.from_pandas(gap, chunksize=1)
    scaffolds = dd.concat([gap, scaffolds], axis="index", interleave_partitions=True).persist()
    if scaffolds.index.name is None:
        assert "scaffold_index" in scaffolds.columns, (scaffolds.columns, scaffolds.head(npartitions=-1)) 
        scaffolds = scaffolds.set_index("scaffold_index")
    # try:
    #     print(scaffolds.loc[-1].compute())
    # except KeyError:
    #     print(scaffolds.head(npartitions=-1))
    #     raise    
    assert scaffolds.loc[-1, "orig_end"].compute().values[0] == gap_size, scaffolds.loc[-1, "orig_end"].compute() 
    assert scaffolds.loc[-1, "orig_length"].compute().values[0] == gap_size, scaffolds.loc[-1, "orig_length"].compute() 
    assert scaffolds.loc[-1, "orig_start"].compute().values[0] == 1, scaffolds.loc[-1, "orig_start"].compute().values[0]
    
    # import sys
    # sys.exit(1)
    assert scaffolds.query("orig_length != orig_length").shape[0].compute() == 0
    data = data.reset_index(drop=False).set_index("scaffold_index").merge(
        scaffolds, how="left", on="scaffold_index")
    assert data.index.name == "scaffold_index"
    assert "index" in data.columns
    data = data.reset_index(drop=False).set_index("index").persist()
    # data["orig_end"] = data["orig_end"].mask(data["scaffold_index"] == -1, gap_size)
    # data["orig_start"] = data["orig_start"].fillna(0).astype(int)
    agp = data[["super", "super_start", "super_end", "super_index", "gap",
                "original_scaffold", "orig_start", "orig_end", "orientation", "alphachr", "cM",
                "scaffold_index", "length", "orig_length"]][:].persist()  # We'll remove the length later
    
            
    for column in ["super_start", "super_end", "super_index", "orig_start", "orig_end", "length", "orig_length"]:
        qstring = f"{column} != {column}"
        nulls = agp.query(qstring).compute()
        if nulls.shape[0] > 0:
            import sys
            print(column, file=sys.stderr)
            print(nulls.shape[0], file=sys.stderr)
            print(nulls.head(), file=sys.stderr)
            sys.exit(1)        
        try:
            agp[column] = agp[column].astype(int)
        except ValueError:
            raise ValueError(column)

    assert agp.query("super_index != super_index").shape[0].compute() == 0
    agp = agp.drop("orig_length", axis=1).persist()
    agp = agp.persist()
    print(agp.head())
    # agp["super"] = agp["super"].astype(str)
    # agp.loc[:, "super"] = "super_" + agp["super"].astype(str)
    # assert agp["super"].dtype == "object", agp["super"].head()
    # agp.loc[data["gap"] == True, "orig_start"] = "scaffold"  # This is the "gap_type" column
    # agp.loc[data["gap"] == True, "orig_end"] = "yes"  # This is the "linkage" column
    # agp["gap"] = agp["gap"].map({True: "U", False: "W"})
    # agp.loc[:, "orientation"] = agp["orientation"].map({1: "+", -1: "-", "paired-ends;map": "paired-ends;map"})
    # agp.loc[agp["ssize"] == 1, "super"] = agp.loc[agp["ssize"] == 1, "original_scaffold"]
    # bait = (agp["ssize"] == 1) & (((agp["orig_start"] != "scaffold") & (agp["orig_start"] != 1)) |
    #                               ((agp["orig_end"] != agp["orig_length"]) & (agp["orig_end"] != "yes")))
    # agp.loc[bait, "super"] = agp.loc[bait, ["original_scaffold", "orig_start", "orig_end"]].astype(str).apply(
    #     lambda row: row["original_scaffold"] + ":" + row["orig_start"] + "-" + row["orig_end"], axis=1)

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
    agp_bed = agp_bed.persist()

    return {"agp": agp, "agp_bed": agp_bed}
