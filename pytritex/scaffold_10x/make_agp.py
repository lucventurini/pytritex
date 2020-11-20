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
                                 chr=np.nan,
                                 alphachr=np.nan,
                                 cM=np.nan,
                                 scaffold_index=-1,
                                 gap=True)
    with_gaps = (initial.groupby("super").size() > 1).compute().to_frame("with_gaps")
    gaps = gaps.merge(with_gaps, on="super", how="left")
    gaps = gaps.query("with_gaps == True")[:].drop("with_gaps", axis=1).set_index("index")
    # Probably "gap" is a reserved keyword
    assert (gaps["scaffold_index"].unique() == [-1]).all()
    super_size = initial.groupby("super").size().to_frame("super_size").compute()
    data = dd.concat([initial, gaps]).reset_index(drop=False).set_index("index")
    data = data.persist()
    out = open("/tmp/outer.txt", "wt")
    print(data.head(), file=out) 
    # data["new_index"] = da.from_array(np.arange(1, data.shape[0].compute() + 1),
    #                               chunks=data.map_partitions(len).compute().values.tolist())
    # data = data.set_index("new_index")
    # data.index = data.index.rename("index")    
    data = data.reset_index(drop=False).set_index("super").merge(super_size, on="super", how="left").persist()    
    data = data.reset_index(drop=False).set_index("index")
    # Assign the super_start position
    
    data["super_end"] = data.groupby("super")["length"].cumsum().compute()
    data["super_start"] = data["super_end"] - data["length"] + 1    
    # assign the final index
    data["super_index"] = (data.groupby("super").cumcount() + 1)    
    data.persist()
    assert "super_size" in data.columns, data.head().to_csv("/tmp/broken.csv", sep="\t")
    scaffolds = info.query("to_use == True")[["orig_scaffold_index", "orig_start", "orig_end", "length"]]
    scaffolds = scaffolds.rename(columns={"length": "orig_length"})
    
    gap = {"scaffold_index": [-1], "orig_start": [1], "orig_end": [gap_size], "orig_length": [gap_size],
           "orig_scaffold_index": [-1]}
    gap = pd.DataFrame(data=gap).set_index("scaffold_index")[scaffolds.columns]
    assert gap.loc[-1, "orig_end"] == gap_size
    assert gap.loc[-1, "orig_length"] == gap_size
    assert gap.loc[-1, "orig_start"] == 1
    scaffolds = dd.concat([gap, scaffolds], axis="index", interleave_partitions=True).persist()
    scaffolds = scaffolds.repartition(npartitions=max(1, int(scaffolds.shape[0].compute() // 5e5)))
    assert scaffolds.loc[-1, "orig_end"].compute().values[0] == gap_size, scaffolds.loc[-1, "orig_end"].compute()
    assert scaffolds.loc[-1, "orig_length"].compute().values[0] == gap_size, scaffolds.loc[-1, "orig_length"].compute()
    gstart = scaffolds.loc[-1, "orig_start"].compute().values[0]
    assert gstart == 1, (gstart, type(gstart))
    assert scaffolds.query("orig_length != orig_length").shape[0].compute() == 0

    data = data.reset_index(drop=False).set_index("scaffold_index").merge(scaffolds, how="left", on="scaffold_index")
    assert data.index.name == "scaffold_index"
    assert "index" in data.columns
    data = data.reset_index(drop=False).set_index("super").persist()
    
    # data["orig_end"] = data["orig_end"].mask(data["scaffold_index"] == -1, gap_size)
    # data["orig_start"] = data["orig_start"].fillna(0).astype(int)
    # map_10x["agp"]["super_size"] = da.from_array(map_10x["agp"].groupby("super")["super_index"].transform("max", meta=int).compute().values,
    #                                                  chunks=map_10x["agp"].map_partitions(len).compute().values.tolist())
    #     map_10x["agp"]["super_name"] = map_10x["agp"]["original_scaffold"].mask(
    #         map_10x["agp"]["super_size"] > 1, "super_" + map_10x["agp"]["super"].astype(str))
    #     map_10x["agp"]["super_length"] = da.from_array(map_10x["agp"].groupby("super")["length"].transform("sum", meta=int).compute().values,
    #                                                    chunks=map_10x["agp"].map_partitions(len).compute().values.tolist())

    data = data.reset_index(drop=False).set_index("super").merge(
        data.groupby("super")["length"].sum().to_frame("super_length"), on="super", how="left")
    data = data.reset_index(drop=False).set_index("index").persist()
    # data["super_name"] = data["original_scaffold"].mask(data["super_size"] > 1, "super_" + data["super"].astype(str))

    agp = data[["super", "super_start", "super_end", "super_index", "gap",
                "orig_scaffold_index", "orig_start", "orig_end", "orientation", "alphachr", "cM",
                "scaffold_index", "length", "super_length", "super_size"]][:].persist()

    scaffolds = info.query("to_use == True")[["scaffold", "orig_scaffold_index"]].reset_index(drop=False)
    scaffolds = scaffolds.set_index("orig_scaffold_index")
    orig_scaffolds = info.query("derived_from_split == False")[["scaffold", "length"]].rename(columns={"scaffold": "orig_scaffold"})
    orig_scaffolds = orig_scaffolds.rename(columns={"length": "orig_scaffold_length"})
    orig_scaffolds.index = orig_scaffolds.index.rename("orig_scaffold_index")
    scaffolds = scaffolds.merge(orig_scaffolds, on="orig_scaffold_index", how="left").compute()
    scaffolds = scaffolds.reset_index(drop=True).set_index("scaffold_index")

    agp = agp.merge(scaffolds, on="scaffold_index", how="left")
    agp["super_name"] = agp["scaffold"].mask(agp["super_size"] > 1, "super_" + agp["super"].astype(str))    
    agp = agp.persist()
    agp = agp.drop("scaffold", axis=1).rename(columns={"orig_scaffold": "scaffold"})
    agp["scaffold"] = agp["scaffold"].mask(agp["scaffold_index"] == -1, "gap")    
    agp = agp.persist()
    
    for column in ["super_start", "super_end", "super_index", "orig_start", "orig_end", "length"]:
        qstring = f"{column} != {column}"
        nulls = agp.query(qstring).compute()
        if nulls.shape[0] > 0:
            h = nulls[["super", "scaffold_index", "scaffold", "orig_start"]].head()
            err = f"""NA values found in {column}
{h}"""            
            raise ValueError(err)
        try:
            agp[column] = agp[column].astype(int)
        except ValueError:
            raise ValueError(column)

    agp = agp.persist()
    agp["orig_scaffold_length"] = agp["orig_scaffold_length"].mask(agp["scaffold_index"] == -1, gap_size).astype(int)
    agp = agp.persist()
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

    agp_bed = agp[:].query("scaffold_index != -1")[["super", "orig_scaffold_index", "orig_start",
                                                    "orig_end", "orientation", "alphachr"]]
    agp_bed = agp_bed.rename(columns={"orig_start": "bed_start", "orig_end": "bed_end", "orientation": "strand",
                                      "orig_scaffold_index": "scaffold"})
    agp_bed["score"] = 1
    agp_bed["bed_start"] -= 1
    # agp_bed["name"] = agp_bed["original_scaffold"].astype(str) + ":" + agp_bed["bed_start"].astype(str) + "-" + agp_bed["bed_end"].astype(str)
    agp_bed["strand"] = agp_bed["strand"].mask(agp_bed["strand"] == 1, "+")
    agp_bed["strand"] = agp_bed["strand"].mask(agp_bed["strand"] == -1, "-")
    agp_bed = agp_bed.persist()

    return {"agp": agp, "agp_bed": agp_bed}
