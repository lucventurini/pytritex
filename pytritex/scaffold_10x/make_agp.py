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
    initial = initial.reset_index(drop=False).set_index("index")
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
    data = dd.concat([initial, gaps]).reset_index(drop=False).set_index("index")
    # Assign the super_start position
    data["super_start"] = (data.groupby("super")["length"].cumsum().shift(1).fillna(0) + 1).astype(int)
    data["super_end"] = data.groupby("super")["length"].cumsum().compute()
    # assign the final index
    data["super_index"] = (data.groupby("super").cumcount() + 1)

    scaffolds = info.query("to_use == True")[["orig_scaffold_index", "orig_start", "orig_end"]]

    gap = {"scaffold_index": [-1], "scaffold": ["gap"], "orig_start": [1], "orig_end": [gap_size],
           "original_scaffold": ["gap"], "orig_length": [gap_size], "orig_scaffold_index": [-1]}
    gap = pd.DataFrame(data=gap).set_index("scaffold_index")[scaffolds.columns]
    assert gap.loc[-1, "orig_end"] == gap_size
    assert gap.loc[-1, "orig_length"] == gap_size
    assert gap.loc[-1, "orig_start"] == 1
    scaffolds = dd.concat([gap, scaffolds], axis="index", interleave_partitions=True).persist()
    scaffolds = scaffolds.repartition(npartitions=max(1, int(scaffolds.shape[0].compute() // 5e5)))
    assert scaffolds.loc[-1, "orig_end"].compute().values[0] == gap_size, scaffolds.loc[-1, "orig_end"].compute()
    assert scaffolds.loc[-1, "orig_length"].compute().values[0] == gap_size, scaffolds.loc[-1, "orig_length"].compute() 
    assert scaffolds.loc[-1, "orig_start"].compute().values[0] == 1, scaffolds.loc[-1, "orig_start"].compute().values[0]
    assert scaffolds.query("orig_length != orig_length").shape[0].compute() == 0

    data = data.reset_index(drop=False).set_index("scaffold_index").merge(scaffolds, how="left", on="scaffold_index")
    assert data.index.name == "scaffold_index"
    assert "index" in data.columns
    data = data.reset_index(drop=False).set_index("super").persist()
    data = data.merge(data.groupby("super").size().to_frame("super_size"), on="super", how="left")
    data = data.merge(data.groupby("super")["length"].sum().to_frame("super_length"), on="super", how="left")
    data = data.reset_index(drop=False).set_index("index")

    agp = data[["super", "super_start", "super_end", "super_index", "gap",
                "orig_scaffold_index", "orig_start", "orig_end", "orientation", "alphachr", "cM",
                "scaffold_index", "length", "super_length", "super_size"]][:].persist()

    for column in ["super_start", "super_end", "super_index", "orig_start", "orig_end", "length"]:
        qstring = f"{column} != {column}"
        nulls = agp.query(qstring).compute()
        if nulls.shape[0] > 0:
            raise ValueError(f"NA values found in {column}")
        try:
            agp[column] = agp[column].astype(int)
        except ValueError:
            raise ValueError(column)

    agp = agp.persist()
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
