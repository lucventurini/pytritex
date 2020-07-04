import pandas as pd
import numpy as np
import dask.dataframe as dd


def make_agp(membership: dd.DataFrame, info: dd.DataFrame, gap_size=100):

    if "scaffold_index" not in membership.columns:
        assert "scaffold_index" == membership.index.name
        data = membership[["length", "super", "bin", "orientation", "rank"]].reset_index(drop=False).compute()
    else:
        assert membership.index.name is None
        data = membership[
            ["scaffold_index", "length", "super", "bin", "orientation", "rank"]
        ].reset_index(drop=True).compute()
    data = data.sort_values(["super", "bin", "length"], ascending=[True, True, False])
    data.loc[:, "index"] = pd.Series(np.arange(1, data.shape[0] + 1) * 2 - 1, index=data.index)
    data.loc[:, "gap"] = False
    # Now create the gap rows
    gaps = pd.DataFrame().assign(scaffold_index=np.nan,
                                 super=data["super"].values,
                                 bin=np.nan,
                                 length=gap_size,
                                 orientation=np.nan,
                                 index=data["index"].values + 1)
    # Probably "gap" is a reserved keyword
    gaps.loc[:, "gap"] = True
    data = pd.concat([data, gaps]).sort_values("index")
    # Remove the last row as it is a dangling gap.
    data = data[data.groupby("super").cumcount(ascending=False) > 0]
    # Assign the super_start position
    data.loc[:, "super_start"] = data.groupby("super")["length"].shift(1).transform("cumsum").fillna(0) + 1
    data.loc[:, "super_end"] = data.groupby("super")["length"].transform("cumsum")
    scaffolds = info[["scaffold"]].compute().reset_index(drop=False)
    scaffolds = scaffolds.append({"scaffold_index": np.nan, "scaffold": "gap"}, ignore_index=True)
    data = data.merge(scaffolds, how="left", on="scaffold_index")
    # Now create the BED file
    agp_bed = pd.DataFrame().assign(
        scaffold_index=data["scaffold_index"].values,
        scaffold=data["scaffold"].values, bed_start=0,
        bed_end=data["length"].values, name=data["scaffold"].values,
        score=1, strand=data["orientation"].map({1: "+", -1: "-", np.nan: "+"}).values,
        super=data["super"].values)

    return {"agp": data, "agp_bed": agp_bed}
