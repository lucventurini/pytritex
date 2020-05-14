import pandas as pd
import numpy as np


def make_agp(membership: pd.DataFrame, gap_size=100):

    pre_agp = membership.loc[:, ["scaffold", "length", "super", "bin", "orientation"]]
    pre_agp = pre_agp.sort_values(["super", "bin", "length"], ascending=[True, True, False])
    pre_agp.loc[:, "index"] = 2 * pd.Series(range(1, pre_agp.shape[0] + 1)) - 1
    pre_agp.loc[:, "gap"] = False
    pre_agp = pd.concat([
        pre_agp,
            pd.DataFrame({"scaffold": "gap", "gap": True, "super": pre_agp["super"], "bin": np.nan,
                      "length": gap_size, "orientation": np.nan, "index": pre_agp["index"] + 1})
    ])
    # Group and sort by iter. Get all rows except the last one (as it is an extra gap).
    # Also store the number of scaffolds in the super-scaffold bin (in the "n" column)
    agp_grouped = pre_agp.sort_values("index").groupby("super")
    agp = pd.concat([agp_grouped.get_group(index).head(n)
                    for index, n in (agp_grouped.size() - 1).iteritems()])
    # Calculate point of insertion of each scaffold/gap within the super-scaffold
    agp.loc[:, "super_start"] = agp.groupby("super")["length"].shift().cumsum().fillna(0) + 1
    agp.loc[:, "super_end"] = agp.groupby("super")["length"].cumsum()
    #   z[, .(scaffold=scaffold, bed_start=0, bed_end=length,
    #   name=scaffold, score=1, strand=ifelse(is.na(orientation) | orientation == 1, "+", "-"), super=super)]->agp_bed
    agp_bed = pre_agp.loc[:, ["scaffold", "length", "orientation", "super"]]
    agp_bed.loc[:, "bed_start"] = 0
    agp_bed.loc[:, "bed_end"] = agp_bed.loc[:, "length"]
    agp_bed.loc[:, "name"] = agp_bed.loc[:, "scaffold"]
    agp_bed.loc[:, "score"] = 1
    agp_bed.loc[(agp_bed["orientation"].isna()) | (agp_bed["orientation"] == 1), "strand"] = "+"
    agp_bed.loc[agp_bed["strand"].isna(), "strand"] = "-"
    del agp_bed["orientation"]
    agp_bed = agp_bed[["scaffold", "bed_start", "bed_end", "name", "score", "strand", "super"]]
    return {"agp": pre_agp, "agp_bed": agp_bed}
