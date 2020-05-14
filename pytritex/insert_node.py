import pandas as pd
import networkit as nk
import numpy as np


def condition_check(df, el):

    """Given a target dataframe (df) and a source dataframe (el) this function
    will check whether there are still clusters which are in bins with distance == 1."""

    bait1 = (el.cluster2.isin(df["cluster"]))
    bait2 = (el.cluster1.isin(
        el.loc[( ~el.cluster1.isin(df["cluster"]) ) & (el.cluster2.isin(df["cluster"])), "cluster1"]
    ))

    __left = el.loc[bait1 & bait2].set_index("cluster2")
    __right = df.loc[:, ["cluster", "bin"]].rename(columns={"cluster": "cluster2"}).set_index("cluster2")
    y = pd.merge(__left, __right, left_index=True, right_index=True, how="outer").loc[
        lambda dataf: ~dataf["cluster1"].isna()].sort_values(["cluster1", "bin"])

    y.loc[:, "dist"] = y.groupby(["cluster1"])["bin"].transform(
        lambda series: np.concatenate([series[1:].values - series[:-1].values, np.array([0])])
        if series.shape[0] > 1 else 0)

    return y


def insert_node(df: pd.DataFrame, el: pd.DataFrame):
    """Function to find the location and insert each additional node in the path.
    Note the "bin" column is ordered from 0 to X."""

    y = condition_check(df, el)
    idx = y["dist"] == 1
    while idx.shape[0] > 0:
        __left = el[["cluster1", "cluster2", "weight"]]
        __left.columns = ["path_node1", "path_node2", "old_path"]
        __left.set_index(["path_node1", "path_node2"])
        __right = pd.DataFrame().assign(
            cluster=y.loc[idx, "cluster1"],
            path_node1=y.loc[idx, "cluster2"],
            path_node2=y.loc[idx + 1, "cluster2"],
            weight=y.loc[idx, "weight"],
            weight2=y.loc[idx + 1, "weight"],
            bin=y.loc[idx, "bin"]).set_index(["path_node1", "path_node2"])
        merged = pd.merge(__left, __right, left_index=True, right_index=True)
        # Now select the ONE node with the least distance
        result = merged.assign(
            diff=merged["weight1"] + merged["weight2"] - merged["old_path"]).sort_values("diff").head(1)
        bin_idx = result["bin"]
        df = pd.concat([
            df.iloc[:bin_idx], pd.DataFrame({"cluster": result["cluster"], "bin": bin_idx + 1}),
            df.iloc[bin_idx:].assign(bin=lambda _df: _df["bin"] + 1)
        ])
        y = condition_check(df, el)
        idx = y["dist"] == 1

    return df
