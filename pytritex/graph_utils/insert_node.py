import pandas as pd
import numpy as np


def condition_check(df, el):

    """Given a target dataframe (df) and a source dataframe (el) this function
    will check whether there are still clusters which are in bins with distance == 1."""

    bait1 = (el.cluster2.isin(df["cluster"]))
    bait2 = (el.cluster1.isin(
        el.loc[(~el.cluster1.isin(df["cluster"])) & (el.cluster2.isin(df["cluster"])), "cluster1"]
    ))

    __left = el.loc[bait1 & bait2].set_index("cluster2")
    __right = df.loc[:, ["cluster", "bin"]].rename(columns={"cluster": "cluster2"}).set_index("cluster2")
    # The "cluster1 != cluster1" is a NumPY idiosyncrasy - NaN is the only value which is
    # always *different* from itself.
    y = pd.merge(__left, __right, left_index=True,
                 right_index=True, how="outer").query("cluster1 == cluster1").sort_values(["cluster1", "bin"])
    # TODO this seems less than ideal.
    grouped = y.groupby(["cluster1"])
    y.loc[:, "dist"] = (grouped["bin"].shift(-1) - grouped["bin"].transform(lambda s: s)).fillna(0)
    y = y.reset_index(drop=False)
    # y.loc[:, "dist"] = grouped["bin"].transform(
    #     lambda series: np.concatenate([series[1:].values - series[:-1].values, np.array([0])])
    #     if series.shape[0] > 1 else 0)
    return y


def insert_node(df: pd.DataFrame, el: pd.DataFrame):
    """Function to find the location and insert each additional node in the path.
    Note the "bin" column is ordered from 0 to X."""

    y = condition_check(df, el)
    idx = y["dist"] == 1
    while idx.any():
        __left = el[["cluster1", "cluster2", "weight"]]
        __left.columns = ["path_node1", "path_node2", "old_path"]
        __left.set_index(["path_node1", "path_node2"], inplace=True)
        path_node2 = y.loc[idx.shift().fillna(False), "cluster2"].values
        weight2 = y.loc[idx.shift().fillna(False), "weight"].values

        values = y.loc[idx, ["cluster1", "cluster2", "weight", "bin"]]

        __right = pd.DataFrame().assign(cluster=values["cluster1"].values)
        __right = __right.assign(path_node1=values["cluster2"].values)
        __right = __right.assign(path_node2=path_node2)
        __right = __right.assign(weight=values["weight"].values)
        __right = __right.assign(weight2=weight2)
        __right = __right.assign(bin=values["bin"].values)
        __right.set_index(["path_node1", "path_node2"], inplace=True)
        merged = pd.merge(__left, __right, left_index=True, right_index=True)
        # Now select the ONE node with the least distance
        result = merged.eval("diff = weight + weight2 - old_path").sort_values("diff").head(1)
        bin_idx = result["bin"].values[0]
        df = pd.concat([
            df.iloc[:bin_idx],
            pd.DataFrame({"cluster": result["cluster"], "bin": bin_idx + 1}),
            df.iloc[bin_idx:].eval("bin = bin + 1")
        ])
        y = condition_check(df, el)
        idx = y["dist"] == 1

    return df
