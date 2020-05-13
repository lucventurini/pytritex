import graph_tool as gt
import graph_tool.topography
import pandas as pd
import numpy as np


def _df_checker(df: pd.DataFrame, el: pd.DataFrame, cluster: pd.Series):
    m = pd.DataFrame().assign(
        cluster1=df.loc[range(df.shape[0] - 1)].reset_index()["cluster"],
        cluster2=df.loc[range(1, df.shape[0])].reset_index()["cluster"])
    # Merge with el based on index "cluster1", "cluster2"; rename weight to weight12; assign a dummy variable to 1
    m = pd.merge(m.set_index(["cluster1", "cluster2"]),
                 el.reset_index().set_index(["cluster1", "cluster2"]), left_index=True, right_index=True
                 )["weight"].reset_index().rename(columns={"weight": "weight12"}).assign(dummy=1)
    m = pd.merge(
        df[:].rename(columns={"cluster": "cluster1", "bin": "bin1"}).set_index("cluster1"),
        m.reset_index().set_index("cluster1"),
        left_index=True, right_index=True).reset_index()
    m = pd.merge(
        df[:].rename(columns={"cluster": "cluster2", "bin": "bin2"}).set_index("cluster2"),
        m.set_index("cluster2"),
        left_index=True, right_index=True).reset_index()
    n = m[:].rename(columns=dict(
        zip(
            ("cluster1", "cluster2", "bin1", "bin2", "weight12"),
            ("cluster3", "cluster4", "bin3", "bin4", "weight34")
        ))
    )
    # Cartesian product
    mn = pd.merge(m, n, how="outer", left_on="dummy", right_on="dummy").loc[
        eval("lambda df: df['bin1'] < df['bin3']")].assign(dummy=np.nan)
    o = el[["cluster1", "cluster2", "weight"]]
    mn = pd.merge(
        o[:].rename(columns={"cluster2": "cluster3", "weight": "weight13"}), mn,
        left_on=["cluster1", "cluster3"], right_on=["cluster1", "cluster3"])
    mn = pd.merge(
        o[:].rename(columns={"cluster1": "cluster2", "cluster2": "cluster4", "weight": "weight24"}), mn,
        left_on=["cluster2", "cluster4"], right_on=["cluster2", "cluster4"])
    mn = mn.assign(old=mn["weight12"] + mn["weight34"], new=mn["weight13"] + mn["weight24"])
    mn = mn.assign(diff=mn["old"] - mn["new"]).sort_values("diff", ascending=False)
    return mn, ((mn["diff"] > 0).shape[0] > 0)


def kopt2(df, el, cluster):
    # Traveling salesman heuristic for Hi-C mapping construction
    el = el.set_index(["cluster1", "cluster2"])
    mn, execute = _df_checker(df, el, cluster)
    while execute is True:
        x = mn.head(1)
        bin1 = df[df["cluster"] == x["cluster1"]]["bin"]
        bin2 = df[df["cluster"] == x["cluster2"]]["bin"]
        bin3 = df[df["cluster"] == x["cluster3"]]["bin"]
        bin4 = df[df["cluster"] == x["cluster4"]]["bin"]
        bin_collection = list(range(1, bin1 + 1)) + list(range(bin3, bin2 + 1)) + list(range(bin4, df.shape[0] + 1))
        df = df.loc[bin_collection]
    return df
