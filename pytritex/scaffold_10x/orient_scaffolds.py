import pandas as pd
import numpy as np
from .make_agp import make_agp


def orient_scaffolds(info: pd.DataFrame, res: pd.DataFrame,
                     membership: pd.DataFrame,
                     link_pos: pd.DataFrame,
                     max_dist_orientation: float, verbose=False):

    if verbose:
        print("Orienting scaffolds")

    # orient scaffolds
    a = membership.loc[
        membership["super_nbin"] > 1, ["scaffold_index", "bin", "super"]].rename(
        columns={"scaffold_index": "scaffold_index1", "bin": "bin1", "super": "super1"}).merge(
        link_pos, on=["scaffold_index1"], how="left")
    a = membership.loc[
        membership["super_nbin"] > 1, ["scaffold_index", "bin", "super"]].rename(
        columns={"scaffold_index": "scaffold_index2", "bin": "bin2", "super": "super2"}).merge(
        a, on=["scaffold2"], how="left")
    a = a.query("(super1 == super2) & (bin1 != bin2)").eval("d = abs(bin2 - bin1)").query(
        "d <= @max_dist_orientation")

    def _previous(series):
        return (series.loc[a["bin2"] < a["bin1"]]).mean()

    def _next(series):
        return (series.loc[a["bin1"] < a["bin2"]]).mean()

    aa = a.rename(columns={"scaffold_index1": "scaffold_index"}).groupby("scaffold_index").agg(
        prv=("pos1", _previous),
        nxt=("pos1", _next))
    aa = info.loc[:, ["scaffold_index", "length"]].set_index("scaffold_index").merge(
        aa, left_index=True, right_index=True, how="right").reset_index(drop=False)

    bait = (~aa["prv"].isna()) & (~aa["nxt"].isna())
    aa.loc[bait, "orientation"] = (aa.loc[bait, "nxt"] >= aa.loc[bait, "prv"]).map({True: 1, False: -1})
    bait = (aa["prv"].isna()) & (~aa["nxt"].isna())
    aa.loc[bait, "orientation"] = (aa.loc[bait, "nxt"] >= aa.loc[bait, "length"] - aa.loc[bait, "nxt"]).map(
        {True: 1, False: -1})
    bait = (~aa["prv"].isna()) & (aa["nxt"].isna())
    aa.loc[bait, "orientation"] = (aa.loc[bait, "prv"] >= aa.loc[bait, "length"] - aa.loc[bait, "prv"]).map(
        {True: 1, False: -1})
    membership = aa.loc[:, ["orientation", "scaffold"]].merge(membership, on=["scaffold"], how="right")
    membership.loc[:, "oriented"] = True
    membership.loc[membership["orientation"].isna(), "oriented"] = False
    membership.loc[membership["orientation"].isna(), "orientation"] = 1
    membership = membership.sort_values(["super", "bin", "rank"])
    # For each group (defined by super) define the super_pos as the previous cumulative length + 1
    # Basically, where we are at in the super scaffold
    membership.loc[:, "super_pos"] = membership.groupby("super")["length"].shift().transform(np.cumsum).fillna(0) + 1
    if verbose is True:
        print("Anchoring super-scaffolds")

    # assign super scaffolds to genetic positions
    __left = membership.merge(
        membership[(~membership["chr"].isna())].groupby(["chr", "super"])["chr"].transform("count").to_frame("nchr"),
    how="right", left_index=True, right_index=True)
    y = __left.merge((__left["nchr"] / __left.groupby(["super"])["nchr"].transform("sum")).to_frame("pchr"), left_index=True,
                 right_index=True).sort_values("nchr", ascending=False).loc[lambda df: ~df["super"].duplicated()]
    yy = membership.loc[~membership["cM"].isna()].merge(y, on=["chr", "super"], how="right")
    res = y.merge(res, on="super", how="right")

    res = yy.groupby("super").agg(cM=("cM", "mean"), min_cM=("cM", "min"),
                                  max_cM=("cM", "max")).merge(res, how="right", left_index=True, right_on="super")
    res = res.sort_values("length", ascending=False)
    agp = make_agp(membership, gap_size=100)

    result = {"membership": membership, "info": res, "agp": agp["agp"], "agp_bed": agp["agp_bed"]}

    return result
