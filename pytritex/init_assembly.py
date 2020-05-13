import pandas as pd


def init_assembly(fai: pd.DataFrame, cssaln: pd.DataFrame, fpairs=None, molecules=None, rename=None) -> dict:
    info = fai.copy()
    info.loc[:, "orig_start"] = 1
    info.loc[:, "orig_end"] = info["length"]
    info.loc[:, "orig_scaffold"] = info["scaffold"]
    if rename is not None:
        assert isinstance(rename, dict)
        info.loc[:, "scaffold"] = info["scaffold"].map(rename)
    z = cssaln.copy()
    z.loc[:, "orig_scaffold"] = z["scaffold"]
    if rename is not None:
        z.loc[:, "scaffold"] = z["scaffold"].map(rename)
    z.loc[:, "orig_pos"] = z["pos"]
    z.loc[:, "orig_scaffold_length"] = z["scaffold_length"]
    if molecules is not None:
        y = molecules.copy()
        y.rename(columns={"scaffold": "orig_scaffold"}, inplace=True)
        y.loc[:, "orig_start"] = y["start"]
        y.loc[:, "orig_end"] = y["end"]
        # TODO: This is absolutely NOT clear. Why should have the scaffold changed in the info?
        # TODO: Why not simpy a column copy?
        y = info[["orig_scaffold", "scaffold"]].merge(y, on="orig_scaffold", how="right")
    else:
        y = pd.DataFrame()

    if fpairs is not None:
        tcc = fpairs.copy()
        tcc.rename(columns={"scaffold1": "orig_scaffold1", "scaffold2": "orig_scaffold2"},
                   inplace=True)
        tcc.loc[:, "orig_pos1"] = tcc["pos1"]
        tcc.loc[:, "orig_pos2"] = tcc["pos2"]
        tcc = pd.merge(info[["orig_scaffold", "scaffold"]].rename(columns={
            "orig_scaffold": "orig_scaffold1", "scaffold": "scaffold1"
        }), tcc, left_on="orig_scaffold1", right_on="orig_scaffold1", how="right")
        tcc = pd.merge(info[["orig_scaffold", "scaffold"]].rename(columns={
            "orig_scaffold": "orig_scaffold2", "scaffold": "scaffold2"
        }), tcc, left_on="orig_scaffold2", right_on="orig_scaffold2", how="right")
    else:
        tcc = pd.DataFrame()
    res = {"info": info, "cssaln": z, "fpairs": tcc, "molecules": y}
    return res
