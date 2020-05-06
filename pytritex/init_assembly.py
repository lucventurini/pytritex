import pandas as pd


def init_assembly(fai: pd.DataFrame, cssaln: pd.DataFrame, fpairs=None, molecules=None, rename=None):
    info = fai[:]
    info["orig_start"] = 1
    info["orig_end"] = info["length"]
    info["orig_scaffold"] = info["scaffold"]
    if rename is not None:
        assert isinstance(rename, dict)
        info["scaffold"] = info["scaffold"].map(rename)
    z = cssaln[:]
    z["orig_scaffold"] = z["scaffold"]
    if rename is not None:
        z["scaffold"] = z["scaffold"].map(rename)
    z["orig_pos"] = z["pos"]
    z["orig_scaffold_length"] = z["scaffold_length"]
    if molecules is not None:
        y = molecules[:]
        y["orig_scaffold"] = y["scaffold"]
        y["orig_start"] = y["start"]
        y["orig_end"] = y["end"]
        y = pd.merge(info[["orig_scaffold", "scaffold"]], y, left_on="orig_scaffold", right_on="orig_scaffold")
    else:
        y = pd.DataFrame()

    if fpairs is not None:
        tcc = fpairs[:]
        tcc["orig_scaffold1"] = tcc["scaffold1"]
        tcc["orig_pos1"] = tcc["pos1"]
        tcc["orig_scaffold2"] = tcc["scaffold2"]
        tcc["orig_pos2"] = tcc["pos2"]
        tcc = pd.merge(info[["orig_scaffold", "scaffold"]].rename(columns={
            "orig_scaffold": "orig_scaffold1", "scaffold": "scaffold1"
        }), tcc, left_on="orig_scaffold1", right_on="orig_scaffold1")
        tcc = pd.merge(info[["orig_scaffold", "scaffold"]].rename(columns={
            "orig_scaffold": "orig_scaffold2", "scaffold": "scaffold2"
        }), tcc, left_on="orig_scaffold2", right_on="orig_scaffold2")
    else:
        tcc = pd.DataFrame()
    res = {"info": info, "cssaln": z, "fpairs": tcc, "molecules": y}
    return res
