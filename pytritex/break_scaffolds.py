import pandas as pd
import numpy as np
from .anchor_scaffolds import anchor_scaffolds
from .add_hic_cov import add_hic_cov
from .add_molecule_cov import add_molecule_cov
from .utils.rolling_join import rolling_join


def _transpose_cssaln(cssaln, fai):
    z = cssaln.copy()
    del z["scaffold_length"]
    del z["scaffold"]
    f = fai.loc[:, ["scaffold", "length", "orig_scaffold", "orig_start"]].rename(
        columns={"length": "scaffold_length"}).assign(orig_pos=lambda df: df["orig_start"])
    assert "scaffold_length" in f.columns, f.columns
    orig_columns = z.columns
    z = f.merge(z, on=["orig_scaffold", "orig_pos"], how="right", suffixes=("", "_z"))
    assert "scaffold_length" in z.columns, z.columns
    missing_bait = z["scaffold_length"].isna()
    z.loc[missing_bait, ["scaffold", "scaffold_length"]] = f.merge(z.loc[missing_bait, orig_columns],
                                                                   on=["orig_scaffold"], how="right",
                                                                   suffixes=("", "_z")
                                                                   )[["scaffold", "scaffold_length"]]
    z.loc[:, "pos"] = z["orig_pos"] - z["orig_start"] + 1
    del z["orig_start"]
    return z


def _transpose_10x_cov(broken_assembly: dict, info, breaks, cores, prefix, regex2) -> dict:
    if "molecule_cov" in assembly and assembly["molecule_cov"].shape[0] > 0:
        print("10X molecule coverage")
        cov = add_molecule_cov(
            broken_assembly, scaffolds=broken_assembly["info"]["scaffold"],
            binsize=broken_assembly["mol_binsize"], cores=cores)
        # #   info[!breaks$scaffold, on="scaffold"]->x
        x = info.loc[~info["scaffold"].isin(breaks["scaffold"]), :].copy()
        x.loc[:, "scaffold"] = prefix + x.loc[:, "scaffold"].str.replace(regex2, "\\2", regex=True)
        x.loc[:, "split"] = True
        print([column for column in cov["info"] if column not in x.columns])
        broken_assembly["info"] = pd.concat([x.loc[:, cov["info"].columns], cov["info"]])
        broken_assembly["mol_binsize"] = assembly["mol_binsize"]
        x = assembly["molecule_cov"].loc[~assembly["molecule_cov"]["scaffold"].isin(breaks["scaffold"]), :].copy()
        x.loc[:, "scaffold"] = prefix + x.loc[:, "scaffold"].str.replace(regex2, "\\2", regex=True)
        if cov["molecule_cov"].shape[0] > 0:
            broken_assembly["molecule_cov"] = pd.concat([x, cov["molecule_cov"]])
        else:
            broken_assembly["molecule_cov"] = x
        print("Molecule cov:", broken_assembly["molecule_cov"].shape[0])
    elif "molecule_cov" in assembly and broken_assembly["molecule_cov"].shape[0] == 0:
        raise ValueError()
    else:
        broken_assembly["molecule_cov"] = pd.DataFrame()
    return broken_assembly


def _transpose_hic_cov(assembly_new, assembly, fai, info, breaks, cores, prefix, regex2):
    if "cov" in assembly and assembly["fpairs"].shape[0] > 0:
        print("Hi-C coverage")
        assert assembly_new["fpairs"]["pos1"].min() >= 0
        cov = add_hic_cov(assembly=assembly_new,
                          scaffolds=fai.loc[fai["split"] == True, "scaffold"],
                          binsize=assembly["binsize"],
                          minNbin=assembly["minNbin"],
                          innerDist=assembly["innerDist"],
                          cores=cores)
        x = assembly["cov"].loc[~assembly["cov"]["scaffold"].isin(breaks["scaffold"]), :].copy()
        x.loc[:, "scaffold"] = prefix + x.loc[:, "scaffold"].str.replace(regex2, "\\2", regex=True)
        if cov["cov"].shape[0] > 0:
            assembly_new["cov"] = pd.concat([x, cov["cov"]])
        else:
            assembly_new["cov"] = x
        x = info.loc[~info["scaffold"].isin(breaks["scaffold"]), :].copy()
        x.loc[:, "scaffold"] = prefix + x.loc[:, "scaffold"].str.replace(regex2, "\\2", regex=True)
        x.loc[:, "split"] = False
        assembly_new["info"] = pd.concat([x.loc[:, cov["info"].columns], cov["info"]])
    else:
        assembly_new["cov"] = pd.DataFrame()
    return assembly_new


def _transpose_fpairs(fpairs: pd.DataFrame, fai: pd.DataFrame) -> pd.DataFrame:
    if fpairs is not None and fpairs.shape[0] > 0:
        assert fpairs["pos1"].min() > 0
        print("Starting to transpose Fpairs")
        z = fpairs.loc[:, ["orig_scaffold1", "orig_scaffold2", "orig_pos1", "orig_pos2"]]
        z.loc[:, "orig_pos1"] = z.loc[:, "orig_pos1"].astype(int)
        z.loc[:, "orig_pos2"] = z.loc[:, "orig_pos2"].astype(int)
        f1 = fai.rename(
            columns={"scaffold": "scaffold1",
                     "orig_scaffold": "orig_scaffold1",
                     "orig_start": "orig_start1"}).assign(
            orig_pos1=fai["orig_start"]).loc[:, ["scaffold1", "orig_scaffold1", "orig_start1", "orig_pos1"]]
        f1.loc[:, "orig_pos1"] = f1.loc[:, "orig_pos1"].astype(int)
        f2 = fai.rename(
            columns={"scaffold": "scaffold2",
                     "orig_scaffold": "orig_scaffold2",
                     "orig_start": "orig_start2"}).assign(
            orig_pos2=fai["orig_start"]).loc[:, ["scaffold2", "orig_scaffold2", "orig_start2", "orig_pos2"]]
        f2.loc[:, "orig_pos2"] = f2.loc[:, "orig_pos2"].astype(int)
        z = rolling_join(f1, z, ["orig_scaffold1"], "orig_pos1", ["orig_scaffold1", "orig_scaffold2"])
        z = rolling_join(f2, z, ["orig_scaffold2"], "orig_pos2", ["orig_scaffold1", "orig_scaffold2"])
        remainder = fpairs.loc[~((fpairs["orig_scaffold1"].isin(z["orig_scaffold1"])) |
                                 (fpairs["orig_scaffold1"].isin(z["orig_scaffold2"])))]
        assert remainder.shape[0] == 0
        z.loc[:, "pos1"] = z["orig_pos1"] - z["orig_start1"] + 1
        assert z["pos1"].min() > 0
        z.loc[:, "pos2"] = z["orig_pos2"] - z["orig_start2"] + 1
        assert z["pos2"].min() > 0
        del z["orig_start1"]
        del z["orig_start2"]
        return z
    else:
        return pd.DataFrame()


def _transpose_molecules(molecules: pd.DataFrame, fai: pd.DataFrame) -> pd.DataFrame:
    if molecules is not None and molecules.shape[0] > 0:
        print("Transpose molecules")
        z = molecules.copy()
        del z["scaffold"]
        #   fai[, .(scaffold, orig_scaffold, orig_start,
        #   s_length=orig_end - orig_start + 1, orig_pos=orig_start)][z, on=c("orig_scaffold", "orig_start"), roll=T]->z
        z = rolling_join(
            fai.loc[:, ["scaffold", "orig_scaffold", "orig_start"]].assign(
                s_length=fai["orig_end"] - fai["orig_start"] + 1,
                orig_pos=fai["orig_start"]
            ),
            z, on=["orig_scaffold"], by="orig_start", right_key=["orig_scaffold"]
        )
        z = z.assign(start=z["orig_start"] - z["orig_pos"] + 1,
                     end=z["orig_end"] - z["orig_pos"] + 1)
        z = z.loc[z.end <= z.s_length, :].copy()
        del z["orig_pos"]
        del z["s_length"]
        print(z.head())
        print(z.shape[0])
        print(z.dropna().shape[0])
        return z.dropna()
    else:
        return pd.DataFrame()


def _concatenate_br_and_fai(br, fai, prefix, regex2):
    assert br[br["length1"].isna()].shape[0] == 0, br[br["length1"].isna()]
    first = br.assign(scaffold=br["scaffold1"],
                      length=br["length1"],
                      orig_start=br["orig_start"] + br["start1"] - 1,
                      orig_end=br["orig_start"] + br["end1"] - 1)[
         ["scaffold", "length", "orig_scaffold", "orig_start", "orig_end", "old_scaffold"]]
    assert first[first["length"].isna()].shape[0] == 0, first[first["length"].isna()]

    second = br.assign(scaffold=br["scaffold2"],
                   length=br["length2"],
                   orig_start=br["orig_start"] + br["start2"] - 1,
                   orig_end=br["orig_start"] + br["end2"] - 1)[
         ["scaffold", "length", "orig_scaffold", "orig_start", "orig_end", "old_scaffold"]]
    assert second[second["length"].isna()].shape[0] == 0, second[second["length"].isna()]

    third = br.assign(scaffold=br["scaffold3"],
                   length=br["length3"],
                   orig_start=br["orig_start"] + br["start3"] - 1,
                   orig_end=br["orig_start"] + br["end3"] - 1)[
         ["scaffold", "length", "orig_scaffold", "orig_start", "orig_end", "old_scaffold"]]
    assert third[third["length"].isna()].shape[0] == 0, third[third["length"].isna()]

    remainder = fai.loc[~fai["scaffold"].isin(br["scaffold"]),
                 ["scaffold", "orig_scaffold", "length", "orig_start", "orig_end", "old_scaffold"]]
    remainder = remainder.assign(
             scaffold=prefix + remainder["scaffold"].str.replace(regex2, "\\2", regex=True)
         )
    assert remainder[remainder["length"].isna()].shape[0] == 0, remainder[remainder["length"].isna()]

    fai = pd.concat([first, second, third, remainder]).reset_index()
    return fai


def _calculate_breaks(br, fai, prefix, regex1, slop, regex2):
    assert br[br["br"].isna()].shape[0] == 0, br
    br = fai.merge(br, how="right", on="scaffold")
    assert br[br["length"] == 0].shape[0] == 0
    br.loc[:, "orig_br"] = br["orig_start"] + br["br"] - 1
    br = br.sort_values(["scaffold", "br"])
    nbr = br.loc[br["scaffold"].duplicated(), :]
    br = br.loc[~br["scaffold"].duplicated(), :]
    maxidx = fai["scaffold"].str.replace(regex1, "\\3", regex=True).astype(int).max()
    br.loc[:, "idx"] = 3 * np.arange(1, br.shape[0] + 1) - 2
    br.loc[:, "scaffold1"] = prefix + pd.Series(maxidx + br["idx"]).astype(str)
    br.loc[:, "start1"] = 1
    assert br[br["br"].isna()].shape[0] == 0, br[br["br"].isna()]
    br.loc[:, "end1"] = np.maximum(0, br["br"] - slop - 1)
    # br[, scaffold2 := paste0(prefix, maxidx + idx + 1)]
    br.loc[:, "scaffold2"] = prefix + pd.Series(maxidx + br["idx"] + 1).astype(str)
    # br[, start2 := pmax(1, br - slop)]
    br.loc[:, "start2"] = np.maximum(1, br["br"] - slop)
    # br[, end2 := pmin(br + slop - 1, length)]
    br.loc[:, "end2"] = np.minimum(br["br"] + slop - 1, br["length"])
    #   br[, scaffold3 := paste0(prefix, maxidx+idx+2)]
    br.loc[:, "scaffold3"] = prefix + pd.Series(maxidx + br["idx"] + 2).astype(str)
    #   br[, start3 := pmin(length + 1, br + slop)]
    br.loc[:, "start3"] = np.minimum(br["length"] + 1, br["br"] + slop)
    #   br[, end3 := length]
    br.loc[:, "end3"] = br.loc[:, "length"]
    #   br[, length1 := 1 + end1 - start1]
    br.loc[:, "length1"] = 1 + br["end1"] - br["start1"]
    #   br[, length2 := 1 + end2 - start2]
    br.loc[:, "length2"] = 1 + br["end2"] - br["start2"]
    #   br[, length3 := 1 + end3 - start3]
    br.loc[:, "length3"] = 1 + br["end3"] - br["start3"]
    fai = _concatenate_br_and_fai(br, fai, prefix, regex2)
    assert fai[fai["length"].isna()].shape[0] == 0, fai[fai["length"].isna()]
    # assert fai[fai["length"] == 0].shape[0] == 0, fai[fai["length"] == 0]
    fai = fai.loc[fai["length"] > 0, :]
    # The "roll" functionality is NOT present in pandas. It is necessary to do a double merge.
    current_cols = ["orig_scaffold", "orig_br"]
    f = fai.loc[:, ["scaffold", "orig_scaffold", "orig_start"]].assign(orig_br=fai["orig_start"])
    nbr = f.merge(nbr.loc[:, current_cols], on=["orig_scaffold", "orig_br"], how="right")
    nbr.loc[nbr["scaffold"].isna(), ["scaffold", "orig_start"]] = f.merge(
        nbr.loc[nbr["scaffold"].isna(), current_cols],
        on=["orig_scaffold"]).loc[:, ["scaffold", "orig_start"]]
    nbr.loc[:, "br"] = nbr["orig_br"] - nbr["orig_start"] + 1
    br = nbr.loc[:, ["scaffold", "br"]]
    return br, nbr, fai


# # Break scaffolds at specified points and lift positional information to updated assembly
def break_scaffolds(breaks: pd.DataFrame, assembly: dict, prefix: str,
                    slop: float, cores=1, species="wheat",
                    regex1=r"(^.*[^-0-9])(([0-9]+)(-[0-9]+)?$)",
                    regex2=r"(^.*[^-0-9])([0-9]+(-[0-9]+)?$)"):

    info, cov, fpairs, cssaln = assembly["info"], assembly["cov"], assembly["fpairs"], assembly["cssaln"]
    # Do we have the HiC coverage? Check
    fpairs_present = (assembly.get("fpairs", None) is not None) and (assembly["fpairs"].shape[0] > 0)
    # Copy the breaks table
    br = breaks.copy()
    del br["length"]
    fai = info.loc[:, ["scaffold", "orig_scaffold", "orig_start", "orig_end", "length"]]
    fai.loc[:, "old_scaffold"] = fai.loc[:, "scaffold"]
    fai.loc[:, "split"] = fai["scaffold"].isin(br["scaffold"])
    stable = fai.loc[~fai["scaffold"].isin(br["scaffold"])]
    fai_to_break = fai.loc[fai["scaffold"].isin(br["scaffold"])]
    print("Split scaffolds")
    j = 0
    while br.shape[0] > 0:
        j += 1
        o = fai_to_break.shape[0]
        br, nbr, fai_to_break = _calculate_breaks(br, fai_to_break, prefix, regex1=regex1, regex2=regex2, slop=slop)
        assert o <= fai_to_break.shape[0], (o, fai_to_break.shape[0])
        print("Iteration", j, "finished.")
        if o < fai_to_break.shape[0]:
            print("The number of scaffolds increased from ", o, " to {}.".format(fai_to_break.shape[0]))

    old_shape = fai.shape[0]
    stable_assembly = {"binsize": assembly["binsize"], "innerDist": assembly["innerDist"],
                       "minNbin": assembly["minNbin"]}
    broken_assembly = {"binsize": assembly["binsize"], "innerDist": assembly["innerDist"],
                       "minNbin": assembly["minNbin"]}
    broken_info = info.loc[info["scaffold"].isin(breaks["scaffold"])]
    stable_assembly["info"] = stable
    broken_assembly["info"] = fai_to_break
    print("Transpose cssaln")
    css_to_break = cssaln["scaffold"].isin(breaks["scaffold"])
    stable_assembly["cssaln"] = cssaln.loc[~css_to_break]
    broken_assembly["cssaln"] = _transpose_cssaln(cssaln.loc[css_to_break], fai_to_break)
    if fpairs_present is True:
        fpairs_to_break = (fpairs["scaffold1"].isin(breaks["scaffold"]) |
                           fpairs["scaffold2"].isin(breaks["scaffold"]))
        stable_assembly["fpairs"] = fpairs.loc[~fpairs_to_break]
        broken_assembly["fpairs"] = _transpose_fpairs(fpairs.loc[fpairs_to_break], broken_assembly["info"])
        assert breaks.shape[0] == 0 or (broken_assembly["fpairs"] is not None and
                                        broken_assembly["fpairs"].shape[0] > 0)
        assert breaks.shape[0] == 0 or broken_assembly["fpairs"]["pos1"].min() > 0
        assert breaks.shape[0] == 0 or broken_assembly["fpairs"]["pos2"].min() > 0
    molecules = assembly.get("molecules", None)
    if molecules is not None:
        molecules_to_break = molecules["scaffold"].isin(breaks["scaffold"])
        stable_assembly["molecules"] = molecules.loc[~molecules_to_break]
        broken_assembly["molecules"] = _transpose_molecules(molecules.loc[molecules_to_break],
                                                            broken_assembly["info"])
    else:
        stable_assembly["molecules"] = None
        broken_assembly["molecules"] = None
    print("Anchor scaffolds")
    broken_assembly = anchor_scaffolds(broken_assembly, popseq=assembly["popseq"], species=species)
    broken_assembly["info"].drop("mr_10x", axis=1, errors="ignore")
    broken_assembly["info"].drop("mr", axis=1, errors="ignore")
    broken_assembly["info"].drop("mri", axis=1, errors="ignore")
    broken_assembly = _transpose_10x_cov(broken_assembly, broken_info,
                                         breaks=breaks, cores=cores, prefix=prefix, regex2=regex2)
    broken_assembly = _transpose_hic_cov(broken_assembly, broken_info,
                                         breaks=breaks, cores=cores, prefix=prefix, regex2=regex2)

    assembly_new = dict()
    for key in ["info", "molecules", "cssaln", ]


    print("Finished anchoring after breaking chimeras")
    return assembly_new
