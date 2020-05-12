import pandas as pd
import numpy as np
from .anchor_scaffolds import anchor_scaffolds
from .add_hic_cov import add_hic_cov
from .add_molecule_cov import add_molecule_cov


def _transpose_10x_cov(assembly: dict, assembly_new: dict, fai, info, breaks, cores, prefix, regex2) -> dict:
    if "molecule_cov" in assembly and assembly["molecule_cov"].shape[0] > 0:
        print("10X molecule coverage")
        cov = add_molecule_cov(
            assembly_new, scaffolds=fai.loc[fai["split"] == True, ["scaffold"]],
            binsize=assembly["mol_binsize"], cores=cores)
        # #   info[!breaks$scaffold, on="scaffold"]->x
        x = info.loc[~info["scaffold"].isin(breaks["scaffold"]), :]
        x.loc[:, "scaffold"] = prefix + x["scaffold"].str.replace(regex2, "\\2", regex=True)
        x.loc[:, "split"] = True
        assembly_new["info"] = pd.concat([x.loc[:, cov["info"].columns], cov["info"]])
        assembly_new["mol_binsize"] = assembly["mol_binsize"]
        x = assembly["molecule_cov"].loc[~assembly["molecule_cov"]["scaffold"].isin(breaks["scaffold"]), :]
        x.loc[:, "scaffold"] = prefix + x["scaffold"].str.replace(regex2, "\\2", regex=True)
        if cov["molecule_cov"].shape[0] > 0:
            assembly_new["molecule_cov"] = pd.concat([x, cov["molecule_cov"]])
        else:
            assembly_new["molecule_cov"] = x
    else:
        assembly_new["molecule_cov"] = pd.DataFrame()
    return assembly_new


def _transpose_hic_cov(assembly_new, assembly, fai, info, breaks, cores, prefix, regex2):
    if "cov" in assembly and assembly["fpairs"].shape[0] > 0:
        print("Hi-C coverage")
        cov = add_hic_cov(assembly=assembly_new,
                          scaffolds=fai.loc[fai["split"] == True, "scaffold"],
                          binsize=assembly["binsize"],
                          minNbin=assembly["minNbin"],
                          innerDist=assembly["innerDist"],
                          cores=cores)
        x = assembly["cov"].loc[~assembly["cov"]["scaffold"].isin(breaks["scaffold"]), :]
        x.loc[:, "scaffold"] = prefix + x["scaffold"].str.replace(regex2, "\\2", regex=True)
        if cov["cov"].shape[0] > 0:
            assembly_new["cov"] = pd.concat(x, cov["cov"])
        else:
            assembly_new["cov"] = x
        x = info.loc[~info["scaffold"].isin(breaks["scaffold"]), :]
        x.loc[:, "scaffold"] = prefix + x["scaffold"].str.replace(regex2, "\\2", regex=True)
        x.loc[:, "split"] = False
        assembly_new["info"] = pd.concat([x.loc[:, cov["info"].columns], cov["info"]])
    else:
        assembly_new["cov"] = pd.DataFrame()
    return assembly_new


def _transpose_fpairs(fpairs: pd.DataFrame, fai: pd.DataFrame) -> pd.DataFrame:
    if fpairs and fpairs.shape[0] > 0:
        z = fpairs.loc[:, ["orig_scaffold1", "orig_scaffold2", "orig_pos1", "orig_pos2"]]
        # TODO: what's the roll?
        #   fai[, .(scaffold1=scaffold, orig_scaffold1=orig_scaffold, orig_start1=orig_start, orig_pos1=orig_start)
        #   ][z, on=c("orig_scaffold1", "orig_pos1"), roll=T]->z
        z = fai.loc[:, ["scaffold", "orig_scaffold", "orig_start"]].assign(orig_pos1=fai["orig_start"]).rename(
            columns={"scaffold": "scaffold1", "orig_scaffold": "orig_scaffold1", "orig_start1": "orig_start"}
        ).merge(z, on=["orig_scaffold1", "orig_pos1"])
        z = fai.loc[:, ["scaffold", "orig_scaffold", "orig_start"]].assign(orig_pos2=fai["orig_start"]).rename(
            columns={"scaffold": "scaffold2", "orig_scaffold": "orig_scaffold2", "orig_start2": "orig_start"}
        ).merge(z, on=["orig_scaffold2", "orig_pos2"])
        z.loc[:, "pos1"] = z["orig_pos1"] - z["orig_start1"] + 1
        z.loc[:, "pos2"] = z["orig_pos2"] - z["orig_start2"] + 1
        z.loc[:, "orig_start1"] = np.nan
        z.loc[:, "orig_start2"] = np.nan
        return z
    else:
        return pd.DataFrame()


def _transpose_molecules(molecules: pd.DataFrame, fai: pd.DataFrame) -> pd.DataFrame:
    if molecules and molecules.shape[0] > 0:
        print("Transpose molecules")
        z = molecules.copy()
        z.loc[:, "scaffold"] = np.nan
        #   fai[, .(scaffold, orig_scaffold, orig_start, s_length=orig_end - orig_start + 1, orig_pos=orig_start)][z, on=c("orig_scaffold", "orig_start"), roll=T]->z

        z.loc[:, "start"] = z.eval("orig_start - orig_pos + 1")
        z.loc[:, "end"] = z.eval("orig_end - orig_pos + 1")
        z = z.loc[z.eval("end <= s_length"), :]  # TODO: is s_length in the table?
        z.loc[:, "orig_pos"] = np.nan
        z.loc[:, "s_length"] = np.nan
        return z
    else:
        return pd.DataFrame()


def _concatenate_br_and_fai(br, fai, prefix, regex2):
    fai = pd.concat(
        #    br[, .(scaffold=scaffold1, length=length1, orig_scaffold, orig_start=orig_start+start1-1, orig_end=orig_start+end1-1, old_scaffold)],
        [br.assign(scaffold=br["scaffold1"],
                   length=br["length1"],
                   orig_start=br["orig_start"] + br["start1"] - 1,
                   orig_end=br["orig_start"] + br["end1"] - 1).loc[:,
         ["scaffold", "length", "orig_scaffold", "orig_start", "orig_end", "old_scaffold"]],
         #    br[, .(scaffold=scaffold2, length=length2, orig_scaffold, orig_start=orig_start+start2-1, orig_end=orig_start+end2-1, old_scaffold)],
         br.assign(scaffold=br["scaffold2"],
                   length=br["length2"],
                   orig_start=br["orig_start"] + br["start2"] - 1,
                   orig_end=br["orig_start"] + br["end2"] - 1).loc[:,
         ["scaffold", "length", "orig_scaffold", "orig_start", "orig_end", "old_scaffold"]],
         #    br[, .(scaffold=scaffold3, length=length3, orig_scaffold, orig_start=orig_start+start3-1, orig_end=orig_start+end3-1, old_scaffold)],
         br.assign(scaffold=br["scaffold3"],
                   length=br["length3"],
                   orig_start=br["orig_start"] + br["start3"] - 1,
                   orig_end=br["orig_start"] + br["end3"] - 1).loc[:,
         ["scaffold", "length", "orig_scaffold", "orig_start", "orig_end", "old_scaffold"]],
         #    fai[!scaffold %in% br$scaffold, .(orig_scaffold, length,  orig_start, orig_end, old_scaffold,
         # 			      scaffold=paste0(prefix, sub(regex2, "\\2", scaffold)))]
         fai.loc[~fai["scaffold"].isin(br["scaffold"]),
                 ["orig_scaffold", "length", "orig_start", "orig_end", "old_scaffold"]].assign(
             scaffold=prefix + fai["scaffold"].str.replace(regex2, "\\2", regex=True)
         )]
    )
    return fai


def _calculate_breaks(br, fai, prefix, regex1, slop, regex2):
    br = fai.merge(br, how="right", on="scaffold")
    br.loc[:, "orig_br"] = br["orig_start"] + br["br"] - 1
    br = br.sort_values(["scaffold", "br"])
    nbr = br.loc[br["scaffold"].duplicated(), :]
    br = br.loc[~br["scaffold"].duplicated(), :]
    maxidx = fai["scaffold"].str.replace(regex1, "\\3", regex=True).astype(int).max()
    br.loc[:, "idx"] = 3 * np.arange(1, br.shape[0] + 1) - 2
    # br[, scaffold1 := paste0(prefix, maxidx + idx)]
    br.loc[:, "scaffold1"] = prefix + pd.Series(maxidx + br["idx"]).astype(str)
    # br[, start1 := 1]
    br.loc[:, "start1"] = 1
    # br[, end1 := pmax(0, br - slop - 1)]
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
    br.loc[:, "end3"] = "length"
    #   br[, length1 := 1 + end1 - start1]
    br.loc[:, "length1"] = 1 + br["end1"] - br["start1"]
    #   br[, length2 := 1 + end2 - start2]
    br.loc[:, "length2"] = 1 + br["end2"] - br["start2"]
    #   br[, length3 := 1 + end3 - start3]
    br.loc[:, "length3"] = 1 + br["end3"] - br["start3"]

    fai = _concatenate_br_and_fai(br, fai, prefix, regex2)
    fai = fai.loc[fai["length"] > 0, :]
    # TODO: I need to understand the roll here
    nbr = fai.loc[:, ["scaffold", "orig_scaffold", "orig_start"]].assign(orig_br=fai["orig_start"]).merge(
        nbr.loc[:, ["orig_scaffold", "orig_br"]], on=["orig_scaffold", "orig_br"], how="right",
        # roll = True?
    )
    nbr.loc[:, "br"] = nbr["orig_br"] - nbr["orig_start"] + 1
    br = nbr.loc[:, ["scaffold", "br"]]
    return br, nbr, fai


# # Break scaffolds at specified points and lift positional information to updated assembly
def break_scaffolds(breaks: pd.DataFrame, assembly: dict, prefix: str, slop: int, cores=1, species="wheat",
                    regex1=r"(^.*[^-0-9])(([0-9]+)(-[0-9]+)?$)",
                    regex2=r"(^.*[^-0-9])([0-9]+(-[0-9]+)?$)"):

    info, cov, fpairs, cssaln = assembly["info"], assembly["cov"], assembly["fpairs"], assembly["cssaln"]
    br = breaks.copy()
    fai = info.loc[:, ["scaffold", "orig_scaffold", "orig_start", "orig_end", "length"]]
    fai.loc[:, "old_scaffold"] = fai.loc[:, "scaffold"]
    print("Split scaffolds")
    j = 0
    while br.shape[0] > 0:
        j += 1
        o = fai.shape[0]
        br, nbr, fai = _calculate_breaks(br, fai, prefix, regex1=regex1, regex2=regex2, slop=slop)
        print("Iteration", j, "finished.")
        print("The number of scaffolds increased from ", o, " to {}.".format(fai.shape[0]))

    fai.loc[:, "split"] = False
    fai.loc[fai["old_scaffold"].isin(breaks["scaffold"]), "split"] = True
    fai.loc[:, "old_scaffold"] = np.nan
    assembly_new = {"info": fai}
    print("Transpose cssaln")
    z = cssaln.copy()
    z.loc[:, "scaffold_length"] = np.nan
    z.loc[:, "scaffold"] = np.nan
    z = fai.loc[:, ["scaffold", "length", "orig_scaffold", "orig_start"]].rename(
        columns={"length": "scaffold_length", "orig_start": "orig_pos"}
    ).merge(z, on=["orig_scaffold", "orig_start"], how="right")
    z.loc[:, "pos"] = z["orig_pos"] - z["orig_start"] + 1
    z.loc[:, "orig_start"] = np.nan
    assembly_new["cssaln"] = z

    assembly_new["fpairs"] = _transpose_fpairs(assembly.get("fpairs", None), fai)
    assembly_new["molecules"] = _transpose_molecules(assembly.get("molecules", None), fai)
    print("Anchor scaffolds")

    assembly_new = anchor_scaffolds(assembly_new, popseq=assembly["popseq"], species=species)
    if "mr_10x" in assembly_new["info"].columns:
        assembly_new["info"].loc[:, "mr_10x"] = np.nan
    if "mr" in assembly_new["info"].columns:
        assembly_new["info"].loc[:, "mr"] = np.nan
        assembly_new["info"].loc[:, "mri"] = np.nan
    assembly_new= _transpose_hic_cov(assembly_new, assembly, fai, info,
                                     breaks=breaks, cores=cores, prefix=prefix, regex2=regex2)
    assembly_new = _transpose_10x_cov(assembly_new, assembly, fai, info,
                                     breaks=breaks, cores=cores, prefix=prefix, regex2=regex2)
    assembly_new["binsize"] = assembly["binsize"]
    assembly_new["innerDist"] = assembly["innerDist"]
    assembly_new["minNbin"] = assembly["minNbin"]
    return assembly_new
