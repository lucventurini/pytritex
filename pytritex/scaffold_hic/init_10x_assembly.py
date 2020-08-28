import dask.dataframe as dd
import numpy as np
import pandas as pd
from ..scaffold_10x.make_agp import make_agp


def transfer_molecules(molecules, map_10x):
    """

    :return:
    """
    # if(molecules){
    #   copy(assembly$molecules)->z
    #   z[, c("orig_scaffold", "orig_start") := list(NULL, NULL)]
    #   super$agp[, .(scaffold, super, super_start, super_end, orientation)][z, on="scaffold"]->z
    #   z[orientation == 1, start := super_start - 1 + start]
    #   z[orientation == 1, end := super_start - 1 + end]
    #   z[orientation == -1, start := super_end - end + 1]
    #   z[orientation == -1, end := super_end - start + 1]
    #   z[, scaffold := NULL]
    #   z[, c("super_start", "super_end", "orientation") := list(NULL, NULL, NULL)]
    #   setnames(z, "super", "scaffold")
    #   z -> s_molecules
    #  } else {
    #   s_molecules <- NULL
    #  }

    super_molecules = molecules[:].drop(["orig_scaffold_index", "orig_start"], axis=1)
    super_molecules = dd.merge(map_10x["agp"][["super", "super_start", "super_end", "orientation"]],
                               super_molecules, on="scaffold_index", how="right")
    super_molecules["start"] = super_molecules["start"].mask(
        super_molecules["orientation"] == 1, super_molecules["super_start"] - 1 + super_molecules["start"])
    super_molecules["end"] = super_molecules["end"].mask(
        super_molecules["orientation"] == 1, super_molecules["super_start"] - 1 + super_molecules["end"])
    super_molecules["start"] = super_molecules["start"].mask(
        super_molecules["orientation"] == -1, super_molecules["super_end"] + 1 + super_molecules["end"])
    super_molecules["end"] = super_molecules["end"].mask(
        super_molecules["orientation"] == -1, super_molecules["super_end"] + 1 - super_molecules["start"])
    super_molecules = super_molecules.reset_index(drop=True).set_index("super").drop(
        ["super_start", "super_end", "orientation"], axis=1)
    return super_molecules



def transfer_cssaln(cssaln, map_10x):
    """"""

    #  copy(assembly$cssaln)->z
    #  z[, orig_scaffold := NULL]
    #  z[, orig_scaffold_length := NULL]
    #  z[, orig_pos := NULL]
    #  z[, scaffold_length := NULL]
    #  super$agp[, .(scaffold, super, super_start, super_end, orientation)][z, on="scaffold"]->z
    #  z[orientation == 1, pos := super_start - 1 + pos]
    #  z[orientation == -1, pos := super_end - pos + 1]
    #  z[, scaffold := NULL]
    #  setnames(z, "super", "scaffold")
    #  z[, c("super_start", "super_end", "orientation") := list(NULL, NULL, NULL)]
    #  super$info[, .(scaffold=super, scaffold_length=length)][z, on="scaffold"]->z
    #  z->s_cssaln

    cssaln = cssaln.drop(["orig_scaffold", "orig_scaffold_length", "orig_pos", "scaffold_length"], axis=1)
    cssaln = dd.merge(map_10x["agp"][["super", "super_start", "super_end", "orientation"]],
                      cssaln, on="scaffold_index")
    cssaln["pos"] = cssaln["pos"].mask(cssaln["orientation"] == 1, cssaln["super_start"] - 1 + cssaln["pos"])
    cssaln["pos"] = cssaln["pos"].mask(cssaln["orientation"] == -1, cssaln["super_end"] + 1 - cssaln["pos"])
    cssaln = cssaln.reset_index(drop=True).set_index("super").drop(
        ["super_start", "super_end", "orientation"], axis=1)
    super_cssaln = map_10x["info"].rename(columns={"length": "super_length"}).merge(cssaln, on="super", how="right")
    return super_cssaln


def transfer_hic(fpairs, map_10x):

    """"""
    # if(!is.null(assembly$fpairs) && nrow(assembly$fpairs) > 0){
    #   copy(assembly$fpairs)->z
    #   z[, c("orig_scaffold1", "orig_pos1") := list(NULL, NULL)]
    #   z[, c("orig_scaffold2", "orig_pos2") := list(NULL, NULL)]
    #   z[, c("chr1", "chr2") := list(NULL, NULL)]
    #   super$agp[, .(scaffold1=scaffold, super, super_start, super_end, orientation)][z, on="scaffold1"]->z
    #   z[orientation == 1, pos1 := super_start - 1 + pos1]
    #   z[orientation == -1, pos1 := super_end - pos1 + 1]
    #   z[, scaffold1 := NULL]
    #   setnames(z, "super", "scaffold1")
    #   z[, c("super_start", "super_end", "orientation") := list(NULL, NULL, NULL)]
    #   super$agp[, .(scaffold2=scaffold, super, super_start, super_end, orientation)][z, on="scaffold2"]->z
    #   z[orientation == 1, pos2 := super_start - 1 + pos2]
    #   z[orientation == -1, pos2 := super_end - pos2 + 1]
    #   z[, scaffold2 := NULL]
    #   setnames(z, "super", "scaffold2")
    #   z[, c("super_start", "super_end", "orientation") := list(NULL, NULL, NULL)]
    #  } else {
    #   z <- NULL
    #  }

    super_fpairs = fpairs[:].drop(["orig_scaffold_index1", "orig_pos1", "orig_scaffold_index2", "orig_pos2",
                                   "chr1", "chr2"], axis=1)
    super_fpairs = super_fpairs.set_index("scaffold_index1").persist()
    left = map_10x["agp"][:][["super", "super_start", "super_end", "orientation"]]
    left.index = left.index.rename("scaffold_index1")
    super_fpairs = dd.merge(left, super_fpairs, on="scaffold_index1", how="right")
    super_fpairs["pos1"] = super_fpairs["pos1"].mask(super_fpairs["orientation"] == 1,
                                                     super_fpairs["super_start"] - 1 + super_fpairs["pos1"])
    super_fpairs["pos1"] = super_fpairs["pos1"].mask(super_fpairs["orientation"] == -1,
                                                     super_fpairs["super_end"] + 1 - super_fpairs["pos1"])
    super_fpairs = super_fpairs.drop(["super_start", "super_end", "orientation"], axis=1)
    super_fpairs = super_fpairs.reset_index(drop=True).rename(columns={"super": "super1"})
    super_fpairs = super_fpairs.set_index("scaffold_index2").persist()
    left.index = left.index.rename("scaffold_index2")
    super_fpairs = dd.merge(left, super_fpairs, on="scaffold_index2", how="right")
    super_fpairs["pos2"] = super_fpairs["pos2"].mask(super_fpairs["orientation"] == 1,
                                                     super_fpairs["super_start"] - 1 + super_fpairs["pos2"])
    super_fpairs["pos2"] = super_fpairs["pos2"].mask(super_fpairs["orientation"] == -1,
                                                     super_fpairs["super_end"] + 1 - super_fpairs["pos2"])
    super_fpairs = super_fpairs.drop(["super_start", "super_end", "orientation"], axis=1)
    super_fpairs = super_fpairs.reset_index(drop=True).rename(columns={"super": "super2"})
    return super_fpairs


def init_10x_assembly(assembly, map_10x, gap_size=100, molecules=False):

    """Function to start the process of using HiC data for 10X-super-scaffolded assemblies.
    :param assembly: the initial assembly from *before* the 10X super-scaffolding.
    :param map_10x: the super-scaffolded assembly
    :param molecules: boolean flag
    :param gap_size: length of gaps in the AGP
    """

    # super <- map_10x

    map_10x.update(make_agp(map_10x["membership"], info=assembly["info_10x"], gap_size=gap_size))
    cssaln = assembly["cssaln"]
    if isinstance(cssaln, str):
        cssaln = dd.read_parquet(cssaln, infer_divisions=True)
    super_cssaln = transfer_cssaln(cssaln, map_10x)
    super_molecules = transfer_molecules(assembly["molecules"], map_10x)
    super_hic = transfer_hic(assembly["fpairs"], map_10x)
    fai = map_10x["agp_bed"].eval("length = bed_end - bed_start")[["length", "super"]]

    # init_assembly(fai=super$agp_bed[, .(length=sum(bed_end - bed_start)), key=.(scaffold=super)],
    # cssaln=s_cssaln, molecules=s_molecules, fpairs=z)
