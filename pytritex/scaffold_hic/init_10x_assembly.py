import dask.dataframe as dd
import numpy as np
import pandas as pd
from ..scaffold_10x.make_agp import make_agp
import os
import time


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
    super_molecules = dd.merge(map_10x["agp"][["scaffold_index", "super", "super_start", "super_end", "orientation"]],
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
        ["super_start", "super_end", "orientation"], axis=1).persist()
    return super_molecules


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
def transfer_cssaln(cssaln, map_10x):
    cssaln = cssaln.drop(["orig_scaffold_index", "orig_pos"], axis=1)
    cssaln = dd.merge(map_10x["agp"][["scaffold_index", "super", "super_start", "super_end", "orientation"]],
                      cssaln, on="scaffold_index")
    cssaln["pos"] = cssaln["pos"].mask(cssaln["orientation"] == 1, cssaln["super_start"] - 1 + cssaln["pos"])
    cssaln["pos"] = cssaln["pos"].mask(cssaln["orientation"] == -1, cssaln["super_end"] + 1 - cssaln["pos"])
    cssaln = cssaln.reset_index(drop=True).set_index("super").drop(["super_start", "super_end", "orientation"], axis=1)
    super_cssaln = map_10x["result"].rename(columns={"length": "super_length"}).merge(cssaln, on="super", how="right")
    super_cssaln = super_cssaln.persist()
    return super_cssaln


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
def transfer_hic(fpairs, map_10x):
    super_fpairs = fpairs[:][["scaffold_index1", "scaffold_index2", "pos1", "pos2", "chr1", "chr2"]][:]
    super_fpairs = super_fpairs.set_index("scaffold_index1").persist()
    left = map_10x["agp"][:][["scaffold_index", "super", "super_start", "super_end", "orientation"]].set_index("scaffold_index")
    left.index = left.index.rename("scaffold_index1")
    super_fpairs = dd.merge(left, super_fpairs, on="scaffold_index1", how="right")
    super_fpairs["pos1"] = super_fpairs["pos1"].mask(super_fpairs["orientation"] == 1, super_fpairs["super_start"] - 1 + super_fpairs["pos1"])
    super_fpairs["pos1"] = super_fpairs["pos1"].mask(super_fpairs["orientation"] == -1, super_fpairs["super_end"] + 1 - super_fpairs["pos1"])
    super_fpairs = super_fpairs.drop(["super_start", "super_end", "orientation"], axis=1)
    super_fpairs = super_fpairs.reset_index(drop=True).rename(columns={"super": "super1"})
    super_fpairs = super_fpairs.set_index("scaffold_index2").persist()
    left.index = left.index.rename("scaffold_index2")
    super_fpairs = dd.merge(left, super_fpairs, on="scaffold_index2", how="right")
    super_fpairs["pos2"] = super_fpairs["pos2"].mask(super_fpairs["orientation"] == 1, super_fpairs["super_start"] - 1 + super_fpairs["pos2"])
    super_fpairs["pos2"] = super_fpairs["pos2"].mask(super_fpairs["orientation"] == -1, super_fpairs["super_end"] + 1 - super_fpairs["pos2"])
    super_fpairs = super_fpairs.drop(["super_start", "super_end", "orientation"], axis=1)
    super_fpairs = super_fpairs.reset_index(drop=True).rename(columns={"super": "super2"}).persist()
    return super_fpairs


# # Initialized an assembly object with genetic map/flowsorting information, 10X and Hi-C links
# init_assembly<-function(fai, cssaln, fpairs=NULL, molecules=NULL, rename=NULL){
#
#  copy(fai)->info
#  info[, orig_start := 1]
#  info[, orig_end := length]
#  if(is.null(rename)) {
#   info[, orig_scaffold := scaffold]
#  } else {
#   setnames(info, "scaffold", "orig_scaffold")
#   rename[info, on="orig_scaffold"] -> info
#  }
#
#  copy(cssaln)->z
#  if(is.null(rename)) {
#   z[, orig_scaffold := scaffold]
#  } else {
#   setnames(z, "scaffold", "orig_scaffold")
#   rename[z, on="orig_scaffold"] -> z
#  }
#  z[, orig_pos := pos]
#  z[, orig_scaffold_length := scaffold_length]
#
#  if(!is.null(molecules)){
#   copy(molecules)->y
#   y[, orig_scaffold := scaffold]
#   y[, orig_start := start]
#   y[, orig_end := end]
#   y[, scaffold := NULL]
#   info[, .(orig_scaffold, scaffold)][y, on="orig_scaffold"]->y
#  } else {
#   y <- data.table()
#  }
#
#  if(!is.null(fpairs)){
#   copy(fpairs)->tcc
#   tcc[, orig_scaffold1 := scaffold1]
#   tcc[, orig_pos1 := pos1]
#   tcc[, orig_scaffold2 := scaffold2]
#   tcc[, orig_pos2 := pos2]
#   tcc[, scaffold1 := NULL]
#   tcc[, scaffold2 := NULL]
#   info[, .(orig_scaffold1=orig_scaffold, scaffold1=scaffold)][tcc, on="orig_scaffold1"]->tcc
#   info[, .(orig_scaffold2=orig_scaffold, scaffold2=scaffold)][tcc, on="orig_scaffold2"]->tcc
#  } else {
#   tcc <- data.table()
#  }
#
#  list(info=info, cssaln=z, fpairs=tcc, molecules=y)
# }


def _init_assembly(fai, cssaln, molecules, fpairs):
    """"""
    info = fai.copy()
    info["orig_start"] = 1
    info["orig_end"] = info["length"]
    info["orig_scaffold"] = info["scaffold"]

    ini_cssaln = cssaln[:].eval("orig_scaffold = scaffold").eval("orig_pos = pos").eval(
        "orig_scaffold_length = scaffold_length")
    if molecules is not None:
        ini_molecules = molecules.eval("""
            orig_scaffold_index = scaffold_index
            orig_start = start
            orig_end = end
        """).drop("scaffold", axis=1)
        ini_molecules = info[["orig_scaffold", "scaffold"]].merge(ini_molecules, on="orig_scaffold")
    else:
        ini_molecules = dd.from_pandas(pd.DataFrame(), chunksize=10)

    if fpairs is not None:
        tcc = fpairs[:].eval("""
            orig_scaffold_index1 = scaffold_index1
            orig_scaffold_index2 = scaffold_index2
            orig_pos1 = pos1
            orig_pos2 = pos2   
        """).drop(["scaffold_index1", "scaffold_index2"], axis=1)
    else:
        tcc = dd.from_pandas(pd.DataFrame(), chunksize=10)
    return {"info": info, "cssaln": ini_cssaln, "fpairs": tcc, "molecules": ini_molecules}


def init_10x_assembly(assembly, map_10x, gap_size=100, molecules=False, save=None):

    """Function to start the process of using HiC data for 10X-super-scaffolded assemblies.

    The purpose is to end up having objects where the super-scaffolds have been renamed as "scaffolds", so that
    we can reuse the code from before.

    :param assembly: the initial assembly from *before* the 10X super-scaffolding.
    :param map_10x: the super-scaffolded assembly
    :param molecules: boolean flag
    :param gap_size: length of gaps in the AGP
    """

    # super <- map_10x

    if "info_10x" not in assembly:
        assert os.path.exists(os.path.join(os.path.dirname(assembly["fai"]), "info_10x"))
        assembly["info_10x"] = dd.read_parquet(os.path.join(os.path.dirname(assembly["fai"]), "info_10x"))
    elif isinstance(assembly["info_10x"], str):
        assembly["info_10x"] = dd.read_parquet(assembly["info_10x"])

    if isinstance(map_10x["membership"], str):
        map_10x["membership"] = dd.read_parquet(map_10x["membership"])
    if isinstance(map_10x["result"], str):
        map_10x["result"] = dd.read_parquet(map_10x["result"])

    if isinstance(assembly["fai"], str):
        fai = dd.read_parquet(assembly["fai"])
    else:
        fai = assembly["fai"]

    print(time.ctime(), "Preparing the AGP")
    map_10x.update(make_agp(map_10x["membership"], info=fai, gap_size=gap_size))
    print(time.ctime(), "Prepared the AGP")    
    cssaln = assembly["cssaln"]
    if isinstance(cssaln, str):
        cssaln = dd.read_parquet(cssaln, infer_divisions=True)

    print(time.ctime(), "Transfering the CSS aln")
    super_cssaln = transfer_cssaln(cssaln, map_10x)
    print(time.ctime(), "Transfered the CSS aln")    
    if molecules is True:
        super_molecules = transfer_molecules(dd.read_parquet(assembly["molecules"]), map_10x)
    else:
        super_molecules = None
    if isinstance(assembly["fpairs"], str):
        assembly["fpairs"] = dd.read_parquet(assembly["fpairs"], infer_divisions=True)

    print(time.ctime(), "Transfering the HiC data")
    super_hic = transfer_hic(assembly["fpairs"], map_10x)
    print(time.ctime(), "Transfered the HiC data")    
    # fai = map_10x["agp_bed"].eval("length = bed_end - bed_start").groupby("super")["length"].sum()
    # fai.index = fai.index.rename("scaffold")
    map_10x["agp"]["super_size"] = map_10x["agp"].groupby("super")["super_index"].transform("max")
    map_10x["agp"]["super_name"] = map_10x["agp"]["super_name"].mask(map_10x["agp"]["super_size"] > 1,
                                                                     "super_" + map_10x["agp"]["super"].astype(str))

    fai = map_10x["agp"].groupby("super").agg(length=pd.NamedAgg("length", "sum"),
                                              scaffold=pd.NamedAgg("super_name", lambda values: values.iloc[0]))
    fai.index = fai.index.rename("scaffold_index")
    fai["start"] = fai["orig_start"] = 1
    fai["end"] = fai["orig_end"] = fai["length"]

    print(time.ctime(), "Initialising the assembly")
    assembly = _init_assembly(fai=fai, cssaln=super_cssaln, molecules=super_molecules, fpairs=super_hic)
    assembly["agp"] = dd.from_pandas(map_10x["agp"], chunksize=int(1e6))
    print(time.ctime(), "Initialised the assembly")    
    if save is not None:
        # return {"info": info, "cssaln": cssaln, "fpairs": tcc, "molecules": ini_molecules}
        for key, item in assembly.items():
            if item is not None:
                path = os.path.join(save, key)
                dd.to_parquet(assembly[key], path, engine="pyarrow", compression="gzip")
                assembly[key] = path

    return assembly
