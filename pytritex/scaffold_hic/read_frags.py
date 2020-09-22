import dask.dataframe as dd
import pandas as pd
import numpy as np
from ..utils import rolling_join
import os
from dask.delayed import delayed

# # Read BED files with positions of restriction fragments on the input assembly, lift coordinates to updated assemblies
# read_fragdata<-function(info, map_10x=NULL, assembly_10x=NULL, file){
#  fragbed<-fread(file, head=F, col.names=c("orig_scaffold", "start", "end"))
#  fragbed[, length := end - start]
#  fragbed[, start := start + 1]
#  info[, .(scaffold, start=orig_start, orig_start, orig_scaffold)][fragbed, on=c("orig_scaffold", "start"), roll=T]->fragbed
#  fragbed[, start := start - orig_start + 1]
#  fragbed[, end := end - orig_start + 1]
#  fragbed[, orig_start := NULL]
#  fragbed[, orig_scaffold := NULL]
#  if(!is.null(assembly_10x)){
#   map_10x$agp[gap == F, .(super, orientation, super_start, super_end, scaffold)][fragbed, on="scaffold"]->fragbed
#   fragbed[orientation == 1, start := super_start - 1 + start]
#   fragbed[orientation == 1, end := super_start - 1 + end]
#   fragbed[orientation == -1, start := super_end - end + 1]
#   fragbed[orientation == -1, end := super_end - start + 1]
#   fragbed[, c("orientation", "super_start", "super_end", "scaffold") := list(NULL, NULL, NULL, NULL)]
#   setnames(fragbed, "super", "orig_scaffold")
#
#   assembly_10x$info[, .(scaffold, start=orig_start, orig_start, orig_scaffold)][fragbed, on=c("orig_scaffold", "start"), roll=T]->fragbed
#   fragbed[, start := start - orig_start + 1]
#   fragbed[, end := end - orig_start + 1]
#   fragbed[, orig_start := NULL]
#   fragbed[, orig_scaffold := NULL]
#
#   fragbed[, .(nfrag = .N), keyby=scaffold][assembly_10x$info, on="scaffold"][is.na(nfrag), nfrag := 0]->z
#  } else {
#   fragbed[, .(nfrag = .N), keyby=scaffold][info, on="scaffold"][is.na(nfrag), nfrag := 0]->z
#  }
#  list(bed=fragbed[], info=z[])
# }


def read_fragdata(fai, fragfile, map_10x, savedir=None):
    
    if fragfile.endswith(".gz"):
        import subprocess as sp
        if savedir is None:
            dirname = "."
        else:
            dirname = savedir
        sp.call("gunzip -c {} > {}".format(fragfile, os.path.join(dirname, "frags.csv")), shell=True)
        fragfile = os.path.join(dirname, "frags.csv")        

    if isinstance(fai, str):
        fai = dd.read_parquet(fai, infer_divisions=True)
    
    left = fai.query("derived_from_split == False")[
        ["scaffold", "orig_scaffold_index", "to_use"]].reset_index(drop=True).rename(
        columns={"scaffold": "orig_scaffold"}).set_index("orig_scaffold")

    frags = dd.read_csv(fragfile, names=["orig_scaffold", "start", "end"], delimiter="\t", header=None,
                        dtype={"orig_scaffold": str, "start": np.int32, "end": np.int32})
    frags["length"] = frags["end"] - frags["start"]
    frags["start"] += 1    
    
    frags = dd.merge(left, frags, on="orig_scaffold").reset_index(drop=True).set_index("orig_scaffold_index")
    if savedir is not None:
        path = dd.to_parquet(frags, os.path.join(savedir, "raw_frags"),
                             compression="gzip", compute=True, engine="pyarrow")
    
    # Select scaffolds that have not been broken up, we'll use them as-is
    frags_to_keep = frags.query("to_use == True")[:].drop("to_use", axis=1)
    frags_to_keep.index = frags_to_keep.index.rename("scaffold_index")
    # Roll-join the others
    frags_to_roll = frags.query("to_use == False").drop("to_use", axis=1)
    left = fai.query("to_use == True & orig_scaffold_index in @rolled",
                     local_dict={
                         "rolled": np.unique(frags_to_roll.index.values.compute()).tolist()})
    left = left[["orig_scaffold_index", "orig_start"]].reset_index(drop=False).set_index("orig_scaffold_index")
    left["start"] = left["orig_start"]
    frags_to_roll = rolling_join(left, frags_to_roll, on="orig_scaffold_index",
                                 by="start").reset_index(drop=False).set_index("scaffold_index")
    frags_to_roll["end"] -= frags_to_roll["orig_start"] + 1
    frags_to_roll["start"] -= frags_to_roll["orig_start"] + 1
    frags_to_roll = frags_to_roll.drop(["orig_start", "orig_scaffold_index"], axis=1)

    # Now concatenate for the final frag file.
    frags = dd.concat([frags_to_keep, frags_to_roll]).astype({"start": np.int32, "end": np.int32})
    left = map_10x["agp"].query("gap == False")[["super", "orientation", "super_start", "super_end"]]
    fragbed = dd.merge(left, frags, on="scaffold_index")
    # map_10x$agp[gap == F, .(super, orientation, super_start, super_end, scaffold)][fragbed, on="scaffold"]->fragbed
    bait = (fragbed["orientation"] == 1)
    fragbed.loc[bait, ["start"]] = fragbed.loc[bait]["super_start"] - 1 + fragbed.loc[bait]["start"]
    fragbed.loc[bait, ["end"]] = fragbed.loc[bait]["super_start"] - 1 + fragbed.loc[bait]["end"]
    bait = (fragbed["orientation"] == -1)
    fragbed.loc[bait, ["start"]] = fragbed.loc[bait]["super_end"] + 1 - fragbed.loc[bait]["end"]
    fragbed.loc[bait, ["end"]] = fragbed.loc[bait]["super_end"] + 1 - fragbed.loc[bait]["start"]
    fragbed = fragbed.drop(["orientation", "super_start", "super_end", "scaffold"], axis=1)
    fragbed = fragbed.rename(columns={"super": "orig_scaffold"})
    # TODO: is there "info" in map_10x?
    left = map_10x["info"][["orig_start", "orig_scaffold_index"]].reset_index(drop=False).set_index(
        "orig_scaffold_index")
    left["start"] = left["orig_start"]
    # assembly_10x$info[,
    # .(scaffold, start=orig_start, orig_start, orig_scaffold)][fragbed, on=c("orig_scaffold", "start"), roll=T]->fragbed
    fragbed = rolling_join(left, fragbed, on="orig_scaffold_index", by="start")
    #   fragbed[, start := start - orig_start + 1]
    #   fragbed[, end := end - orig_start + 1]
    #   fragbed[, orig_start := NULL]
    #   fragbed[, orig_scaffold := NULL]
    #
    #   fragbed[, .(nfrag = .N), keyby=scaffold][assembly_10x$info, on="scaffold"][is.na(nfrag), nfrag := 0]->z
    fragbed = fragbed.eval("""start = start - orig_start + 1
    end = end - orig_start + 1""").drop(["orig_start"]).reset_index(drop=True).set_index("scaffold_index")
    info = fragbed.groupby("scaffold_index").size().to_frame("nfrag")
    info = info.merge(map_10x["info"], on="scaffold_index", how="right")
    info["nfrag"] = info["nfrag"].fillna(0)
    if savedir is not None:
        fragbed_name = os.path.join(savedir, "fragbed")
        dd.to_parquet(fragbed, fragbed_name, compression="gzip", engine="pyarrow")
        fraginfo_name = os.path.join(savedir, "frag_info")
        dd.to_parquet(info, fraginfo_name, compression="gzip", engine="pyarrow")
        return {"bed": fragbed_name, "info": fraginfo_name}
    else:
        return {"bed": fragbed, "info": info}
