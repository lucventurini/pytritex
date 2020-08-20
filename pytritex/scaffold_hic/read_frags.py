import dask.dataframe as dd
import pandas as pd
import numpy as np
from ..utils import rolling_join


def read_fragdata(fai, fragfile, info, assembly_10x=None):

    frags = dd.read_csv(fragfile, names=["orig_scaffold", "start", "end"])
    frags["length"] = frags.eval("end - start")
    frags["start"] += 1
    left = fai.query("derived_from_split == False")[
        ["scaffold", "orig_scaffold_index", "to_use"]].reset_index(drop=True).rename(
        columns={"scaffold": "orig_scaffold"}).set_index("orig_scaffold")
    frags = dd.merge(left, frags, on="orig_scaffold").reset_index(drop=True).set_index("orig_scaffold_index")
    # Select scaffolds that have not been broken up, we'll use them as-is
    frags_to_keep = frags.query("to_use == True")[:].drop("to_use", axis=1)
    frags_to_keep.index = frags_to_keep.index.rename("scaffold_index")
    # Roll-join the others
    frags_to_roll = frags.query("to_use == False").drop("to_use", axis=1)
    left = fai.query("to_use == True & orig_scaffold_index in @rolled",
                     local_dict={
                         "rolled": frags_to_roll["orig_scaffold_index"].unique().values.compute().tolist()})
    left = left[["orig_scaffold_index", "orig_start"]].reset_index(drop=False).set_index("orig_scaffold_index")
    left["start"] = left["orig_start"]
    frags_to_roll = rolling_join(left, frags_to_roll, on="orig_scaffold_index",
                                 by="start").reset_index(drop=False).set_index("scaffold_index")
    frags_to_roll["end"] -= frags_to_roll["orig_start"] + 1
    frags_to_roll["start"] -= frags_to_roll["orig_start"] + 1
    frags_to_roll = frags_to_roll.drop(["orig_start", "orig_scaffold_index"], axis=1)

    # Now concatenate for the final frag file.
    frags = dd.concat([frags_to_keep, frags_to_roll]).astype({"start": np.int32, "end": np.int32})
    if assembly_10x is not None:
        raise NotImplementedError()
    else:
        nfrags = frags.groupby("scaffold_index").size()
        