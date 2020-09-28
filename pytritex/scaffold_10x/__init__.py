from ..graph_utils.branch_remover import _initial_branch_remover
from .link_finder import _initial_link_finder
from .orient_scaffolds import orient_scaffolds
from pytritex.graph_utils.tip_removal.tip_remover import remove_tips
from .scaffold_unanchored import _scaffold_unanchored
import time
from joblib import Memory
import os
import hashlib
from dask.distributed import Client
import logging
logger = logging.getLogger("dask")
import json
import dask.dataframe as dd
from sys import stdout
from io import FileIO
import pandas as pd


def scaffold_10x(assembly: dict, memory: Memory, save_dir: str,
                 client: Client,
                 min_npairs=5, max_dist=1e5, min_nmol=6,
                 min_nsample=2, popseq_dist=5, max_dist_orientation=5,
                 ncores=1, verbose=True, raw=False, unanchored=True):
    """Construct super-scaffold from 10X links. Use heuristics based on genetic map information
    to prune erroneous edges in the scaffold graph."""

    info, molecules, fai = assembly["info"], assembly["molecules"], assembly["fai"]
    # print(assembly.keys())
    print(time.ctime(), "Finding initial links")

    sha = hashlib.sha256()
    params = (min_npairs, max_dist, min_nmol, min_nsample, popseq_dist, max_dist_orientation, raw)
    sha.update(str(params).encode())
    hash_string = sha.hexdigest()
    del sha
    logger.warning("%s Starting parameters %s, digest folder %s",
                                      time.ctime(), params, hash_string)
    folder = os.path.join(save_dir, "joblib", "pytritex", "scaffold_10x", hash_string)
    os.makedirs(folder, exist_ok=True)
    # Dump the parameters
    params_dict = {"min_npairs": int(min_npairs), "max_dist": int(max_dist), "min_nmol": int(min_nmol),
                   "min_nsample": int(min_nsample), "popseq_dist": float(popseq_dist), "max_dist_orientation": int(max_dist_orientation),
                   "raw": bool(raw), "info": str(info), "molecules": str(molecules), "fai": str(fai), "unanchored": bool(unanchored)}
    with open(os.path.join(folder, "meta.json"), "wt") as meta_out:
        json.dump(params_dict, meta_out)
    
    sample_count, links, link_pos = memory.cache(_initial_link_finder,
                                                 ignore=["client"])(
        info=info, molecules=molecules,
        save_dir=folder,
        client=client,
        fai=fai, verbose=verbose, popseq_dist=popseq_dist,
        min_npairs=min_npairs, max_dist=max_dist, min_nmol=min_nmol,
        min_nsample=min_nsample)
    print(time.ctime(), "Found initial links")
    excluded = set()
    print(time.ctime(), "Starting initial pruning")
    out, excluded = memory.cache(_initial_branch_remover, ignore=["ncores", "client"])(
        client=client, save_dir=folder, links=links, info=info, excluded=excluded, ncores=ncores)
    res = out["info"]
    membership = out["membership"]
    print(time.ctime(), "Finished initial pruning")

    if verbose is True:
        print("Finding initial super-scaffolds")

    if raw is False:

        membership, res, excluded = memory.cache(remove_tips, ignore=["client", "ncores", "verbose"])(
            links=links, excluded=excluded,
            out=out, info=info,
            client=client,
            save_dir=os.path.join(folder, "tip_removal"),
            ncores=ncores, verbose=verbose,
            min_dist=1e4)
        if popseq_dist > 0 and unanchored is True:

            membership, res = memory.cache(_scaffold_unanchored,
                                           ignore=["client", "ncores", "verbose"])(
                links,
                excluded, membership, info,
                sample_count,
                popseq_dist=popseq_dist,
                save_dir=os.path.join(folder, "unanchored"), client=client,
                ncores=1, verbose=False)

    membership, result = memory.cache(orient_scaffolds, ignore=["client"])(
        info=info, res=res, membership=membership,
        link_pos=link_pos, max_dist_orientation=max_dist_orientation,
        save_dir=os.path.join(folder, "orientation"),
        client=client)

    return membership, result


def print_agp(agp, fai, out=None):

    # agp = data[["super", "super_start", "super_end", "super_index", "gap",
    #                 "orig_scaffold_index", "orig_start", "orig_end", "orientation", "alphachr", "cM",
    #                 "scaffold_index", "length", "super_length", "super_size"]][:].persist()

    if isinstance(fai, str):
        fai = dd.read_parquet(fai)

    if isinstance(agp, str):
        agp = dd.read_parquet(agp)

    orig_scaffolds = fai.query("derived_from_split == False")[["scaffold"]]
    orig_scaffolds.index = orig_scaffolds.index.rename("orig_scaffold_index")

    # Columns:
    # 1 Name, 2 begin, 3 end, 4 part num, 5 Type, 6 [[ID OR gap length]],
    # 7 [[Comp beg OR Gap_type (scaffold)]], 8 [[End OR evidence (yes/no)]],
    # 9 [[Orientation OR gap evidence (paired-ends)]],
    # 10 chromosome, 11 cM

    agp = agp.merge(orig_scaffolds, on="orig_scaffold_index", how="left").persist()
    agp["type"] = agp["gap"].map({True: "U", False: "W"})
    agp["super_name"] = agp["scaffold"].mask(agp["super_size"] > 1,
                                             "super_" + agp["super"].astype(str))
    agp["start"] = agp["orig_start"].mask(agp["scaffold_index"] == -1, "scaffold")
    agp["end"] = agp["orig_end"].mask(agp["scaffold_index"] == -1, "yes")
    agp["strand"] = agp["orientation"].mask(agp["scaffold_index"] == -1, "paired-ends")
    agp["strand"] = agp["orientation"].map({1: "+", -1: "-", "paired-ends": "paired-ends"})

    to_print = agp[["super_name", "super_start", "super_end", "super_index", "type", "scaffold",
                    "start", "end", "strand", "alphachr", "cM"]]
    if out is None:
        out = stdout
    elif isinstance(out, str):
        out = open(out, "wt")
    else:
        assert isinstance(out, FileIO)

    print("##agp-version	2.0", file=out)
    to_print = to_print.compute()
    to_print.to_csv(out, sep="\t", index=False)
    return
