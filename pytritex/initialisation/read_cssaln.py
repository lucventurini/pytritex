import pandas as pd
from pytritex.utils.chrnames import chrNames
import subprocess as sp
import dask.dataframe as dd
import numpy as np
import os
from pafpy import PafRecord
import gzip
import functools


def read_cssaln(bam: str, popseq: pd.DataFrame,
                fai: pd.DataFrame, minqual=30, minlen=10000):
    wheatchr = chrNames(species="wheat")
    # data = {"css_contig": [], "scaffold": [], "pos": []}
    proc = "samtools view -q {minqual} -F260 {bam} | cut -f 1,3,4"
    cssaln = dd.read_csv(sp.Popen(proc.format(**locals()), shell=True, stdout=sp.PIPE).stdout,
                         names=["css_contig", "scaffold", "pos"], sep="\t")
    # cssaln = pd.DataFrame(data)
    cssaln["sorted_lib"] = cssaln.css_contig.str.replace(
        "ta_contig_", "").replace("_[0-9]+$", "", regex=True)
    cssaln["sorted_alphacr"] = cssaln["sorted_lib"].str.replace("(L|S)$", "", regex=True)
    cssaln["sorted_subgenome"] = cssaln["sorted_alphacr"].str.replace("^[1-7]", "", regex=True)
    cssaln["sorted_arm"] = cssaln["sorted_lib"].str.replace("[1-7][A-D]", "", regex=True)
    cssaln = dd.merge(popseq, cssaln, left_on=["css_contig"], right_on=["css_contig"])
    cssaln = dd.merge(cssaln, wheatchr.copy().rename(
        columns={wheatchr.columns[0]: "popseq_alphachr",
                 wheatchr.columns[1]: "popseq_chr"}), left_on="popseq_alphachr", right_on="popseq_alphachr"
    )
    cssaln = dd.merge(cssaln, wheatchr.copy().rename(
        columns={wheatchr.columns[0]: "sorted_alphachr",
                 wheatchr.columns[1]: "sorted_chr"}), left_on="sorted_alphachr", right_on="sorted_alphachr"
                      )
    cssaln = dd.merge(fai[["scaffold", "length"]].rename({"length": "scaffold_length"}), cssaln,
                      left_on="scaffold", right_on="scaffold")
    cssaln = cssaln[cssaln.contig_length >= minlen]
    return cssaln


def read_morexaln_minimap(paf: str,
                          popseq: pd.DataFrame,
                          fai: pd.DataFrame,
                          save_dir: str,
                          minqual=30,
                          minlen=500, ref=None):

    if paf.endswith("gz"):
        opener = functools.partial(gzip.open, mode="rt")
    else:
        opener = functools.partial(open, mode="rt")

    css_contig, scaffold, css_length, pos = [], [], [], []

    with opener(paf) as pafbuf:
        for line in pafbuf:
            line = PafRecord.from_str(line)
            if line.mapq < minqual or line.is_primary() is False:
                continue
            if ref is None:
                if line.qname in popseq.index.values:
                    ref = True
                elif line.tname in popseq.index.values:
                    ref = False
                else:
                    raise KeyError("Neither query nor target are in the PopSeq. Aborting.")
            if ref is True:
                css_contig.append(line.qname)
                css_length.append(line.qlen)
                scaffold.append(line.tname)
                pos.append(line.tstart)
            elif ref is False:
                scaffold.append(line.qname)
                css_contig.append(line.tname)
                css_length.append(line.tlen)
                pos.append(line.qstart)

    morex = dd.from_pandas(pd.DataFrame().assign(**{
        "scaffold": scaffold,
        "css_contig": css_contig,
        "pos": pos
    }), chunksize=int(1e6))

    morex["pos"] += 1
    morex = dd.merge(popseq.set_index("css_contig"), morex, on="css_contig", how="right")
    morex = morex.reset_index(drop=True).set_index("scaffold")
    morex = dd.merge(fai[["scaffold", "length"]].reset_index(drop=False).set_index("scaffold"),
                     morex, on=["scaffold"], how="right").reset_index(drop=True).set_index("scaffold_index")
    morex = morex.query("length >= @minlen & css_length >= @minlen",
                        local_dict={"minlen": minlen})
    morex = morex[~morex.popseq_index.isna()]
    morex["orig_scaffold_index"] = morex.index.values
    morex["orig_pos"] = morex["pos"]
    fname = os.path.join(save_dir, "cssaln")
    dd.to_parquet(morex, fname, compression="gzip", compute=True, engine="pyarrow")
    return fname
