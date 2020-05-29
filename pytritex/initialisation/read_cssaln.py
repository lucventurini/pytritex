import pandas as pd
from pytritex.utils.chrnames import chrNames
import subprocess as sp
import dask.dataframe as dd
import numpy as np


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
                          minqual=30, minlen=500,
                          ref=True,
                          prefix=True):
    command = "zgrep tp:A:P {paf} | awk -v l={minlen} -v q={minqual} '$2>=l && $12>=q'"
    if ref is True:
        names = ["css_contig", "scaffold", "pos"]
        cols = [0, 5, 7]
    else:
        names= ["scaffold", "pos", "css_contig"]
        cols = [0, 2, 5]
    dtypes = {"css_contig": str, "scaffold": str, "pos": np.int32}
    morex = dd.from_pandas(pd.read_csv(
        sp.Popen(command.format(paf=paf, minqual=minqual, minlen=minlen), shell=True, stdout=sp.PIPE).stdout,
        sep="\t", names=names, usecols=cols, dtype=dtypes),
        npartitions=10)
    if prefix is True:
        morex["css_contig"] = morex.css_contig.str.replace("^", "morex_")
    morex["pos"] = pd.to_numeric(morex["pos"].compute() + 1, downcast="signed")
    morex = dd.merge(popseq, morex, on="css_contig", how="right").drop("css_contig", axis=1)
    morex = dd.merge(fai[["scaffold", "scaffold_index", "length"]],
        morex, left_on=["scaffold"], right_on=["scaffold"], how="right").drop("scaffold", axis=1)
    morex["orig_scaffold_index"] = morex["scaffold_index"]
    morex["orig_pos"] = morex["pos"]
    return morex
