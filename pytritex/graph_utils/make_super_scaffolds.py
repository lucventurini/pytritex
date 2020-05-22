import pandas as pd
from pytritex.graph_utils.make_super import make_super


def make_super_scaffolds(links, info: pd.DataFrame, excluded=pd.Series([]), ncores=1, prefix=None):
    info2 = info.loc[:, ["scaffold_index", "popseq_chr", "popseq_cM", "length"]].rename(
        columns={"pospeq_chr": "chr", "popseq_cM": "cM"})
    if excluded is not None:
        excluded_scaffolds = excluded.copy()
    else:
        excluded = pd.Series([], name="scaffold_index")
        excluded_scaffolds = excluded.copy()
    input_df = info2.assign(excluded=info2["scaffold_index"].isin(excluded_scaffolds)).rename(
        columns={"scaffold_index": "cluster"})
    hl = links.copy().rename(columns={"scaffold_index1": "cluster1", "scaffold_index2": "cluster2"})
    super_scaffolds = make_super(
        hl=hl,
        cluster_info=input_df,
        verbose=False, prefix=prefix, cores=ncores,
        paths=True, path_max=0, known_ends=False, maxiter=100)
    mem_copy = super_scaffolds["membership"].copy().rename(columns={"cluster": "scaffold_index"})
    maxidx = super_scaffolds["super_info"]["super"].astype(
        str).str.replace("{}_".format(prefix), "", regex=True).astype(int).max()

    bait = ~info.scaffold_index.isin(mem_copy["scaffold_index"])
    _to_concatenate = info.loc[bait, ["scaffold_index", "popseq_chr", "popseq_cM", "length"]].rename(
        columns={"popseq_chr": "chr", "popseq_cM": "cM"}).assign(
        bin=1, rank=0, backbone=True,
        excluded=info.loc[~info.scaffold_index.isin(
            mem_copy["scaffold_index"]), "scaffold_index"].isin(excluded_scaffolds),
        super="{}_".format(prefix) + (maxidx + pd.Series(range(1, info.loc[bait].shape[0] + 1))).astype(str)
    )
    mem_copy = pd.concat([mem_copy, _to_concatenate])
    res = mem_copy.groupby("super").agg(n=("scaffold_index", "size"),
                                        nbin=("bin", "max"),
                                        max_rank=("rank", "max"),
                                        length=("length", "sum"))
    mem_copy = res.loc[:, ["n", "nbin"]].rename(columns={"n": "super_size", "nbin": "super_nbin"}).merge(
        mem_copy, left_index=True, right_on="super", how="right")

    return {"membership": mem_copy, "info": res}