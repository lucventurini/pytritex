import pandas as pd
from pytritex.graph_utils.make_super import make_super
import dask.dataframe as dd


def make_super_scaffolds(links: str,
                         info: str,
                         save_dir: str,
                         excluded=pd.Series([]), ncores=1):
    links = dd.read_parquet(links, infer_divisions=True)
    info = dd.read_parquet(info, infer_divisions=True)
    info2 = info.loc[:, ["popseq_chr", "popseq_cM", "length"]].rename(
        columns={"popseq_chr": "chr", "popseq_cM": "cM"})
    assert "popseq_chr" not in info2.columns and "chr" in info2.columns
    assert "popseq_cM" not in info2.columns and "cM" in info2.columns
    if excluded is not None:
        excluded_scaffolds = excluded.copy()
    else:
        excluded = pd.Series([], name="scaffold_index")
        excluded_scaffolds = excluded.copy()
    input_df = info2.assign(excluded=info.index.isin(excluded_scaffolds))
    input_df.index = input_df.index.rename("cluster")
    hl = links.copy().rename(columns={"scaffold_index1": "cluster1", "scaffold_index2": "cluster2"})
    super_scaffolds = make_super(
        hl=hl,
        cluster_info=input_df,
        verbose=False, cores=ncores,
        paths=True, path_max=0, known_ends=False, maxiter=100)
    mem_copy = super_scaffolds["membership"].copy().rename(columns={"cluster": "scaffold_index"})
    maxidx = super_scaffolds["super_info"]["super"].max()

    bait = ~info.scaffold_index.isin(mem_copy["scaffold_index"])
    _to_concatenate = info.loc[bait, ["scaffold_index", "popseq_chr", "popseq_cM", "length"]].rename(
        columns={"popseq_chr": "chr", "popseq_cM": "cM"}).assign(
        bin=1, rank=0, backbone=True,
        excluded=info.loc[~info.scaffold_index.isin(
            mem_copy["scaffold_index"]), "scaffold_index"].isin(excluded_scaffolds),
        super=(maxidx + pd.Series(range(1, info.loc[bait].shape[0] + 1)))
    )
    mem_copy = pd.concat([mem_copy, _to_concatenate])
    res = mem_copy.groupby("super").agg(n=("scaffold_index", "size"),
                                        nbin=("bin", "max"),
                                        max_rank=("rank", "max"),
                                        length=("length", "sum"))
    mem_copy = res.loc[:, ["n", "nbin"]].rename(columns={"n": "super_size", "nbin": "super_nbin"}).merge(
        mem_copy, left_index=True, right_on="super", how="right")

    return {"membership": mem_copy, "info": res}
