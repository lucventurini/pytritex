import pandas as pd


def make_super_scaffolds(links, info: pd.DataFrame, excluded=(), ncores=1, prefix=None):
    info2 = info[["scaffold", "popseq_chr", "popseq_cM", "length"]][:].rename(columns={
        "popseq_chr": "chr", "popseq_cM": "cM"})
    excluded_scaffolds = excluded
    input_df = info2[info2.scaffold.isin(excluded_scaffolds)].rename(
        columns={"scaffold": "cluster"}
    )
    hl = links.copy().rename(columns={"scaffold1": "cluster1", "scaffold2": "cluster2"})
    s = make_super(hl=hl, cluster_info=input_df, verbose=False, prefix=prefix, cores=ncores,
                   paths=True, path_max=0, known_ends=False, maxiter=100)
    m = s["membership"][:].rename(columns={"cluster": "scaffold"})
    maxidx = s["super_info"]["super"].astype(str).str.replace("^{}_".format(prefix), "", regex=True).astype(int).max()
    binder= info[~info["scaffold"].isin(m["scaffold"])][:][["scaffold", "length" "popseq_chr", "popseq_cM"]].rename(
        columns={"popseq_chr": "chr", "popseq_cM": "cM"})
    binder["bin"], binder["rank"], binder["backbone"] = 1, 0, True
    binder["excuded"] = binder["scaffold"].isin(excluded_scaffolds)
    binder["super"] = pd.Series(list(range(maxidx + 1, maxidx + 1 + binder.shape[0]))).astype(str).replace(
        "^(.)", "{}_\\1".format(prefix), regex=True)
    m = m.append(binder)
    res = m[:].groupby("super").agg({
    "bin": ["count", "max"], "rank": "max", "length": "sum"}).rename(
        columns={("bin", "count"): "n", ("bin", "max"): "nbin", "rank": "max_rank"}
    ).reset_index()
    m = pd.merge(m, res[["super", "n", "nbin"]].rename(columns={"n": "super_size", "nbin": "super_nbin"}),
                 left_on="super", right_on="super")
    return {"membership": m, "info": res}


