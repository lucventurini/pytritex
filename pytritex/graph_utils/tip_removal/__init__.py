import pandas as pd
import dask.dataframe as dd

minimum_distance = 1e4


def _calculate_degree(links: str, excluded) -> pd.DataFrame:
    links = dd.read_parquet(links, infer_divisions=True)
    degree = links[~links["scaffold_index2"].isin(excluded)].groupby(
        "scaffold_index1").size().to_frame("degree").compute()
    degree.index = degree.index.rename("scaffold_index")
    return degree
