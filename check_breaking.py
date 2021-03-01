import sys
from pytritex.chimera_breaking.calculate_broken_scaffolds import calculate_broken_scaffolds
from pytritex.chimera_breaking.find_10x_breaks import find_10x_breaks
import dask.dataframe as dd


coverage = "assembly/joblib/pytritex/chimera_breaking/1/molecule_cov_10x"
fai = "assembly/joblib/pytritex/chimera_breaking/1/fai/"

ratio=-3
interval=50000
minNbin=20
dist=2000.0
slop=200
cores=4

breaks = find_10x_breaks(dd.read_parquet(coverage, infer_divisions=True), interval=interval, minNbin=minNbin, dist=dist, ratio=ratio)
dd.to_parquet(dd.from_pandas(breaks, chunksize=10**4), "assembly/joblib/pytritex/chimera_breaking/2/breaks", compute=True, compression="gzip")
calculate_broken_scaffolds(breaks, fai=fai, slop=slop, save_dir="assembly/joblib/pytritex/chimera_breaking/2/", cores=cores)
