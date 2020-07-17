# PyTritex
Port to Python of Tritex. The original code was written by Martin Mascher as an R pipeline. Due to the limitations of the original code (in particular its very high memory requirements) we are now presenting a rewrite of the code in Python/Cython, using [Pandas](https://github.com/pandas-dev/pandas/) / [Dask](https://github.com/dask/dask) as our library of choice for serialising the data.

This port is not yet complete. Specifically, we are currently missing:
- the plots
- the last part of the pipeline, when HiC data is used to create the final pseudomolecules.

The latter part is due to our 

## General structure

The 
