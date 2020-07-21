# PyTritex
Port to Python of Tritex. The original code was written by Martin Mascher as an R pipeline. Due to the limitations of the original code (in particular its very high memory requirements) we are now presenting a rewrite of the code in Python/Cython, using [Pandas](https://github.com/pandas-dev/pandas/) / [Dask](https://github.com/dask/dask) as our library of choice for serialising the data.

This port is not yet complete. Specifically, we are currently missing:
- the plots
- the last part of the pipeline, when HiC data is used to create the final pseudomolecules.

We might finalise the porting of the latter part of the pipeline later; for the time being, we instead recommend using remapping the scaffolds onto a reference genome (e.g. using [ragtag](https://github.com/malonge/ragtag)) and then using standard methods to find and correct misplacements, e.g. [Juicer](https://github.com/aidenlab/juicer) or [Cooltools](https://github.com/mirnylab/cooler) (with alignments provided by [Distiller](https://github.com/mirnylab/distiller-nf)).

## General structure

The pipeline has the following structure:

1. Data from the original assembly, 10X alignments, HiC alignments, and the CSS alignments is hoovered into parquet data stores (`initialisation`).
1. Using the data coming from the CSS alignments and the PopSeq reference data, each scaffold is assigned to a different chromosomal location if possible. If a scaffold is compatible with more than one location, PyTritex will keep count of the best location, the runner-up, and the proportion of both the best location towards the total, as well as the runner-up vs the best location (`anchoring`).
1. In the third step, the 10X and HiC data are analysed to derive coverage statistics for each scaffold. At this stage we are interested in obtaining a per-bin coverage: that is, determine in a windowed analysis the coverage of the data (`sequencing_coverage`).
1. In the fourth step, pytritex will try to find and break chimeric scaffolds. This is an **iterative procedure**, which is repeated until no novel breakpoints are found. Chimeric breakpoints are identified by looking for significant drops in coverage compared to the windowed baseline (obtained in point 3)
1. In the fifth step, PyTritex merges the superscaffolds according to the 10X links, filtered by two paramers:
   1. the scaffolds must have been assigned to the same chromosome
   1. the genetic distance (in *cM*) must be lower than `max_dist` (default 5 *cM*)

## Initialisation



  - An important secondary step will be to assign to each scaffold a "distance" `d` value, ie its (binned) distance from the nearest scaffold edge. The 10X and HiC coverage will be averaged for all scaffolds *based on the distance value*.
  - Using this value, we also calculate the ratio of the observed coverage vs average coverage for bins at that distance. The deviation is normalised as a `log2(ratio)` and assigned to the `r` column.
- In the fourth step, pytritex will try to find and break chimeric scaffolds. This is an **iterative procedure**, which is repeated until no novel breakpoints are found:
  - First, chimeric breakpoints are identified as follows (`find_10x_breaks`):
    - Assign each binned coverage value to a super-bin (using the `interval` parameter).
    - Exclude scaffolds that have less than `minNbin` binned values (these are deemed to be too small to infer anything).
    - Identify bins where the ratio `r` (see above) is **`ratio` (default -3)** or lower - ie, **bins where the coverage is at least 8 times lower than expected.
    - Exclude bins that are located at less than `dist` basepairs from the scaffold edge (default 5000, for wheat we use **2000**)
    - **The remaining bins with low coverage are the chimeric breakpoints.**
  - Secondly, the scaffolds are broken up using the breakpoints as a guide.
    - If we are in the second (or further) iteration, breakpoints in scaffolds that had already been broken up previously will be reconducted to the **original** scaffolds. This ensures that when the cycle ends, the ultimate scaffolds will be directly trackable to the original assembly.
    