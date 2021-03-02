import itertools
import pandas as pd


def chrNames(agp=False, species="wheat"):
    if species == "wheat":
        chrom_list = list(itertools.chain(*["{chrom}A {chrom}B {chrom}D".format(
                    chrom=chrom).split() for chrom in range(1, 8)]))
    elif species == "barley":
        chrom_list = list(itertools.chain(*["{chrom}H".format(
                    chrom=chrom) for chrom in range(1, 8)]))
    elif species == "rye":
        chrom_list = list(itertools.chain(*["{chrom}R".format(
            chrom=chrom) for chrom in range(1, 8)]))
    elif species == "lolium":
        chrom_list = list(itertools.chain(*["{chrom}".format(
            chrom=chrom) for chrom in range(1, 8)]))
    elif species == "sharonensis":
        chrom_list = list(itertools.chain(*["{chrom}S".format(
            chrom=chrom) for chrom in range(1, 8)]))
    elif species == "oats":
        chrom_list = list(itertools.chain(*["{chrom}M".format(
            chrom=chrom) for chrom in range(1, 8)]))
    else:
        raise KeyError(
            "Parameter 'species' is not valid. Please set 'species' to one of"
"\"wheat\", \"barley\", \"oats\", \"lolium\", \"sharonensis\" or \"rye\".")
    if agp is True:
        chrom_list = ["Un"] + chrom_list
        idx_start = 0
    else:
        idx_start = 1
    chroms = pd.DataFrame({
        "chr": range(idx_start, len(chrom_list) + idx_start),
        "alphachr": chrom_list
    })
    return chroms


wheatchr = chrNames
