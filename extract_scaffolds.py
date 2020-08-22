#!/usr/bin/env python

import sys
import pandas as pd
import pysam
import argparse
import joblib
import itertools
from pytritex.utils.chrnames import chrNames
from numpy import nan
import textwrap
import subprocess as sp
from time import ctime


pop_columns = ["popseq_index", "sorted_alphachr",  "sorted_chr", "popseq_alphachr",  "popseq_chr",
               "css_contig",  "css_contig_length",  "popseq_cM", "sorted_lib",
               "sorted_subgenome", "sorted_arm"]

valid_chroms = set(["".join([str(_) for _ in a]) for a in itertools.product(list(range(1, 8)), ["A", "B", "D"])])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-ot", "--out-table", required=True)
    parser.add_argument("-og", "--out-genome", required=True)
    parser.add_argument("agp", type=argparse.FileType("rt"))
    parser.add_argument("genome", type=pysam.FastaFile)
    args = parser.parse_args()

    sorted_alphachr, css_contig, css_contig_length, popseq_cM, sorted_subgenome = [], [], [], [], []
    cs42_chrom, chrom_starts, chrom_ends = [], [], []

    # out_genome = pysam.BGZFile(args.out_genome, mode="wb")
    out_genome = open(args.out_genome, mode="wb")
    genome_lines = []

    current_seq, current_chrom = None, None

    total = 0
    for line in args.agp:
        if line[0] == "#":
            continue
        fields = line.rstrip().split("\t")
        if fields[4] != "W":
            continue
        # chr1A   1       1437363 1       W       chr1A_super1    scaffold38755|1 1       1437363 ?       1A      1.20302702702703
        chrom, cstart, cend, sup_ind, ttype, name, scaf, start, end, strand, pop_chrom, pop_cM = fields
        cstart, cend = int(cstart), int(cend)
        if chrom != current_chrom:
            current_chrom = chrom[:]
            current_seq = args.genome.fetch(chrom)
            print(ctime(), "Getting", chrom, file=sys.stderr)

        cstart, cend = int(cstart), int(cend)
        seq = textwrap.wrap(args.genome.fetch(chrom, cstart - 1, cend), 60)
        seq = current_seq[cstart - 1:cend]
        assert len(seq) == cend - cstart + 1
        ls = len(seq)
        seq = [seq[x:max(x + 60, ls)] for x in range(0, ls, 60)]
        seq  = [seq]
        total += ls
        print(ctime(), "Extracted", scaf, file=sys.stderr)
        genome_lines.append(">{}".format(scaf))
        genome_lines.extend(seq)
        if total >= 5 * 10**6:
            print(ctime(), "Writing", len(genome_lines), "to the output file", file=sys.stderr)
            genome_lines = "\n".join(genome_lines) + "\n"
            out_genome.write(genome_lines.encode())
            out_genome.flush()
            print(ctime(), "Written", len(genome_lines), "to the output file", file=sys.stderr)
            genome_lines = []
            total = 0

        if pop_chrom in ("NA", "?"):
            pop_chrom = nan
            subgenome = nan
        else:
            if pop_chrom not in valid_chroms:
                pop_chrom = nan
                subgenome = nan
            else:
                subgenome = pop_chrom[1]
                assert subgenome in ("A", "B", "D")
        if pop_cM in ("NA", "?"):
            pop_cM = nan
        else:
            try:
                pop_cM = float(pop_cM)
            except ValueError:
                pop_cM = nan
        sorted_alphachr.append(pop_chrom)
        css_contig.append(name)
        css_contig_length.append(int(end) - int(start) + 1)
        popseq_cM.append(pop_cM)
        sorted_subgenome.append(subgenome)
        cs42_chrom.append(chrom)
        chrom_starts.append(cstart)
        chrom_ends.append(cend)

    genome_lines = ("\n".join(genome_lines) + "\n")
    out_genome.write(genome_lines.encode())
    out_genome.flush()
    out_genome.close()
    popseq_df = pd.DataFrame().assign(sorted_alphachr=sorted_alphachr,
                                      css_contig=css_contig, css_contig_length=css_contig_length,
                                      popseq_cM=popseq_cM,
                                      sorted_subgenome=sorted_subgenome,
                                      cs42_chrom=cs42_chrom, chrom_start=chrom_starts, chrom_end=chrom_ends)
    
    popseq_df.loc[:, "popseq_alphachr"] = popseq_df.loc[:, "sorted_alphachr"].values
    names = chrNames()
    names = names.rename(columns={"chr": "sorted_chr", "alphachr": "sorted_alphachr"})
    popseq_df = popseq_df.merge(names, on="sorted_alphachr", how="left")
    popseq_df.loc[:, "popseq_chr"] = popseq_df.loc[:, "sorted_chr"].values
    popseq_df.loc[:, "sorted_lib"] = popseq_df.loc[:, "popseq_alphachr"].values
    popseq_df.loc[:, "sorted_arm"] = nan
    popseq_df.loc[:, "popseq_index"] = list(range(1, popseq_df.shape[0] + 1))
    popseq_df = popseq_df[pop_columns + ["cs42_chrom", "chrom_start", "chrom_end"]]
    popseq_df = popseq_df.sort_values(["popseq_chr", "popseq_cM"])
    joblib.dump(popseq_df, args.out_table, compress=("gzip", 6))
    return


main()
