#!/usr/bin/env python3
import pysam
import sys
import argparse
import textwrap
from functools import partial
import gzip
from itertools import islice
import bgzip
import time


complements = str.maketrans("ACGTNURYSWKMBVDHacgtnuryswkmbvdh", "TGCANAYRSWMKVBHDtgcanayrswmkvbhd")


def reverse_complement(seq):
    """
    Reverse complement a nucleotide sequence.
    :param seq: Sequence to be reverse complemented
    :return: A reverse complemented sequence
    """
    return seq.translate(complements)[::-1]


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("genome", type=pysam.FastaFile)
    parser.add_argument("agp", type=argparse.FileType("rt"))
    parser.add_argument("out", default=sys.stdout, type=argparse.FileType("wt"), nargs="?")
    args = parser.parse_args()

    current = None
    curr_seq = ""
    if args.out == sys.stdout:
        out = sys.stdout
    else:
        out = args.out
        # out = bgzip.BGZipWriter(args.out)
    lines = []
    cMs = []
    chrom = None
    total = 0
    for line in args.agp:
        if line[0] == "#":  # Skip comment lines
            continue
        fields = line.rstrip().split("\t")
        # super_0 1       116519  1       W       Triticum_aestivum_hereward_EIv1_scaffold_030317 1       116519  +       4A      61.282326
        ssuper, ssuper_start, ssuper_end, index, stype, reference, ref_start, ref_end, orientation = fields[:9]
        if current is None:
            current = ssuper
            if len(fields) >= 10:
                chrom, cM = fields[9], float(fields[10])
                cMs = [cM]
            else:
                current_chrom, current_cM = None, None
                cMs = []
        elif current != ssuper:
            # # Print out the scaffold            
            # lines.append(f">{current}".encode())
            # lines.extend([line.encode() for line in textwrap.wrap(curr_seq, width=60)])
            if len(cMs) > 0:
                assert chrom is not None
                cMs = sorted(cMs)
                cm_start, cm_end = cMs[0], cMs[-1]
                lines.append(">{current}_{chrom}_{cm_start}_{cm_end}".format(**locals()))
            else:
                lines.append(f">{current}")
            new = [curr_seq[n:n + 60] for n in range(0, len(curr_seq), 60)]
            # print(time.ctime(), "Adding", len(new), " lines to the total - previously", total, file=sys.stderr)
            total += len(new)
            lines.extend(new)
            if total >= 5000000:
                # print(time.ctime(), "Printing out", total, "lines", file=sys.stderr)
                total = 0
                lines = "\n".join(lines) + "\n"
                out.write(lines)
                # print(time.ctime(), "Printed out", total, "lines", file=sys.stderr)                
                lines = []
            curr_seq = ""
            current = ssuper
            chrom = current_chrom
            if len(fields) >= 10:
                chrom, cM = fields[9], float(fields[10])
                cMs = [cM]
            else:
                chrom, cMs = None, []

        if stype == "U":
            curr_seq += "N" * int(float(reference))
        elif stype == "N":  # Gap of known length
            curr_seq += "N" * int(reference)
        else:
            if len(fields) >= 10:
                current_chrom, current_cM = fields[9], float(fields[10])
                assert chrom == current_chrom
                cMs.append(current_cM)
            try:
                seq = args.genome.fetch(reference, int(ref_start) - 1, int(ref_end))
            except ValueError:
                raise ValueError("\n".join([line, ssuper, ssuper_start, ssuper_end, index, stype, reference, ref_start, ref_end, orientation]))
            if orientation == "+":
                pass
            elif orientation == "-":
                seq = reverse_complement(seq)
            else:
                raise ValueError(f"Invalid orientation: {orientation}")
            curr_seq += seq

    if current is not None:
        lines.append(f">{current}")
        lines.extend([curr_seq[n:n + 60] for n in range(0, len(curr_seq), 60)])        

    out.write(("\n".join(lines) + "\n"))
    return


if __name__ == "__main__":
    main()
