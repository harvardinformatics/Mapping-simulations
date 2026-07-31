#!/usr/bin/env python3
"""
consensus_diff_snps.py -- emit the SNP differences between an iteration consensus and
the reference as a VCF, for the s01 baseline-vs-iterative-mapping analysis.

In the SNP-only baseline the consensus is the SAME LENGTH as the reference and shares
its coordinates (substitutions do not shift positions -- verified: the bcftools
consensus chain is a single gapless block). So the pipeline's recovered divergence at
any iteration is simply the set of positions where the consensus base differs from the
reference base -- an EXACT position-wise comparison, with no alignment and no second
variant caller involved.

Each differing position is written as a homozygous-alt SNP (GT 1/1), matching how the
haploid consensus represents a fixed difference, so the output can be scored directly
against the truth with compare_vcfs.py.

Both FASTAs must contain a single sequence of identical length; the script errors out
otherwise (a length mismatch would mean an indel crept in and coordinates are no longer
identity -- which must be handled by liftover, not this script).
"""

import sys
import argparse
import pysam


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--reference", required=True, help="reference FASTA (single contig)")
    p.add_argument("--consensus", required=True, help="iteration consensus FASTA (single contig)")
    p.add_argument("--sample-name", default="consensus")
    p.add_argument("--out", required=True, help="output VCF path (plain text)")
    return p.parse_args()


def load_single(path):
    fa = pysam.FastaFile(path)
    if len(fa.references) != 1:
        sys.exit("ERROR: %s must contain exactly one sequence (found %d)"
                 % (path, len(fa.references)))
    name = fa.references[0]
    return name, fa.fetch(name).upper()


def main():
    args = parse_args()
    ref_name, ref = load_single(args.reference)
    con_name, con = load_single(args.consensus)

    if len(ref) != len(con):
        sys.exit("ERROR: reference (%d) and consensus (%d) differ in length -- an indel "
                 "is present, so coordinates are NOT identity and this SNP-only diff is "
                 "invalid. Use liftover instead." % (len(ref), len(con)))

    bases = set("ACGT")
    with open(args.out, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write("##source=consensus_diff_snps.py (s01)\n")
        out.write("##contig=<ID=%s,length=%d>\n" % (ref_name, len(ref)))
        out.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n"
                  % args.sample_name)
        n = 0
        for i, (r, c) in enumerate(zip(ref, con)):
            # only call a SNP where both bases are ACGT and they differ
            if r != c and r in bases and c in bases:
                out.write("%s\t%d\t.\t%s\t%s\t.\tPASS\t.\tGT\t1/1\n"
                          % (ref_name, i + 1, r, c))
                n += 1
    sys.stderr.write("# %d SNP differences (consensus vs reference)\n" % n)


if __name__ == "__main__":
    main()
