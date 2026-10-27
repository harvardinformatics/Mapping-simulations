#!/usr/bin/env python3
"""
RETIRED 2026-10-26 -- kept for the record, not part of any current analysis.
The "phantom" sites this tool explains do not exist. They came from the first version of
effective_truth.py, which read bases through mason's -oa CIGAR; mason writes dense SNP clusters there as
fake 1D..1I gap pairs, so real alt sites looked all-reference. The corrected scan (analyses/lib/
effective_truth.py, reading bases at each read's haplotype position oP) finds no phantom sites. The
"~1 bp from another variant" result below therefore describes where the CIGAR artifact occurs (dense
clusters), not why mason drops variants, which it does not. Its outputs moved to
analyses/<section>/results/retired/phantom-adjacency.csv.

phantom_adjacency.py -- WHY are truth SNPs phantom (never simulated into the reads)?

mason_simulator with -iv only APPLIES the input VCF (it does not generate de novo variants, so nothing
mutates a site "back to reference"), and the truth VCF has no same-position duplicate records. Instead,
phantom sites are almost always IMMEDIATELY ADJACENT to another variant: when two variants land ~1 bp
apart, mason's haplotype materialization drops one of them, leaving that position as reference. This
tool quantifies that by comparing, for phantom vs non-phantom SNP sites, the distance to the nearest
OTHER truth variant. Emits one row per divergence (phantom vs non-phantom).

Usage (templates use {D}):
  python phantom_adjacency.py --truth "data/.../{D}/0.005/42/mm39-truth.norm.vcf.gz" \
      --phantom "analyses/.../effective-truth/eff-{D}-phantom.tsv" --out results/.../phantom-adjacency.csv \
      --divs 0.04,0.06,0.08,0.10
"""
import argparse
import bisect
import gzip
import random
import statistics

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))   # analyses/lib
import provenance                       # shared provenance block for the data dictionaries

# the phantom flags this tool reads are an artifact (see the module docstring)
ARTIFACT_NOTE = ("KNOWN ARTIFACT: the phantom flags this table is built from came from reading bases "
                 "through mason's -oa CIGAR (dense SNP clusters written as fake 1D..1I gap pairs); the "
                 "corrected scan (effective_truth.py, 2026-10-26) finds no phantom sites. Retired.")


def write_dictionary(a):
    """Dictionary for the --out table (backfill under the ORIGINAL command with --dicts-only)."""
    kw = dict(command=a.dicts_only, backfill=True) if a.dicts_only else {}
    with open(a.out + ".dict.txt", "w") as d:
        d.write(provenance.header(a.out, notes=[ARTIFACT_NOTE], **kw))
        d.write("Distance from each phantom-flagged hom-div SNP to the nearest OTHER truth variant (any\n"
                "type), vs a random sample of non-phantom hom-div SNPs. Phantom flags from %s;\n"
                "truth VCF %s. Up to %d sites per group (random.seed(1)).\n\n"
                "One row per (divergence, group).\n\nFields:\n"
                "  divergence      : simulated divergence D\n"
                "  group           : phantom | non_phantom\n"
                "  n               : sites measured\n"
                "  median_nn_bp    : median distance (bp) to the nearest other truth variant\n"
                "  mean_nn_bp      : mean of that distance\n"
                "  frac_within_2bp : fraction of sites with a truth neighbour <= 2 bp away\n"
                % (a.phantom, a.truth, a.sample))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", required=True, help="truth VCF template; {D}")
    ap.add_argument("--phantom", required=True, help="phantom TSV template (from effective_truth.py); {D}")
    ap.add_argument("--out", required=True)
    ap.add_argument("--divs", default="0.04,0.06,0.08,0.10")
    ap.add_argument("--sample", type=int, default=3000)
    ap.add_argument("--dicts-only", metavar="ORIGINAL_COMMAND",
                    help="write ONLY the data dictionary for an existing --out table, recording "
                         "ORIGINAL_COMMAND as provenance (backfill)")
    a = ap.parse_args()
    if a.dicts_only:
        write_dictionary(a)
        return
    random.seed(1)

    rows = []
    for D in a.divs.split(","):
        allpos = []                                     # every truth variant position (any type)
        homdiv = []                                     # hom-div biallelic SNP positions
        for line in gzip.open(a.truth.format(D=D), "rt"):
            if line[0] == "#":
                continue
            c = line.rstrip("\n").split("\t")
            p = int(c[1])
            allpos.append(p)
            if len(c[3]) == 1 and len(c[4]) == 1 and c[9].split(":")[0] in ("1|1", "1/1"):
                homdiv.append(p)
        allpos.sort()
        phantom = set()
        for line in open(a.phantom.format(D=D)):
            if line.startswith("pos"):
                continue
            phantom.add(int(line.split("\t")[0]))

        def nn(p):                                      # distance to nearest OTHER truth variant
            i = bisect.bisect_left(allpos, p)
            best = 10 ** 9
            for j in (i - 1, i + 1):
                if 0 <= j < len(allpos) and allpos[j] != p:
                    best = min(best, abs(allpos[j] - p))
            return best

        ph = [p for p in homdiv if p in phantom]
        nonph = [p for p in random.sample(homdiv, min(a.sample, len(homdiv))) if p not in phantom]
        for label, pts in (("phantom", ph[:a.sample]), ("non_phantom", nonph)):
            d = [nn(p) for p in pts]
            if not d:
                continue
            rows.append((D, label, len(d), statistics.median(d), round(statistics.mean(d), 2),
                         round(sum(x <= 2 for x in d) / len(d), 3)))

    with open(a.out, "w") as out:
        out.write("divergence,group,n,median_nn_bp,mean_nn_bp,frac_within_2bp\n")
        for r in rows:
            out.write(",".join(str(x) for x in r) + "\n")
    write_dictionary(a)
    print(f"wrote {a.out}")
    print(f"{'D':>5} {'group':>12} {'median_nn':>10} {'frac<=2bp':>10}")
    for D, g, n, md, mn, fr in rows:
        print(f"{D:>5} {g:>12} {md:>10} {fr:>10}")


if __name__ == "__main__":
    main()
