#!/usr/bin/env python3
"""
allele_balance.py -- reference-allele bias at heterozygous sites.

The canonical reference-bias measurement: at a true het site (0/1) the alt-bearing reads
map worse than the ref-bearing reads, so the allele balance skews toward the reference
(reference-allele fraction RAF = ref_reads/(ref_reads+alt_reads) rises above the true 0.5),
and the site drops out or gets miscalled. This is a MAPPING-side measurement taken straight
from the pileup, independent of the variant caller.

For each het SNP in the truth VCF, we count reads supporting the ref vs alt base in the
mapped BAM (via pysam count_coverage, MAPQ- and base-quality-filtered). Iteration 1 maps to
the reference (het positions are reference coordinates). For a length-shifted consensus
(s02 iter>1) pass --chain (the ref->consensus chain) to lift each het position into the
BAM's coordinate system; the ref/alt bases are unchanged because het sites are never baked
into the consensus.

Output: <out_prefix>-ab-summary.csv -- one row of summary stats + a 10-bin RAF histogram.
        <out_prefix>-ab-summary.dict.txt -- data dictionary.
"""

import sys
import gzip
import bisect
import argparse
import pysam
import provenance                       # shared provenance block for the data dictionaries

BASE_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}
NBINS = 10
# columns of the -ab-summary.csv table
COLUMNS = (("region", "coverage", "divergence", "heterozygosity", "seed", "caller", "iteration",
            "n_het", "n_covered", "n_nocov", "lift_fail", "mean_depth",
            "mean_raf", "frac_ref_biased", "frac_alt_dropout")
           + tuple("raf_b%d" % i for i in range(NBINS)))


def field_definitions(min_mapq=None, min_baseq=None, dropout_alt_frac=None):
    """Column definitions of the -ab-summary.csv table. A function (not a constant) because the
    thresholds are run arguments; combine_summaries.py calls it without them, for the combined
    (committed) tables, and the text then names the arguments instead of their values."""
    mapq = "%d" % min_mapq if min_mapq is not None else "--min-mapq (default 20)"
    baseq = "%d" % min_baseq if min_baseq is not None else "--min-baseq (default 13)"
    drop = ("%.2f" % dropout_alt_frac if dropout_alt_frac is not None
            else "--dropout-alt-frac (default 0.2)")
    return ("Reference-allele bias at true heterozygous SNP sites, from the pileup.\n"
            "RAF = ref_reads/(ref_reads+alt_reads); truth = 0.5, reference bias => RAF > 0.5.\n"
            "MAPQ >= %s, base quality >= %s; iter>1 positions lifted via --chain.\n\n" % (mapq, baseq) +
            "region,coverage,divergence,heterozygosity,seed : simulation condition\n"
            "caller          : variant caller whose consensus this iteration mapped to\n"
            "iteration       : mapping iteration (1 = single-shot map to the reference)\n"
            "n_het           : het SNP sites in the truth\n"
            "n_covered       : het sites with >=1 informative (ref or alt) read (RAF computed)\n"
            "n_nocov         : het sites with 0 informative reads (fully uncovered)\n"
            "lift_fail       : het positions that did not lift (chain gap); 0 without --chain\n"
            "mean_depth      : mean informative depth at covered het sites\n"
            "mean_raf        : mean reference-allele fraction (0.5 = unbiased, >0.5 = ref bias)\n"
            "frac_ref_biased : fraction of covered het sites with RAF > 0.5\n"
            "frac_alt_dropout: fraction of covered het sites with alt fraction < %s (near-lost)\n" % drop +
            "raf_b0..raf_b9  : histogram of RAF over covered sites, bins [0,0.1),...,[0.9,1.0]\n")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--truth-vcf", required=True, help="truth VCF (mason format ok)")
    p.add_argument("--bam", required=True, help="mapped BAM for this iteration")
    p.add_argument("--chain", default=None,
                   help="optional ref->consensus chain to lift het positions into the BAM's "
                        "coordinates (s02 iter>1); omit for identity coords (s01, or iter1)")
    p.add_argument("--region", required=True)
    p.add_argument("--coverage", required=True)
    p.add_argument("--divergence", required=True)
    p.add_argument("--heterozygosity", required=True)
    p.add_argument("--seed", required=True)
    p.add_argument("--caller", required=True)
    p.add_argument("--iteration", required=True)
    p.add_argument("--min-mapq", type=int, default=20)
    p.add_argument("--min-baseq", type=int, default=13)
    p.add_argument("--dropout-alt-frac", type=float, default=0.2,
                   help="a covered het site with alt fraction below this is 'dropped'")
    p.add_argument("--out-prefix", required=True)
    return p.parse_args()


def load_chain(path):
    """ref->query(consensus) point lift from a UCSC chain (single contig, same strand)."""
    starts, ends, offs = [], [], []
    t = q = 0
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith("chain"):
                f = line.split(); t, q = int(f[5]), int(f[10]); continue
            p = line.split(); size = int(p[0])
            starts.append(t); ends.append(t + size); offs.append(q - t)
            if len(p) == 3:
                t += size + int(p[1]); q += size + int(p[2])

    def lift(x):
        i = bisect.bisect_right(starts, x) - 1
        return x + offs[i] if i >= 0 and x < ends[i] else None
    return lift


def load_het_snps(path):
    """Return [(pos0, ref, alt)] for heterozygous SNP truth sites (GT 0/1 or 1/0)."""
    sites = []
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 10:
                continue
            ref, alt, gt = c[3], c[4], c[9]
            gt = gt.split(":")[0]
            if len(ref) == 1 and len(alt) == 1 and gt in ("0|1", "1|0", "0/1", "1/0"):
                sites.append((int(c[1]) - 1, ref.upper(), alt.upper()))
    return sites


def main():
    args = parse_args()
    lift = load_chain(args.chain) if args.chain else None
    het = load_het_snps(args.truth_vcf)
    sys.stderr.write("# %d het SNP truth sites\n" % len(het))

    bam = pysam.AlignmentFile(args.bam, "rb")
    clen = bam.get_reference_length(args.region)

    def keep(r):
        return (not r.is_secondary and not r.is_supplementary and not r.is_duplicate
                and not r.is_unmapped and r.mapping_quality >= args.min_mapq)

    # per-base coverage over the whole contig (A,C,G,T), base-quality filtered
    cov = bam.count_coverage(args.region, 0, clen,
                             quality_threshold=args.min_baseq, read_callback=keep)

    n_het = len(het)
    n_cov = n_nocov = 0
    sum_raf = 0.0
    sum_depth = 0
    ref_biased = 0          # RAF > 0.5 (skewed toward reference)
    dropout = 0             # alt fraction below --dropout-alt-frac
    lift_fail = 0
    bins = [0] * NBINS
    for pos0, ref, alt in het:
        p = lift(pos0) if lift is not None else pos0
        if p is None:
            lift_fail += 1
            continue
        if p < 0 or p >= clen or ref not in BASE_IDX or alt not in BASE_IDX:
            continue
        rc = cov[BASE_IDX[ref]][p]
        ac = cov[BASE_IDX[alt]][p]
        d = rc + ac
        if d == 0:
            n_nocov += 1
            continue
        n_cov += 1
        raf = rc / d
        sum_raf += raf
        sum_depth += d
        bins[min(int(raf * NBINS), NBINS - 1)] += 1
        if raf > 0.5:
            ref_biased += 1
        if (ac / d) < args.dropout_alt_frac:
            dropout += 1

    mean_raf = sum_raf / n_cov if n_cov else 0.0
    mean_depth = sum_depth / n_cov if n_cov else 0.0
    frac_ref_biased = ref_biased / n_cov if n_cov else 0.0
    frac_dropout = dropout / n_cov if n_cov else 0.0

    head = list(COLUMNS)
    row = ([args.region, args.coverage, args.divergence, args.heterozygosity, args.seed,
            args.caller, args.iteration, n_het, n_cov, n_nocov, lift_fail,
            "%.2f" % mean_depth, "%.6f" % mean_raf, "%.6f" % frac_ref_biased, "%.6f" % frac_dropout]
           + [bins[i] for i in range(NBINS)])

    csv_path = args.out_prefix + "-ab-summary.csv"
    with open(csv_path, "w") as out:
        out.write(",".join(head) + "\n")
        out.write(",".join(str(x) for x in row) + "\n")

    with open(args.out_prefix + "-ab-summary.dict.txt", "w") as d:
        d.write(provenance.header(csv_path))
        d.write(field_definitions(args.min_mapq, args.min_baseq, args.dropout_alt_frac))

    print("caller=%s iter=%s div=%s het=%s: n_covered=%d mean_RAF=%.4f frac_ref_biased=%.3f "
          "frac_dropout=%.3f mean_depth=%.1f" %
          (args.caller, args.iteration, args.divergence, args.heterozygosity,
           n_cov, mean_raf, frac_ref_biased, frac_dropout, mean_depth))


if __name__ == "__main__":
    main()
