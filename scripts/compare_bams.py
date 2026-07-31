#!/usr/bin/env python3
"""
compare_bams.py -- read-mapping accuracy for the v2 mason reference-bias pipeline.

Compares where reads were SIMULATED FROM (mason truth BAM, already in reference
coordinates) against where `bwa mem` actually PLACED them (mapped BAM). Because the
mason `-oa` truth alignment is in reference coordinates, NO CrossMap / liftover is
needed -- this is the key simplification over the v1 pipeline (whose liftover step
produced a constant, spurious "missing reads" artifact even at zero divergence).

Reads are matched by (fragment id, mate). mason names reads 'simulated.<N>', so we
encode (N, mate) as a single integer key to keep the truth lookup table compact.

Each PRIMARY read is placed in exactly one class, comparing by reference *name*
(so a truth BAM and a mapped BAM with different headers still compare correctly):

  exact       - bwa placed it at the exact true reference position (delta = 0)
  close       - within TOLERANCE bp of the true position (one read length)
  mismapped   - same reference contig but farther than TOLERANCE from truth
  diff_chr    - mapped to a different contig than truth
  unmapped    - present in truth but flagged unmapped by bwa
  missing     - present in truth but absent from the mapped BAM (v1 liftover artifact;
                expected to be 0 with mason truth)

Only primary alignments are counted (secondary/supplementary skipped). Positions use
pysam `reference_start` (0-based in both BAMs), so they are directly comparable.

Outputs:
  <out_prefix>-bam-summary.csv       : one row of class counts + proportions
  <out_prefix>-bam-summary.dict.txt  : data dictionary for that CSV
"""

import sys
import argparse
import pysam


# --- parameters -------------------------------------------------------------
TOLERANCE = 150  # bp; within one read length of truth counts as "close"


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--truth-bam", required=True, help="mason -oa truth alignment BAM")
    p.add_argument("--mapped-bam", required=True, help="bwa mem mapped BAM")
    p.add_argument("--region", required=True, help="region / chromosome label")
    p.add_argument("--coverage", required=True, help="coverage label")
    p.add_argument("--divergence", required=True, help="divergence label")
    p.add_argument("--heterozygosity", required=True, help="heterozygosity label")
    p.add_argument("--seed", required=True, help="replicate seed label")
    p.add_argument("--out-prefix", required=True, help="output path prefix")
    return p.parse_args()


def read_key(aln):
    """Encode (fragment number, mate) as one integer key.

    mason read names look like 'simulated.<N>'. Read 1 -> even, read 2 -> odd.
    """
    n = int(aln.query_name.rsplit(".", 1)[1])
    return n * 2 + (0 if aln.is_read1 else 1)


def main():
    args = parse_args()

    # --- 1. load truth: key -> (reference_name, reference_start) -----------
    # Compare by contig NAME (not tid) so differing BAM headers still line up.
    truth = {}
    with pysam.AlignmentFile(args.truth_bam, "rb") as tf:
        for a in tf:
            if a.is_secondary or a.is_supplementary:
                continue
            truth[read_key(a)] = (a.reference_name, a.reference_start)
    sys.stderr.write("# loaded %d truth reads\n" % len(truth))

    # --- 2. stream mapped BAM and classify --------------------------------
    # "close" is split into 'shifted' (whole read moved <=TOLERANCE bp, fully aligned)
    # vs 'clipped' (soft-clipped near truth). A soft-clip is the mechanistically
    # important case -- clipped bases (which carry the divergent alleles) drop out of
    # the pileup -- and it is tracked separately from a simple positional shift.
    counts = {k: 0 for k in
              ["exact", "shifted", "clipped", "mismapped", "diff_chr",
               "unmapped", "missing", "not_in_truth"]}
    with pysam.AlignmentFile(args.mapped_bam, "rb") as mf:
        for a in mf:
            if a.is_secondary or a.is_supplementary:
                continue
            t = truth.pop(read_key(a), None)   # pop => leftovers are 'missing'
            if t is None:
                counts["not_in_truth"] += 1
                continue
            if a.is_unmapped:
                counts["unmapped"] += 1
                continue
            truth_ref, truth_pos = t
            if a.reference_name != truth_ref:
                counts["diff_chr"] += 1
                continue
            delta = abs(a.reference_start - truth_pos)
            has_clip = a.cigartuples is not None and any(op == 4 for op, _ in a.cigartuples)
            if delta > TOLERANCE:                 # wrong location dominates
                counts["mismapped"] += 1
            elif has_clip:                        # correct locus but soft-clipped
                counts["clipped"] += 1
            elif delta == 0:                      # perfect, fully aligned
                counts["exact"] += 1
            else:                                 # whole read shifted a little, no clip
                counts["shifted"] += 1
    counts["missing"] = len(truth)  # never seen in mapped BAM

    # --- 3. write summary CSV ---------------------------------------------
    classes = ["exact", "shifted", "clipped", "mismapped", "diff_chr", "unmapped", "missing"]
    total = sum(counts[c] for c in classes)

    header = (["region", "coverage", "divergence", "heterozygosity", "seed", "total_reads"]
              + classes + [c + "_prop" for c in classes])
    row = ([args.region, args.coverage, args.divergence, args.heterozygosity, args.seed, total]
           + [counts[c] for c in classes]
           + ["%.6f" % (counts[c] / total) if total else "0" for c in classes])

    csv_path = args.out_prefix + "-bam-summary.csv"
    with open(csv_path, "w") as out:
        out.write(",".join(header) + "\n")
        out.write(",".join(str(x) for x in row) + "\n")

    sys.stderr.write("# not_in_truth (mapped reads with no truth match): %d\n"
                     % counts["not_in_truth"])

    # --- 4. data dictionary -----------------------------------------------
    with open(args.out_prefix + "-bam-summary.dict.txt", "w") as d:
        d.write("Data dictionary for %s\n" % csv_path)
        d.write("Read-mapping accuracy: bwa placement vs mason truth (reference coords).\n")
        d.write("Reads matched by (fragment id, mate); tolerance = %d bp; contigs compared by name.\n\n" % TOLERANCE)
        d.write("region         : reference region / chromosome label\n")
        d.write("coverage       : simulated fold coverage\n")
        d.write("divergence     : simulated fixed divergence (D) label\n")
        d.write("heterozygosity : simulated heterozygosity (H) label\n")
        d.write("seed           : replicate random seed\n")
        d.write("total_reads    : primary reads compared (sum of the placement classes)\n")
        d.write("exact          : exact true position (delta = 0), fully aligned, NO soft-clip\n")
        d.write("shifted        : correct locus, whole read moved 1-%d bp, no soft-clip\n" % TOLERANCE)
        d.write("clipped        : correct locus (<=%d bp) but soft-clipped (alt-bearing ends dropped);\n" % TOLERANCE)
        d.write("                 includes 3'-clipped reads that sit at the exact start position\n")
        d.write("mismapped      : same contig but > %d bp from the true position\n" % TOLERANCE)
        d.write("diff_chr       : mapped to a different contig than truth\n")
        d.write("unmapped       : present in truth but unmapped by bwa\n")
        d.write("missing        : present in truth but absent from the mapped BAM (v1 liftover artifact)\n")
        d.write("*_prop         : each count as a proportion of total_reads\n")

    # --- 5. echo to log ---------------------------------------------------
    print("region=%s cov=%s div=%s het=%s total=%d"
          % (args.region, args.coverage, args.divergence, args.heterozygosity, total))
    for c in classes:
        print("  %-10s %12d  (%.4f)" % (c, counts[c], counts[c] / total if total else 0))
    print("wrote %s" % csv_path)


if __name__ == "__main__":
    main()
