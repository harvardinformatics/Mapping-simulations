#!/usr/bin/env python3
"""
Position-accuracy comparison for the mason2 POC (Phase 1) -- previews the v2
`compare_bams` logic.

For every simulated read we compare:
  * TRUTH position  : where the read really came from, taken from the mason
                      `-oa` alignment BAM (already in REFERENCE coordinates -- no
                      CrossMap / liftover needed).
  * MAPPED position : where `bwa mem` actually placed the read.

Because the POC reference is a single contig (chr19), a read falls into one of:
  exact       - bwa placed it at the exact true reference position
  close       - within TOLERANCE bp of the true position (one read length)
  mismapped   - same chromosome but farther than TOLERANCE from truth
  diff_chr    - different chromosome (impossible here; kept for generality)
  unmapped    - bwa could not place the read
  missing     - read present in truth but absent from the mapped BAM
                (this is the v1 CrossMap artifact; expected to be 0 here)

Reads are matched by (fragment id, mate). mason names reads 'simulated.<N>', so
we encode (N, mate) as a single integer key to keep the truth lookup table small.

Only PRIMARY alignments are counted (secondary/supplementary are skipped).
Positions use pysam `reference_start`, which is 0-based for both BAMs and thus
directly comparable.

Outputs:
  <out_prefix>-summary.csv        : one row of class counts + proportions
  <out_prefix>-summary.dict.txt   : data dictionary for the CSV columns
"""

import sys
import argparse
import pysam
from collections import Counter

# --- parameters -------------------------------------------------------------

TOLERANCE = 150  # bp; a read landing within one read-length of truth is "close"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--truth-bam", required=True, help="mason -oa truth alignment BAM")
    p.add_argument("--mapped-bam", required=True, help="bwa mem mapped BAM")
    p.add_argument("--region", default="19", help="region label for the summary row")
    p.add_argument("--coverage", default="30", help="coverage label for the summary row")
    p.add_argument("--divergence", default="0.04", help="divergence label for the summary row")
    p.add_argument("--out-prefix", required=True, help="output path prefix")
    return p.parse_args()


def read_key(aln):
    """Encode (fragment number, mate) as one integer.

    mason read names look like 'simulated.<N>'. Read 1 -> even key, read 2 -> odd.
    """
    n = int(aln.query_name.rsplit(".", 1)[1])
    return n * 2 + (0 if aln.is_read1 else 1)


def main():
    args = parse_args()

    # --- 1. load truth positions: key -> (tid, pos) -------------------------
    # Single pass over the (smaller, all-mapped) truth BAM builds the lookup table.
    truth = {}
    with pysam.AlignmentFile(args.truth_bam, "rb") as tf:
        for a in tf:
            if a.is_secondary or a.is_supplementary:
                continue
            truth[read_key(a)] = (a.reference_id, a.reference_start)
    sys.stderr.write("# loaded %d truth reads\n" % len(truth))

    # --- 2. stream mapped BAM and classify each read ------------------------
    counts = Counter()
    with pysam.AlignmentFile(args.mapped_bam, "rb") as mf:
        for a in mf:
            if a.is_secondary or a.is_supplementary:
                continue
            # pop so that whatever remains in `truth` at the end is 'missing'
            t = truth.pop(read_key(a), None)
            if t is None:
                counts["not_in_truth"] += 1
                continue
            if a.is_unmapped:
                counts["unmapped"] += 1
                continue
            truth_tid, truth_pos = t
            if a.reference_id != truth_tid:
                counts["diff_chr"] += 1
                continue
            delta = abs(a.reference_start - truth_pos)
            if delta == 0:
                counts["exact"] += 1
            elif delta <= TOLERANCE:
                counts["close"] += 1
            else:
                counts["mismapped"] += 1

    # reads still in truth were never seen in the mapped BAM
    counts["missing"] = len(truth)

    # --- 3. write summary CSV ----------------------------------------------
    classes = ["exact", "close", "mismapped", "diff_chr", "unmapped", "missing"]
    total = sum(counts[c] for c in classes)

    header = (["region", "coverage", "divergence", "total_reads"]
              + classes
              + [c + "_prop" for c in classes])
    row = ([args.region, args.coverage, args.divergence, total]
           + [counts[c] for c in classes]
           + ["%.6f" % (counts[c] / total) if total else "0" for c in classes])

    csv_path = args.out_prefix + "-summary.csv"
    with open(csv_path, "w") as out:
        out.write(",".join(header) + "\n")
        out.write(",".join(str(x) for x in row) + "\n")

    # also report anything that showed up in mapped but not truth (should be 0)
    sys.stderr.write("# not_in_truth (mapped reads with no truth match): %d\n"
                     % counts["not_in_truth"])

    # --- 4. write the data dictionary --------------------------------------
    dict_path = args.out_prefix + "-summary.dict.txt"
    with open(dict_path, "w") as d:
        d.write("Data dictionary for %s\n" % csv_path)
        d.write("Position-accuracy of bwa-mapped reads vs mason truth (ref coords).\n")
        d.write("Read match = same (fragment id, mate); tolerance = %d bp.\n\n" % TOLERANCE)
        d.write("region       : reference region / chromosome label\n")
        d.write("coverage     : simulated fold coverage\n")
        d.write("divergence   : simulated sample-vs-reference divergence (target)\n")
        d.write("total_reads  : primary reads compared (sum of the six classes)\n")
        d.write("exact        : mapped to the exact true reference position (delta = 0)\n")
        d.write("close        : mapped within %d bp of the true position\n" % TOLERANCE)
        d.write("mismapped    : same chromosome but > %d bp from the true position\n" % TOLERANCE)
        d.write("diff_chr     : mapped to a different chromosome than truth\n")
        d.write("unmapped     : read present in truth but unmapped by bwa\n")
        d.write("missing      : read present in truth but absent from the mapped BAM\n")
        d.write("*_prop       : the above counts as a proportion of total_reads\n")

    # --- 5. echo the result to the log -------------------------------------
    print("region=%s coverage=%s divergence=%s total_reads=%d"
          % (args.region, args.coverage, args.divergence, total))
    for c in classes:
        print("  %-10s %12d  (%.4f)" % (c, counts[c], counts[c] / total if total else 0))
    print("wrote %s" % csv_path)
    print("wrote %s" % dict_path)


if __name__ == "__main__":
    main()
