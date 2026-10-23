#!/usr/bin/env python3
"""
compare_bams.py -- read-mapping accuracy for the mason reference-bias pipeline.

Compares where reads were SIMULATED FROM (mason truth BAM, already in reference
coordinates) against where `bwa mem` actually PLACED them (mapped BAM). Because the
mason `-oa` truth alignment is in reference coordinates, NO CrossMap / liftover is
needed.

This is the analyses/lib copy used by s01's read-classification-by-iteration step.
It adds two columns over the exploratory scripts/compare_bams.py version -- `caller`
and `iteration` -- so per-iteration classifications from all three callers can be
concatenated into one table. It is otherwise identical.

Reads are matched by (fragment id, mate). mason names reads 'simulated.<N>', so we
encode (N, mate) as a single integer key to keep the truth lookup table compact.

Each PRIMARY read is placed in exactly one class, comparing by reference *name*
(so a truth BAM in reference coordinates and a BAM mapped to an equal-length,
identity-coordinate consensus still compare correctly at every iteration):

  exact       - bwa placed it at the exact true reference position (delta = 0)
  shifted     - correct locus, whole read moved 1..TOLERANCE bp, no soft-clip
  clipped     - correct locus (<=TOLERANCE bp) but soft-clipped (alt-bearing ends drop)
  mismapped   - same reference contig but farther than TOLERANCE from truth
  diff_chr    - mapped to a different contig than truth
  unmapped    - present in truth but flagged unmapped by bwa
  missing     - present in truth but absent from the mapped BAM (expected 0 with mason)

Only primary alignments are counted (secondary/supplementary skipped). Positions use
pysam `reference_start` (0-based in both BAMs), so they are directly comparable.

Outputs:
  <out_prefix>-bam-summary.csv       : one row of class counts + proportions
  <out_prefix>-bam-summary.dict.txt  : data dictionary for that CSV
"""

import sys
import bisect
import argparse
import pysam


# --- parameters -------------------------------------------------------------
TOLERANCE = 150  # bp; within one read length of truth counts as "close"
# We also point-lift the read END and check the reference footprint: a read placed at
# the right START but spanning the WRONG indels (wrong footprint) is 'span_off', not
# exact. END_TOLERANCE absorbs the few-bp indel-boundary jitter (measured p95 ~5 bp).
# For SNP-only (identity coords, no --chain) start and end shift together, so span_off
# is ~0 and behavior matches the start-only version.
END_TOLERANCE = 10


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
    p.add_argument("--caller", required=True, help="variant caller whose iteration this BAM is from")
    p.add_argument("--iteration", required=True, help="pseudo-it iteration (1 = single-shot)")
    p.add_argument("--chain", default=None,
                   help="optional UCSC chain (target=mapped-BAM coords, query=reference); "
                        "point-lifts each mapped read's start into reference coords. Use for "
                        "iterations whose consensus is length-shifted by indels; omit for "
                        "identity-coordinate (SNP-only) iterations.")
    p.add_argument("--out-prefix", required=True, help="output path prefix")
    return p.parse_args()


def read_key(aln):
    """Encode (fragment number, mate) as one integer key.

    mason read names look like 'simulated.<N>'. Read 1 -> even, read 2 -> odd.
    """
    n = int(aln.query_name.rsplit(".", 1)[1])
    return n * 2 + (0 if aln.is_read1 else 1)


def make_lifter(chain_path):
    """Point-liftover from a UCSC chain (single contig, same strand -- as emitted by
    `bcftools consensus -c` then chainSwap'd to target=consensus, query=reference).
    Returns lift(pos) -> pos-in-reference, or None if pos falls in a chain gap (the read
    starts inside a deletion). Only the read START is lifted; unlike a full alignment
    re-projection (CrossMap bam), a point lift does NOT drop indel-spanning reads --
    essential at high indel density, where most reads cross an indel."""
    starts, ends, offs = [], [], []
    t = q = 0
    with open(chain_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith("chain"):
                f = line.split()
                t, q = int(f[5]), int(f[10])          # tStart (source), qStart (dest)
                continue
            p = line.split()
            size = int(p[0])
            starts.append(t); ends.append(t + size); offs.append(q - t)
            if len(p) == 3:                            # gap sizes to the next block
                t += size + int(p[1])
                q += size + int(p[2])

    def lift(pos):
        i = bisect.bisect_right(starts, pos) - 1
        if i >= 0 and pos < ends[i]:
            return pos + offs[i]
        return None
    return lift


def main():
    args = parse_args()

    # --- 1. load truth: key -> (reference_name, reference_start) -----------
    # Compare by contig NAME (not tid) so differing BAM headers still line up.
    truth = {}
    with pysam.AlignmentFile(args.truth_bam, "rb") as tf:
        for a in tf:
            if a.is_secondary or a.is_supplementary:
                continue
            truth[read_key(a)] = (a.reference_name, a.reference_start, a.reference_end)
    sys.stderr.write("# loaded %d truth reads\n" % len(truth))
    lift = make_lifter(args.chain) if args.chain else None

    # --- 2. stream mapped BAM and classify --------------------------------
    # "close" is split into 'shifted' (whole read moved <=TOLERANCE bp, fully aligned)
    # vs 'clipped' (soft-clipped near truth). A soft-clip is the mechanistically
    # important case -- clipped bases (which carry the divergent alleles) drop out of
    # the pileup -- and it is tracked separately from a simple positional shift.
    counts = {k: 0 for k in
              ["exact", "shifted", "clipped", "span_off", "mismapped", "diff_chr",
               "unmapped", "missing", "not_in_truth", "lift_fail"]}
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
            truth_ref, truth_start, truth_end = t
            if a.reference_name != truth_ref:
                counts["diff_chr"] += 1
                continue
            mapped_start = a.reference_start
            end_lifted   = a.reference_end - 1          # last aligned ref base (inclusive)
            if lift is not None:                        # consensus coords -> reference coords
                mapped_start = lift(mapped_start)
                if mapped_start is None:                # START in a chain gap -> can't place at all
                    counts["lift_fail"] += 1
                    continue
                end_lifted = lift(end_lifted)           # may be None: read ENDS inside an insertion
            start_delta = abs(mapped_start - truth_start)
            has_clip = a.cigartuples is not None and any(op == 4 for op, _ in a.cigartuples)
            if start_delta > TOLERANCE:                 # wrong location dominates
                counts["mismapped"] += 1
            elif has_clip:                              # correct locus but soft-clipped
                counts["clipped"] += 1
            elif start_delta != 0:                      # whole read start shifted a little, no clip
                counts["shifted"] += 1
            elif end_lifted is not None and abs(end_lifted - (truth_end - 1)) > END_TOLERANCE:
                counts["span_off"] += 1                 # right start, wrong reference footprint
            else:                                       # start exact + footprint ok (or end in an insertion)
                counts["exact"] += 1
    counts["missing"] = len(truth)  # never seen in mapped BAM

    # --- 3. write summary CSV ---------------------------------------------
    classes = ["exact", "shifted", "clipped", "span_off", "mismapped", "diff_chr", "unmapped", "missing"]
    total = sum(counts[c] for c in classes)

    header = (["region", "coverage", "divergence", "heterozygosity", "seed",
               "caller", "iteration", "total_reads"]
              + classes + [c + "_prop" for c in classes])
    row = ([args.region, args.coverage, args.divergence, args.heterozygosity, args.seed,
            args.caller, args.iteration, total]
           + [counts[c] for c in classes]
           + ["%.6f" % (counts[c] / total) if total else "0" for c in classes])

    csv_path = args.out_prefix + "-bam-summary.csv"
    with open(csv_path, "w") as out:
        out.write(",".join(header) + "\n")
        out.write(",".join(str(x) for x in row) + "\n")

    sys.stderr.write("# not_in_truth (mapped reads with no truth match): %d\n"
                     % counts["not_in_truth"])
    sys.stderr.write("# lift_fail (start in a chain gap, point-lift): %d\n"
                     % counts["lift_fail"])

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
        d.write("caller         : variant caller whose consensus this iteration mapped to\n")
        d.write("iteration      : pseudo-it iteration (1 = single-shot map to the reference)\n")
        d.write("total_reads    : primary reads compared (sum of the placement classes)\n")
        d.write("exact          : true START (delta = 0) AND reference footprint matches truth\n")
        d.write("                 (lifted END within %d bp), fully aligned, NO soft-clip\n" % END_TOLERANCE)
        d.write("shifted        : correct locus, whole read START moved 1-%d bp, no soft-clip\n" % TOLERANCE)
        d.write("clipped        : correct locus (<=%d bp) but soft-clipped (alt-bearing ends dropped);\n" % TOLERANCE)
        d.write("                 includes 3'-clipped reads that sit at the exact start position\n")
        d.write("span_off       : correct START but reference footprint wrong by > %d bp -- the read\n" % END_TOLERANCE)
        d.write("                 spans the wrong indel structure (only distinguishable with --chain)\n")
        d.write("mismapped      : same contig but > %d bp from the true position\n" % TOLERANCE)
        d.write("diff_chr       : mapped to a different contig than truth\n")
        d.write("unmapped       : present in truth but unmapped by bwa\n")
        d.write("missing        : present in truth but absent from the mapped BAM (expected 0)\n")
        d.write("*_prop         : each count as a proportion of total_reads\n")

    # --- 5. echo to log ---------------------------------------------------
    print("caller=%s iter=%s region=%s cov=%s div=%s het=%s total=%d"
          % (args.caller, args.iteration, args.region, args.coverage,
             args.divergence, args.heterozygosity, total))
    for c in classes:
        print("  %-10s %12d  (%.4f)" % (c, counts[c], counts[c] / total if total else 0))
    print("wrote %s" % csv_path)


if __name__ == "__main__":
    main()
