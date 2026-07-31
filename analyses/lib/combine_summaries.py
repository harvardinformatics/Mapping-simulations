#!/usr/bin/env python3
"""
combine_summaries.py -- concatenate per-condition summary CSVs (one row each, from
compare_bams.py or compare_vcfs.py across coverage/divergence combinations) into a
single combined table for downstream analysis.

All inputs must share an identical header (they do, since they come from the same
generating script). The header is written once, followed by the data row(s) of every
input file in the order given. A sidecar data dictionary is written that points back
to the per-file dictionary describing the columns.
"""

import sys
import argparse


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--inputs", nargs="+", required=True,
                   help="per-condition summary CSVs to concatenate")
    p.add_argument("--out", required=True, help="combined output CSV path")
    p.add_argument("--kind", default="summary",
                   help="label for the data dictionary (e.g. 'bam' or 'vcf')")
    return p.parse_args()


def main():
    args = parse_args()

    header = None
    rows = []
    for path in args.inputs:
        with open(path) as fh:
            lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
        if not lines:
            sys.exit("ERROR: empty summary file: %s" % path)
        this_header, this_rows = lines[0], lines[1:]
        if header is None:
            header = this_header
        elif this_header != header:
            sys.exit("ERROR: header mismatch in %s\n  expected: %s\n  found:    %s"
                     % (path, header, this_header))
        rows.extend(this_rows)

    with open(args.out, "w") as out:
        out.write(header + "\n")
        for r in rows:
            out.write(r + "\n")

    with open(args.out + ".dict.txt", "w") as d:
        d.write("Data dictionary for %s\n" % args.out)
        d.write("Combined %s summary: one row per (coverage, divergence) condition.\n" % args.kind)
        d.write("Columns are identical to the per-condition files produced by the\n")
        d.write("compare_%ss.py script; see that script's *-summary.dict.txt for full\n" % args.kind)
        d.write("column definitions.\n")
        d.write("Header: %s\n" % header)

    print("combined %d input files -> %s (%d data rows)"
          % (len(args.inputs), args.out, len(rows)))


if __name__ == "__main__":
    main()
