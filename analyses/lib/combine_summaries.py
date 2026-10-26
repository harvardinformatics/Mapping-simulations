#!/usr/bin/env python3
"""
combine_summaries.py -- concatenate per-condition summary CSVs (one row each, from
compare_bams.py, compare_vcfs.py or allele_balance.py across coverage/divergence
combinations) into a single combined table for downstream analysis.

All inputs must share an identical header (they do, since they come from the same
generating script). The header is written once, followed by the data row(s) of every
input file in the order given.

A sidecar data dictionary is written with the FULL column definitions of the generating
script embedded (taken from that script, so the two never drift apart). This matters
because the combined tables are committed to git, while the per-condition dictionaries
live in the git-ignored data/ tree.

--dicts-only documents an existing combined table without re-reading the inputs (backfill).
"""

import sys
import argparse
import provenance                       # shared provenance block for the data dictionaries
import compare_vcfs                     # FIELD_DEFINITIONS / BASE_COLUMNS of -vcf-summary.csv
import compare_bams                     # FIELD_DEFINITIONS of -bam-summary.csv
import allele_balance                   # field_definitions() of -ab-summary.csv

# --kind -> (generating script, what one input row is)
KINDS = {"vcf": ("compare_vcfs.py", "variant-calling accuracy"),
         "bam": ("compare_bams.py", "read-mapping accuracy"),
         "readclass-bam": ("compare_bams.py", "read-mapping accuracy (read classification)"),
         "ab": ("allele_balance.py", "het-site allele balance")}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--inputs", nargs="+",
                   help="per-condition summary CSVs to concatenate (required unless --dicts-only)")
    p.add_argument("--out", required=True, help="combined output CSV path")
    p.add_argument("--kind", required=True, choices=sorted(KINDS),
                   help="which script made the inputs; selects the column definitions")
    p.add_argument("--dicts-only", metavar="ORIGINAL_COMMAND",
                   help="write ONLY the data dictionary for an existing --out table, recording "
                        "ORIGINAL_COMMAND as provenance (backfill)")
    p.add_argument("--dict-note", action="append", default=[],
                   help="extra provenance note for the dictionary (repeatable)")
    a = p.parse_args()
    if not a.dicts_only and not a.inputs:
        p.error("--inputs is required unless --dicts-only is given")
    return a


def definitions(kind, header):
    """Column-definition text for a combined table of this kind and header.

    The header is reconciled with the generating script's CURRENT columns:
      * columns the script writes but this table lacks are listed (the table was made by an
        earlier version of the script, before those columns existed);
      * vcf tables may end with compare_vcfs.py --extra condition columns, which are listed;
      * any other unknown column is an error -- the definitions would not describe the table."""
    script = KINDS[kind][0]
    cols = header.split(",")
    if kind == "vcf":
        known, text, extras_ok = compare_vcfs.BASE_COLUMNS, compare_vcfs.FIELD_DEFINITIONS, True
    elif kind in ("bam", "readclass-bam"):
        known, text, extras_ok = compare_bams.COLUMNS, compare_bams.FIELD_DEFINITIONS, False
    elif kind == "ab":
        known, text, extras_ok = allele_balance.COLUMNS, allele_balance.field_definitions(), False
    else:
        sys.exit("ERROR: unknown --kind %s" % kind)

    unknown = [c for c in cols if c not in known]
    if unknown and not extras_ok:
        sys.exit("ERROR: columns not written by %s: %s" % (script, ", ".join(unknown)))
    if unknown and cols[-len(unknown):] != unknown:
        sys.exit("ERROR: unknown columns are not trailing --extra columns: %s" % ", ".join(unknown))
    absent = [c for c in known if c not in cols]

    if absent:
        text += ("\nNOT IN THIS TABLE: %s -- defined above but absent here, because the table\n"
                 "was made by an earlier version of %s that did not write them yet.\n"
                 % (", ".join(absent), script))
    for name in unknown:
        text += "%-17s : extra condition column (compare_vcfs.py --extra %s=VALUE)\n" % (name, name)
    return text


def write_dictionary(a, header, n_inputs=None):
    script, what = KINDS[a.kind]
    if a.dicts_only:
        hdr = provenance.header(a.out, command=a.dicts_only, backfill=True, notes=a.dict_note)
    else:
        hdr = provenance.header(a.out, notes=a.dict_note)
    with open(a.out + ".dict.txt", "w") as d:
        d.write(hdr)
        d.write("Combined %s summary: the one-row-per-condition summary CSVs written by\n" % what)
        d.write("analyses/lib/%s, concatenated%s. One row per input (condition).\n\n"
                % (script, "" if n_inputs is None else " (%d inputs)" % n_inputs))
        d.write("Column definitions (from %s):\n\n" % script)
        d.write(definitions(a.kind, header))
        d.write("\nHeader: %s\n" % header)


def main():
    args = parse_args()

    # --dicts-only: document an existing combined table from its own header
    if args.dicts_only:
        with open(args.out) as fh:
            header = fh.readline().rstrip("\n")
        write_dictionary(args, header)
        return

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

    write_dictionary(args, header, n_inputs=len(args.inputs))

    print("combined %d input files -> %s (%d data rows)"
          % (len(args.inputs), args.out, len(rows)))


if __name__ == "__main__":
    main()
