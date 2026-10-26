#!/usr/bin/env python3
"""
backfill_combined_summaries.py -- rewrite the dictionaries of the combined summary tables (written by
analyses/lib/combine_summaries.py inside the section workflows) so they carry the FULL column
definitions and their producing command. Companion to backfill_dictionaries.sh.

The old dictionaries of these tables only said "see that script's *-summary.dict.txt", which lives in
the git-ignored data/ tree. For every committed dictionary still saying that, this script
  1. finds the table's entry in its section's MANIFEST.yaml (the producing command is taken from
     there, with the entry's {placeholders} filled from the table's path),
  2. infers the --kind from the table's header (compare_vcfs / compare_bams / allele_balance columns),
  3. runs combine_summaries.py --dicts-only, which rewrites the dictionary as a backfill.
No table is regenerated. Fails loudly if a table has no manifest entry or an unknown header.

Run from the project root:  conda run -n mapping-sims python analyses/lib/backfill/backfill_combined_summaries.py
"""
import os
import subprocess
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "analyses", "lib"))
import check_manifest as cm                    # expand() / fill() / manifest layout
import yaml

OLD_MARK = "see that script's"                 # the pointer line of the old dictionaries


def kind_of(header):
    """--kind for combine_summaries.py, from the table's columns."""
    cols = header.split(",")
    if "total_reads" in cols:
        return "readclass-bam"
    if "mean_raf" in cols:
        return "ab"
    if "div_tp" in cols:
        return "vcf"
    sys.exit("ERROR: cannot tell the kind of a table with header %s" % header)


def main():
    os.chdir(ROOT)
    tracked = subprocess.run(["git", "ls-files", "analyses/*.dict.txt"], check=True,
                             capture_output=True, text=True).stdout.split()
    targets = [d[:-len(".dict.txt")] for d in tracked if OLD_MARK in open(d).read()]

    # manifest entries by expanded path
    entries = {}
    for s in sorted(os.listdir("analyses")):
        mp = os.path.join("analyses", s, "MANIFEST.yaml")
        if not os.path.exists(mp):
            continue
        for e in yaml.safe_load(open(mp)).get("outputs") or []:
            if e.get("path"):
                for f in cm.expand(e["path"], e.get("values")):
                    entries[f] = e

    for table in targets:
        e = entries.get(table)
        if e is None:
            sys.exit("ERROR: %s has no MANIFEST.yaml entry (write the manifest first)" % table)
        command = cm.fill(e["command"], e["path"], table)
        header = open(table).readline().rstrip("\n")
        # a non-current manifest status (e.g. a table mixing two tool versions) is carried into the
        # dictionary as a provenance note, so the caveat travels with the table
        notes = []
        if e.get("status", "current") != "current" or e.get("note"):
            notes = ["--dict-note", "%s%s" % ("" if e.get("status", "current") == "current"
                                              else e["status"].upper() + ": ", e.get("note", ""))]
        subprocess.run([sys.executable, "analyses/lib/combine_summaries.py", "--out", table,
                        "--kind", kind_of(header), "--dicts-only", command] + notes, check=True)
        print("rewrote %s.dict.txt  (kind %s)" % (table, kind_of(header)))
    print("%d combined-summary dictionaries rewritten" % len(targets))


if __name__ == "__main__":
    main()
