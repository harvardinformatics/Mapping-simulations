#!/usr/bin/env python3
"""
mapping_vs_calling.py -- the s05 mechanism table: split every UNRECOVERED divergent
SNP into a MAPPING failure vs a CALLING failure, as a function of local divergence.

For each divergent truth site (from variant_read_fate.py --pileup output) we know, per
caller x iteration:
  * recovered (0/1),
  * alt_p        = alt-supporting reads in the caller's pileup AT the site (evidence
                   actually present for the caller to see), and
  * the truth-centric read fate (ex_alt/sh_alt = alt-carrying reads that arrived;
    ex_clip/sh_clip = arrived but soft-clipped over p; mm_*/dc/un = mismapped away).

Classification of each site:
  * recovered        -- called (in the consensus diff).
  * calling-limited  -- NOT recovered, but the alt allele IS present in the pileup with
                        >= --alt-min supporting reads: the evidence was there and the
                        caller did not emit it (assembler/genotyper failure).
  * mapping-limited  -- NOT recovered, and the alt allele is absent/too weak at p
                        (< --alt-min): the alt-carrying reads never arrived (clipped,
                        mismapped, or unmapped), so no caller could have called it.

Corroborating quantities per bin: mean f_alt = (ex_alt+sh_alt)/n_true_alt (alt-support
survival -- should be low where mapping-limited), and the split of the mapping failure
into "clipped over p" vs "mismapped/absent".

Output: one tidy CSV binned by local divergence (read scale, +/-150 bp), with a data
dictionary. This is the section's mechanism figure: it shows AT WHAT local divergence
each caller flips from calling-limited to mapping-limited.
"""

import os
import sys
import csv
import argparse
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "lib"))
import provenance                       # shared provenance block for the data dictionaries

# same local-divergence bin edges as local_recall.py (keep the axis consistent)
BIN_EDGES = [0.0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.13, 0.16, 0.20, 0.25, 0.30,
             0.40, 0.60, 1.01]


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--vrf-dir", required=True,
                   help="dir with vrf-pileup-<caller>.csv (from run_vrf_pileup.sh)")
    p.add_argument("--callers", nargs="+", default=["gatk", "bcftools", "freebayes"])
    p.add_argument("--localdiv", required=True,
                   help="mm39-truth-locald.tsv for this config (POS -> local_d_150)")
    p.add_argument("--window", default="150", help="local-D window column to bin on")
    p.add_argument("--alt-min", type=int, default=2,
                   help="min alt reads in the pileup to count the evidence as PRESENT "
                        "(>= this and uncalled = calling-limited). Default 2.")
    p.add_argument("--out", required=True)
    return p.parse_args()


def load_localdiv(path, window):
    """POS(1-based, int) -> local divergence at the requested window."""
    d = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ipos = header.index("POS")
        icol = header.index("local_d_%s" % window)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            d[int(f[ipos])] = float(f[icol])
    return d


def main():
    a = parse_args()
    ld = load_localdiv(a.localdiv, a.window)
    edges = np.array(BIN_EDGES)
    nb = len(edges) - 1
    rows_out = []

    for caller in a.callers:
        path = os.path.join(a.vrf_dir, "vrf-pileup-%s.csv" % caller)
        if not os.path.isfile(path):
            sys.stderr.write("WARN: missing %s (skipped)\n" % path)
            continue

        # accumulators keyed by (iteration, bin): counts + sums for means
        acc = {}
        with open(path) as fh:
            r = csv.DictReader(fh)
            for row in r:
                pos = int(row["snp_pos"])
                if pos not in ld:
                    continue                       # not a scored div site (shouldn't happen)
                b = int(np.clip(np.digitize(ld[pos], edges) - 1, 0, nb - 1))
                it = row["iteration"]
                rec = int(row["recovered"])
                alt_p = int(row["alt_p"])
                n_true = int(row["n_true_alt"]) or 1
                f_alt = (int(row["ex_alt"]) + int(row["sh_alt"])) / n_true
                clipped = (int(row["ex_clip"]) + int(row["sh_clip"])) / n_true
                mismapped = (int(row["mm_clip"]) + int(row["mm_noclip"])
                             + int(row["dc"]) + int(row["un"])) / n_true

                if rec:
                    cls = "recovered"
                elif alt_p >= a.alt_min:
                    cls = "calling"
                else:
                    cls = "mapping"

                key = (it, b)
                s = acc.setdefault(key, dict(n=0, recovered=0, calling=0, mapping=0,
                                             f_alt=0.0, clipped=0.0, mismapped=0.0,
                                             alt_p=0.0, depth_p=0.0))
                s["n"] += 1
                s[cls] += 1
                s["f_alt"] += f_alt
                s["clipped"] += clipped
                s["mismapped"] += mismapped
                s["alt_p"] += alt_p
                s["depth_p"] += int(row["depth_p"])

        for (it, b), s in sorted(acc.items()):
            n = s["n"]
            rows_out.append((
                caller, it, a.window, "%.4g" % edges[b], "%.4g" % edges[b + 1],
                "%.4g" % ((edges[b] + edges[b + 1]) / 2.0), n,
                s["recovered"], s["calling"], s["mapping"],
                s["recovered"] / n, s["calling"] / n, s["mapping"] / n,
                s["f_alt"] / n, s["clipped"] / n, s["mismapped"] / n,
                s["alt_p"] / n, s["depth_p"] / n))

    with open(a.out, "w") as out:
        out.write("caller,iteration,window,locald_lo,locald_hi,locald_mid,n_sites,"
                  "n_recovered,n_calling,n_mapping,"
                  "frac_recovered,frac_calling,frac_mapping,"
                  "mean_f_alt,mean_clipped_frac,mean_mismapped_frac,"
                  "mean_alt_p,mean_depth_p\n")
        for r in rows_out:
            out.write(("%s,%s,%s,%s,%s,%s,%d,%d,%d,%d,"
                       "%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n") % r)

    with open(a.out + ".dict.txt", "w") as d:
        d.write(provenance.header(a.out))
        d.write("s05 mapping-vs-calling decomposition, binned by local divergence.\n")
        d.write("One row per (caller, iteration, local-D bin). Alt-evidence threshold\n")
        d.write("--alt-min = %d reads.\n\n" % a.alt_min)
        d.write("Fields:\n")
        d.write("  caller, iteration : condition\n")
        d.write("  window            : local-D window half-size (bp) used for binning\n")
        d.write("  locald_lo/hi/mid  : local-divergence bin edges + midpoint (sites/bp)\n")
        d.write("  n_sites           : divergent truth SNPs sampled in this bin\n")
        d.write("  n_recovered/calling/mapping : site counts by class\n")
        d.write("  frac_recovered    : recovered / n_sites (== recall)\n")
        d.write("  frac_calling      : CALLING-limited (alt present >= alt-min but uncalled) / n\n")
        d.write("  frac_mapping      : MAPPING-limited (alt absent at p) / n\n")
        d.write("  mean_f_alt        : mean (ex_alt+sh_alt)/n_true_alt = alt-support survival\n")
        d.write("                      (low where mapping-limited: alt reads did not arrive)\n")
        d.write("  mean_clipped_frac : mean fraction of alt-carrying reads that arrived but\n")
        d.write("                      were soft-clipped over p\n")
        d.write("  mean_mismapped_frac : mean fraction mismapped/diff-chr/unmapped\n")
        d.write("  mean_alt_p/depth_p: mean alt-supporting reads / total depth in the pileup at p\n")
        d.write("\nfrac_recovered + frac_calling + frac_mapping = 1 within each bin.\n")

    sys.stderr.write("# wrote %d rows to %s\n" % (len(rows_out), a.out))


if __name__ == "__main__":
    main()
