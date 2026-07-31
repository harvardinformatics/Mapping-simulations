#!/usr/bin/env python3
"""
make_loss_decomposition.py -- decompose SNP reference-bias loss into its two
mechanisms (caller-limited vs mapping-limited) as a function of divergence, for the
loss-mechanism heatmap.

Both inputs come from the SAME (SNP-only) run so the decomposition is internally
consistent:
  * overall SNP loss rate per divergence   -> per-condition compare_vcfs summaries
    (div_ref_or_absent / n_truth_div)
  * split of lost sites into caller- vs mapping-limited -> the site-linkage table,
    using effective ALT fraction (MAPQ>=20): >=50% = evidence present = caller-limited
    (a pileup caller recovers these); <50% = evidence degraded = mapping-limited.

Output: a long-format CSV (one row per divergence x mechanism) with the count and the
fraction of ALL true divergence sites lost to that mechanism, plus a data dictionary.
"""

import csv
import glob
import os
from collections import defaultdict

BASE = "data/mm39-mason-chr19"
LINK = "summary-data/mm39-mason-chr19-site-linkage.csv"
OUT = "summary-data/mm39-mason-chr19-loss-decomposition.csv"
DIVS = ["0.00", "0.02", "0.04", "0.06", "0.08", "0.10"]
VAF_CUT = 0.5  # effective ALT fraction separating caller- vs mapping-limited loss


def overall_loss():
    """divergence -> (n_truth_div, n_lost) from the SNP-only per-condition summaries."""
    out = {}
    for d in DIVS:
        f = os.path.join(BASE, "summary-files", "30X", "19", d, "0.005", "mm39-vcf-summary.csv")
        if not os.path.isfile(f):
            out[d] = (0, 0)
            continue
        row = list(csv.DictReader(open(f)))[0]
        out[d] = (int(row["n_truth_div"]), int(row["div_ref_or_absent"]))
    return out


def caller_fraction():
    """divergence -> fraction of LOST sites that are caller-limited (ALT frac >= cut)."""
    present = defaultdict(int)
    total = defaultdict(int)
    for r in csv.DictReader(open(LINK)):
        if r["outcome"] != "lost":
            continue
        d = r["divergence"]
        total[d] += 1
        if float(r["pileup_alt_frac_hq"]) >= VAF_CUT:
            present[d] += 1
    return {d: (present[d] / total[d] if total[d] else 0.0) for d in total}


def main():
    loss = overall_loss()
    cfrac = caller_fraction()

    rows = []
    for d in DIVS:
        n_truth, n_lost = loss[d]
        cf = cfrac.get(d, 0.0)
        caller_lost = round(n_lost * cf)
        mapping_lost = n_lost - caller_lost
        for mech, n in [("caller-limited", caller_lost), ("mapping-limited", mapping_lost)]:
            frac = n / n_truth if n_truth else 0.0
            rows.append((d, mech, n, "%.6f" % frac))

    with open(OUT, "w") as out:
        out.write("divergence,mechanism,n_sites_lost,frac_of_truth_lost\n")
        for r in rows:
            out.write(",".join(str(x) for x in r) + "\n")

    with open(OUT + ".dict.txt", "w") as dct:
        dct.write("Data dictionary for %s\n" % OUT)
        dct.write("SNP reference-bias loss decomposed by mechanism vs divergence (SNP-only run).\n\n")
        dct.write("divergence        : simulated fixed SNP divergence D\n")
        dct.write("mechanism         : caller-limited (alt evidence present, ALT frac >= %.2f;\n" % VAF_CUT)
        dct.write("                    a pileup caller recovers these) vs mapping-limited (alt\n")
        dct.write("                    evidence degraded, ALT frac < %.2f; needs better mapping)\n" % VAF_CUT)
        dct.write("n_sites_lost      : number of true divergence sites lost to that mechanism\n")
        dct.write("                    (overall lost count x that mechanism's share of lost sites)\n")
        dct.write("frac_of_truth_lost: n_sites_lost / n_truth_div (fraction of ALL true div sites)\n")
        dct.write("\nProvenance: overall loss from per-condition compare_vcfs summaries; the\n")
        dct.write("caller/mapping split from the site-linkage table (both SNP-only run).\n")

    print("wrote %s" % OUT)
    for r in rows:
        print("  ", r)


if __name__ == "__main__":
    main()
