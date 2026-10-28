#!/usr/bin/env python3
"""
cluster_recall.py -- recall of divergent truth SNPs binned by how CLUSTERED they are: the number of
OTHER divergent truth SNPs within +/-W bp of each site.

WHY THIS EXISTS
---------------
GATK HaplotypeCaller re-assembles each active region from k-mers (default k = 10 and 25) anchored on
the reference. Close-packed SNPs leave few k-mers shared with the reference, so its assembly -- and
the call -- can fail even when the reads carrying the alleles are present. bcftools and freebayes
call from the aligned reads and should be insensitive to clustering until the reads themselves stop
aligning. Binning recall by the local SNP count at the k-mer scale tests this directly, at a FIXED
genome-wide divergence (uniform simulations), where the cluster count varies only by chance.

DEFINITIONS
-----------
  site        : a divergent truth SNP -- a record of the normalized truth VCF whose ID starts with
                --site-id (default "div_snp", the homozygous 1|1 divergence SNPs)
  n_near      : number of OTHER sites (same definition) whose position is within +/-W bp of the site
  recovered   : the site's (POS, ALT) is in that iteration's recovered VCF (the cumulative set of
                SNPs written into the consensus, the same set the section summary scores as div_tp)
  recall      : recovered sites / sites, per (caller, divergence, seed, iteration, window, n_near bin)

--check-summary verifies, for every condition, that the recovered count summed over the bins equals
div_tp in the section's summary table (and the site count equals n_truth_div), and stops otherwise.

Usage (from the project root; see analyses/s01-baseline-iterative/run_cluster_recall.sh):
  python analyses/lib/cluster_recall.py \
      --truth "data/.../sim/30X/19/{D}/0.005/{seed}/mm39-truth.norm.vcf.gz" \
      --recovered "data/.../iter/{caller}/30X/19/{D}/0.005/{seed}/iter{k}/mm39_recovered.vcf.gz" \
      --divs 0.02,0.04 --seeds 42,43,44 --iterations 1,2,3 --windows 10,25 \
      --check-summary analyses/.../s01-base-iterative-summary.csv --out OUT.csv
"""
import argparse
import gzip
import os
import sys
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import provenance                       # shared provenance block for the data dictionaries

# --- n_near bins (upper edges inclusive; the last bin is open-ended) -------------------------------
BIN_EDGES = (0, 1, 2, 3, 5)             # bins: 0 | 1 | 2 | 3 | 4-5 | 6+
BIN_LABELS = ("0", "1", "2", "3", "4-5", "6+")
BASES = {"A": 0, "C": 1, "G": 2, "T": 3}

COLUMNS = ("caller", "divergence", "seed", "iteration", "window_bp", "n_near_bin", "n_sites",
           "n_recovered", "recall")


def snp_keys(path, site_id=None):
    """Sorted positions and (position, alt) keys of the single-base biallelic records of a VCF.
    With site_id, only records whose ID starts with it. Fails on an unknown ALT base."""
    pos, alt = [], []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            c = line.split("\t", 6)
            if site_id and not c[2].startswith(site_id):
                continue
            if len(c[3]) != 1 or len(c[4]) != 1:
                continue                # indels / MNPs: not SNP sites
            if c[4] not in BASES:
                raise ValueError("%s: ALT %r at %s is not a single base" % (path, c[4], c[1]))
            pos.append(int(c[1])); alt.append(BASES[c[4]])
    pos, alt = np.array(pos, dtype=np.int64), np.array(alt, dtype=np.int64)
    if np.any(np.diff(pos) < 0):
        raise ValueError("%s: records not sorted by position" % path)
    return pos, pos * 4 + alt


def n_near(pos, w):
    """For each site, the number of OTHER sites within +/-w bp."""
    return np.searchsorted(pos, pos + w, "right") - np.searchsorted(pos, pos - w, "left") - 1


def bin_index(n):
    """Bin index per site (0 .. len(BIN_LABELS) - 1)."""
    return np.searchsorted(np.array(BIN_EDGES), n, "left")


def one_condition(job):
    """Rows for one (divergence, seed): every caller x iteration x window."""
    a, d, seed = job
    pos, keys = snp_keys(a.truth.format(D=d, seed=seed), a.site_id)
    if len(pos) == 0:
        raise ValueError("no sites in the truth VCF for D=%s seed=%s" % (d, seed))
    bins = {w: bin_index(n_near(pos, w)) for w in a.windows}
    rows = []
    for caller in a.callers:
        for k in a.iterations:
            _, rec = snp_keys(a.recovered.format(caller=caller, D=d, seed=seed, k=k))
            hit = np.isin(keys, rec)
            for w in a.windows:
                n = np.bincount(bins[w], minlength=len(BIN_LABELS))
                r = np.bincount(bins[w], weights=hit, minlength=len(BIN_LABELS)).astype(np.int64)
                for i, lab in enumerate(BIN_LABELS):
                    rows.append((caller, d, seed, k, w, lab, int(n[i]), int(r[i]),
                                 "%.6f" % (r[i] / n[i]) if n[i] else "NA"))
    sys.stderr.write("# D=%s seed=%s: %d sites\n" % (d, seed, len(pos)))
    return rows


def check_summary(rows, path):
    """Total sites / recovered per condition must equal n_truth_div / div_tp in the summary table."""
    import csv
    want = {}
    with open(path) as fh:
        for r in csv.DictReader(fh):
            want[(r["caller"], r["divergence"], r["seed"], r["iteration"])] = (int(r["n_truth_div"]),
                                                                                int(r["div_tp"]))
    got = {}
    for caller, d, seed, k, w, _, n, rec, _ in rows:
        key = (caller, d, seed, str(k), w)
        s = got.setdefault(key, [0, 0]); s[0] += n; s[1] += rec
    bad = []
    for (caller, d, seed, k, w), (n, rec) in sorted(got.items()):
        exp = want.get((caller, d, seed, k))
        if exp is None or exp != (n, rec):
            bad.append("%s D=%s seed=%s iter=%s W=%d: sites/recovered %d/%d, summary %s"
                       % (caller, d, seed, k, w, n, rec, exp))
    if bad:
        raise SystemExit("--check-summary FAILED:\n  " + "\n  ".join(bad))
    sys.stderr.write("# --check-summary: all %d conditions match %s\n" % (len(got), path))


def write_dictionary(a):
    body = ("Recall of divergent truth SNPs binned by the number of OTHER divergent truth SNPs within\n"
            "+/-W bp (see the cluster_recall.py docstring for the rationale).\n\n"
            "Sites      : records of the normalized truth VCF (template %s) whose ID starts\n"
            "             with %r and REF/ALT are single bases.\n"
            "Recovered  : the site's (POS, ALT) is in that iteration's recovered VCF (template %s),\n"
            "             the cumulative set of SNPs written into the consensus.%s\n"
            "Windows    : W = %s bp (half-window; n_near counts other sites in [pos-W, pos+W]).\n"
            "Bins       : n_near = %s (the last bin is open-ended).\n"
            "Conditions : callers %s; divergences %s; seeds %s; iterations %s.\n\n"
            "Fields:\n"
            "  caller      : variant caller (gatk | bcftools | freebayes)\n"
            "  divergence  : simulated divergence D (label as in the data paths)\n"
            "  seed        : simulation seed\n"
            "  iteration   : mapping iteration (1 = single-shot)\n"
            "  window_bp   : half-window W (bp)\n"
            "  n_near_bin  : bin of n_near, the number of OTHER sites within +/-W bp\n"
            "  n_sites     : sites in the bin\n"
            "  n_recovered : of those, recovered at this iteration\n"
            "  recall      : n_recovered / n_sites (NA when the bin is empty)\n"
            % (a.truth, a.site_id, a.recovered,
               "\n             Checked against %s: per condition, the summed sites and\n"
               "             recovered equal its n_truth_div and div_tp." % a.check_summary
               if a.check_summary else "",
               ", ".join(str(w) for w in a.windows), " | ".join(BIN_LABELS), ", ".join(a.callers),
               ", ".join(a.divs), ", ".join(a.seeds), ", ".join(str(k) for k in a.iterations)))
    provenance.write_dictionary(a.out, body)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", required=True, help="normalized truth VCF template; {D}, {seed}")
    ap.add_argument("--recovered", required=True, help="recovered VCF template; {caller}, {D}, {seed}, {k}")
    ap.add_argument("--divs", required=True, help="comma-separated divergence labels")
    ap.add_argument("--seeds", required=True, help="comma-separated seeds")
    ap.add_argument("--callers", default="gatk,bcftools,freebayes")
    ap.add_argument("--iterations", default="1,2,3")
    ap.add_argument("--windows", default="10,25", help="comma-separated half-windows W (bp)")
    ap.add_argument("--site-id", default="div_snp", help="ID prefix of the truth records used as sites")
    ap.add_argument("--check-summary", help="section summary table with n_truth_div / div_tp to verify against")
    ap.add_argument("--threads", type=int, default=8, help="(divergence, seed) conditions run in parallel")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    a.divs, a.seeds, a.callers = a.divs.split(","), a.seeds.split(","), a.callers.split(",")
    a.iterations = [int(k) for k in a.iterations.split(",")]
    a.windows = [int(w) for w in a.windows.split(",")]

    jobs = [(a, d, s) for d in a.divs for s in a.seeds]
    with Pool(min(a.threads, len(jobs))) as pool:
        rows = [r for part in pool.map(one_condition, jobs) for r in part]
    if a.check_summary:
        check_summary(rows, a.check_summary)

    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w") as out:
        out.write(",".join(COLUMNS) + "\n")
        for r in rows:
            out.write(",".join(str(x) for x in r) + "\n")
    write_dictionary(a)
    sys.stderr.write("# wrote %d rows to %s\n" % (len(rows), a.out))


if __name__ == "__main__":
    main()
