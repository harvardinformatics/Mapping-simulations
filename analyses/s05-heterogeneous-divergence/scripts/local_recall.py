#!/usr/bin/env python3
"""
local_recall.py -- the s05 headline table: SNP recovery binned by MEASURED LOCAL
divergence, for every field configuration x caller x iteration.

For each divergent (1|1) truth SNP we know (a) its local divergence at the read
(150 bp) and fragment (500 bp) window scales (from annotate_local_divergence.py's
mm39-truth-locald.tsv) and (b) whether the pipeline recovered it at a given
iteration (its position appears in that iteration's mm39_recovered.vcf.gz, the
cumulative consensus-vs-reference SNP diff). Recovery is a pure position match:
s05 is SNP-only so the consensus shares reference coordinates (no liftover), and
the recovered set is exactly the truth positions the consensus fixed.

This script scans the s05 data tree and emits ONE tidy CSV:
  field, corr_len, sigma, caller, iteration, window, locald_lo, locald_hi,
  locald_mid, n_sites, n_recovered, recall
plus a whole-config (unbinned) recall row (window="all", locald bins = NA) as a
cross-check against the aggregate summary CSVs.

Binning is EXPLICIT (edges below), finer where recovery changes fastest. Every
site falls in exactly one bin; sites at local D == the top edge go in the last bin.
"""

import os
import sys
import gzip
import argparse
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "lib"))
import provenance                       # shared provenance block for the data dictionaries

# Local-divergence bin edges (sites/bp). Finer through the 0.05-0.30 transition
# where recovery falls off; the last bin catches the diverged cores.
BIN_EDGES = [0.0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.13, 0.16, 0.20, 0.25, 0.30,
             0.40, 0.60, 1.01]


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data-dir", required=True, help="s05 data root (contains sim/ and iter/)")
    p.add_argument("--ref-str", default="mm39")
    p.add_argument("--region", default="19")
    p.add_argument("--cov", default="30")
    p.add_argument("--div", default="0.06")
    p.add_argument("--het", default="0.005")
    p.add_argument("--seed", default="42")
    p.add_argument("--configs", nargs="+", required=True, help="field tags, e.g. L20000-S1.0")
    p.add_argument("--callers", nargs="+", default=["gatk", "bcftools", "freebayes"])
    p.add_argument("--iters", nargs="+", default=["1", "2", "3"])
    p.add_argument("--windows", nargs="+", default=["150", "500"],
                   help="local-D window columns to bin on (must exist in the locald TSV)")
    p.add_argument("--out", required=True, help="output tidy CSV")
    return p.parse_args()


def parse_field(tag):
    # 'L20000-S1.0' or 'L20000-S1.0-C0.20' -> (corr_len, sigma); cap (if any) ignored here
    parts = tag.split("-")
    return parts[0][1:], parts[1][1:]


def load_locald(path, windows):
    """Return (pos int64 array, {window: local_d float array})."""
    pos = []
    cols = {w: [] for w in windows}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {w: header.index("local_d_%s" % w) for w in windows}
        ipos = header.index("POS")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            pos.append(int(f[ipos]))
            for w in windows:
                cols[w].append(float(f[idx[w]]))
    return (np.array(pos, dtype=np.int64),
            {w: np.array(cols[w], dtype=np.float64) for w in windows})


def load_recovered_positions(path):
    """Sorted int64 array of 1-based POS in a bgzipped recovered VCF."""
    rec = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            rec.append(int(line.split("\t", 2)[1]))
    return np.sort(np.array(rec, dtype=np.int64))


def main():
    args = parse_args()
    edges = np.array(BIN_EDGES)
    rows = []

    for cfg in args.configs:
        corr, sigma = parse_field(cfg)
        sim = os.path.join(args.data_dir, "sim", args.cov + "X", args.region,
                           args.div, args.het, cfg, args.seed,
                           args.ref_str + "-truth-locald.tsv")
        if not os.path.isfile(sim):
            sys.exit("ERROR: missing local-D table: %s" % sim)
        pos, local = load_locald(sim, args.windows)
        # precompute bin index per window (independent of caller/iter)
        binidx = {w: np.clip(np.digitize(local[w], edges) - 1, 0, len(edges) - 2)
                  for w in args.windows}

        for caller in args.callers:
            for it in args.iters:
                rvcf = os.path.join(args.data_dir, "iter", caller, args.cov + "X",
                                    args.region, args.div, args.het, cfg, args.seed,
                                    "iter" + it, args.ref_str + "_recovered.vcf.gz")
                if not os.path.isfile(rvcf):
                    sys.stderr.write("WARN: missing recovered VCF (skipped): %s\n" % rvcf)
                    continue
                rec = load_recovered_positions(rvcf)
                recovered = np.isin(pos, rec, assume_unique=True)  # vectorized membership

                # whole-config recall (cross-check row)
                rows.append((cfg, corr, sigma, caller, it, "all", "", "", "",
                             pos.size, int(recovered.sum()),
                             recovered.mean()))

                # binned recall per window
                for w in args.windows:
                    bi = binidx[w]
                    for b in range(len(edges) - 1):
                        m = bi == b
                        n = int(m.sum())
                        if n == 0:
                            continue
                        nrec = int(recovered[m].sum())
                        rows.append((cfg, corr, sigma, caller, it, w,
                                     "%.4g" % edges[b], "%.4g" % edges[b + 1],
                                     "%.4g" % ((edges[b] + edges[b + 1]) / 2.0),
                                     n, nrec, nrec / n))

    with open(args.out, "w") as out:
        out.write("field,corr_len,sigma,caller,iteration,window,"
                  "locald_lo,locald_hi,locald_mid,n_sites,n_recovered,recall\n")
        for r in rows:
            out.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%d,%d,%.6g\n" % r)

    with open(args.out + ".dict.txt", "w") as d:
        d.write(provenance.header(args.out))
        d.write("s05 SNP recovery binned by measured local divergence.\n")
        d.write("One row per (field config, caller, iteration, window, local-D bin), plus\n")
        d.write("a window='all' whole-config cross-check row per (config, caller, iteration).\n\n")
        d.write("Fields:\n")
        d.write("  field       : divergence-field tag L{corr_len}-S{sigma}\n")
        d.write("  corr_len    : field correlation length (bp); scale of heterogeneity\n")
        d.write("  sigma       : field log-amplitude; 0 = homogeneous control\n")
        d.write("  caller      : gatk | bcftools | freebayes\n")
        d.write("  iteration   : mapping iteration (1 = single-shot; 2..N iterative)\n")
        d.write("  window      : local-D window half-size in bp (150 read, 500 fragment); 'all' = whole config\n")
        d.write("  locald_lo/hi: local-divergence bin edges (sites/bp); blank for window='all'\n")
        d.write("  locald_mid  : bin midpoint (plotting x); blank for window='all'\n")
        d.write("  n_sites     : truth divergent (1|1) SNP sites in the bin\n")
        d.write("  n_recovered : of those, recovered at this iteration (in the consensus diff)\n")
        d.write("  recall      : n_recovered / n_sites\n\n")
        d.write("Local divergence per site = divergent-site density in +/- window bp (from\n")
        d.write("annotate_local_divergence.py). Mean divergence is FIXED at div=%s across all\n" % args.div)
        d.write("configs, so differences between configs are due to spatial STRUCTURE only.\n")
        d.write("Bin edges: %s\n" % ", ".join("%g" % e for e in BIN_EDGES))

    sys.stderr.write("# wrote %d rows to %s\n" % (len(rows), args.out))


if __name__ == "__main__":
    main()
