#!/usr/bin/env python3
"""
mappability_recovery.py -- the s08 multi-mapping confounder: small-variant recovery
binned by local REFERENCE mappability.

Unlike an SV (a sample-side change), the "can't uniquely place a read" failure is a
property of the REFERENCE: where chr19 contains repeated/near-identical sequence, a
150 bp read fits several places, bwa assigns low MAPQ, the MQ>=20 filter drops it, and
the variant there cannot be called -- no matter how similar the sample is. This is the
mechanism behind segmental-duplication / difficult-region variant-calling failures
(the sim analog of GIAB's difficulty stratification).

We measure it with the real repeat structure of chr19: GenMap mappability (K=100, E=2 =
1 / number of near-identical 100-mer copies within chr19; 1.0 = unique, 0.5 = 2 copies,
0.33 = 3, ...). Each divergent SNP truth site (from s01, SNP-only, so the divergence
axis is clean) is labeled with the mappability at its position and its recovery is
scored against the recovered set, per caller x divergence x iteration.

Prediction: recovery is at baseline in unique sequence and collapses as mappability
falls -- a divergence-INDEPENDENT loss (it is about read placement ambiguity, not
sample similarity), stacking on top of the divergence-driven reference bias.
"""

import os
import sys
import gzip
import argparse
import numpy as np

# 3 monotone bins: unique, moderately repetitive, highly repetitive. (Finer splits of the
# sparse middle are non-monotone -- the exactly-2-copy class behaves oddly -- so we keep
# the robust three-way contrast.)
MAP_EDGES = [0.0, 0.2, 0.999, 1.001]
MAP_LABELS = ["highly repetitive (<0.2)", "repetitive (0.2-1.0)", "unique (1.0)"]


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bedgraph", required=True, help="GenMap mappability bedgraph (chr19)")
    p.add_argument("--data-dir", required=True, help="s01 data root")
    p.add_argument("--ref-str", default="mm39")
    p.add_argument("--region", default="19")
    p.add_argument("--cov", default="30")
    p.add_argument("--het", default="0.005")
    p.add_argument("--seed", default="42")
    p.add_argument("--divs", nargs="+", required=True)
    p.add_argument("--callers", nargs="+", default=["gatk", "bcftools", "freebayes"])
    p.add_argument("--iters", nargs="+", default=["1", "3"])
    p.add_argument("--out", required=True)
    return p.parse_args()


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def load_bedgraph(path, region):
    """Return sorted (starts, ends, vals) int64/float arrays for the region (0-based)."""
    starts, ends, vals = [], [], []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] != region:
                continue
            starts.append(int(f[1])); ends.append(int(f[2])); vals.append(float(f[3]))
    order = np.argsort(np.array(starts))
    return (np.array(starts)[order], np.array(ends)[order], np.array(vals)[order])


def mappability_at(pos0, starts, ends, vals):
    """Mappability value at each 0-based position (np.nan if uncovered, e.g. N runs)."""
    idx = np.clip(np.searchsorted(starts, pos0, side="right") - 1, 0, len(starts) - 1)
    inside = (starts[idx] <= pos0) & (pos0 < ends[idx])
    out = np.where(inside, vals[idx], np.nan)
    return out


def load_div_snp_pos(truth_norm):
    """1-based POS of GT 1|1 SNP truth sites (s01 is SNP-only)."""
    pos = []
    with opener(truth_norm) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f[3]) == 1 and len(f[4]) == 1 and f[9].split(":")[0] in ("1|1", "1/1"):
                pos.append(int(f[1]))
    return np.array(pos, dtype=np.int64)


def load_recovered_pos(path):
    rec = []
    with opener(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            rec.append(int(line.split("\t", 2)[1]))
    return np.sort(np.array(rec, dtype=np.int64))


def main():
    a = parse_args()
    starts, ends, vals = load_bedgraph(a.bedgraph, a.region)
    rows = []

    for div in a.divs:
        sim = os.path.join(a.data_dir, "sim", a.cov + "X", a.region, div, a.het, a.seed)
        truth = os.path.join(sim, a.ref_str + "-truth.norm.vcf.gz")
        if not os.path.isfile(truth):
            sys.stderr.write("WARN: missing truth %s\n" % truth); continue
        pos = load_div_snp_pos(truth)
        mp = mappability_at(pos - 1, starts, ends, vals)      # 0-based lookup
        valid = ~np.isnan(mp)
        binidx = np.clip(np.digitize(mp, np.array(MAP_EDGES)) - 1, 0, len(MAP_EDGES) - 2)

        for caller in a.callers:
            for it in a.iters:
                rvcf = os.path.join(a.data_dir, "iter", caller, a.cov + "X", a.region, div,
                                    a.het, a.seed, "iter" + it, a.ref_str + "_recovered.vcf.gz")
                if not os.path.isfile(rvcf):
                    sys.stderr.write("WARN: missing recovered %s\n" % rvcf); continue
                rec = load_recovered_pos(rvcf)
                recovered = np.isin(pos, rec, assume_unique=True)
                for b in range(len(MAP_EDGES) - 1):
                    m = valid & (binidx == b)
                    n = int(m.sum())
                    if n == 0:
                        continue
                    nrec = int(recovered[m].sum())
                    rows.append((caller, div, it, MAP_LABELS[b], b, n, nrec, nrec / n))

    with open(a.out, "w") as out:
        out.write("caller,divergence,iteration,map_bin,map_rank,n_sites,n_recovered,recall\n")
        for r in rows:
            out.write("%s,%s,%s,%s,%d,%d,%d,%.6g\n" % r)

    with open(a.out + ".dict.txt", "w") as d:
        d.write("Data dictionary for %s\n\n" % a.out)
        d.write("s08 multi-mapping confounder: SNP recovery binned by local reference\n")
        d.write("mappability (GenMap K=100 E=2 on chr19). One row per (caller, divergence,\n")
        d.write("iteration, mappability bin). s01 SNP-only data.\n\n")
        d.write("Fields:\n")
        d.write("  caller, divergence, iteration : condition\n")
        d.write("  map_bin   : mappability bin (1.0 = unique; 0.5 ~ 2 near-identical copies\n")
        d.write("              within chr19; <0.2 ~ >=5 copies = highly repetitive)\n")
        d.write("  map_rank  : 0..4 ordinal of the bin (0 = most repetitive, 4 = unique)\n")
        d.write("  n_sites, n_recovered, recall : truth SNPs / recovered / recall in the bin\n\n")
        d.write("Mappability is the value at the SNP position (the 100-mer starting there).\n")
        d.write("Low mappability = reads there multi-map -> low MAPQ -> dropped -> uncallable,\n")
        d.write("independent of divergence.\n")

    sys.stderr.write("# %d rows; bedgraph intervals=%d\n" % (len(rows), starts.size))


if __name__ == "__main__":
    main()
