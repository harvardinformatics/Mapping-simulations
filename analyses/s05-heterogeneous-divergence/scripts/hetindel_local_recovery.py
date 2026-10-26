#!/usr/bin/env python3
"""
hetindel_local_recovery.py -- the s05 het-indels separability test: SNP AND indel
recovery binned by local divergence, for every caller x iteration.

The separability question: heterogeneity relocates sites along the local-divergence
axis, and recovery-vs-local-divergence is (for SNPs) a property of the mapper/caller.
So the null is that heterogeneous recovery at local divergence x equals HOMOGENEOUS
(s02) recovery at MEAN divergence x. If the heterogeneous curve sits BELOW the
homogeneous one at matched local divergence, the effect is super-additive -- expected
for indels, since co-clustering makes hot regions locally indel-dense in a way s02's
uniform indels never are.

Everything is computed in truth.norm coordinates (the same normalized space the
accumulated recovered set `acc` lives in), so indel matching is exact:
  * div SNP/indel truth sites = normalized truth records with GT 1|1;
  * local divergence = divergent-SNP density (the field) at each site's position
    (SNP: self-inclusive window; indel: div-SNP density around the indel);
  * recovered = the site's (POS, REF, ALT) appears in that iteration's acc VCF
    (both truth and acc are bcftools norm -m- --atomize, so exact key match is valid).

Output: one tidy CSV (caller, iteration, vtype, local-D bin, n_sites, n_recovered,
recall) + a dictionary. Overlay against s02's homogeneous recall-vs-divergence.
"""

import os
import sys
import gzip
import argparse
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "lib"))
import divergence_field as dfield
import provenance                       # shared provenance block for the data dictionaries

BIN_EDGES = [0.0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.13, 0.16, 0.20, 0.25, 0.30,
             0.40, 0.60, 1.01]


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data-dir", required=True, help="hetindel data root (has sim/ and iter/)")
    p.add_argument("--ref-str", default="mm39")
    p.add_argument("--region", default="19")
    p.add_argument("--cov", default="30")
    p.add_argument("--div", default="0.06")
    p.add_argument("--het", default="0.005")
    p.add_argument("--seed", default="42")
    p.add_argument("--callers", nargs="+", default=["gatk", "bcftools", "freebayes"])
    p.add_argument("--iters", nargs="+", default=["1", "2", "3"])
    p.add_argument("--window", type=int, default=150, help="local-D window half-size (bp)")
    p.add_argument("--out", required=True)
    return p.parse_args()


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def load_truth(vcf):
    """Return dict: {(pos,ref,alt): (vtype, 0-based pos)} for GT 1|1 truth sites, and the
    sorted 0-based array of divergent-SNP positions (the field density source)."""
    sites = {}
    snp_pos = []
    with opener(vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            gt = f[9].split(":")[0]
            if gt not in ("1|1", "1/1"):
                continue
            pos, ref, alt = int(f[1]), f[3].upper(), f[4].upper()
            vtype = "snp" if (len(ref) == 1 and len(alt) == 1) else "indel"
            sites[(pos, ref, alt)] = (vtype, pos - 1)
            if vtype == "snp":
                snp_pos.append(pos - 1)
    return sites, np.sort(np.array(snp_pos, dtype=np.int64))


def load_acc(vcf):
    """Set of (pos, ref, alt) in an accumulated recovered VCF."""
    acc = set()
    with opener(vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            acc.add((int(f[1]), f[3].upper(), f[4].upper()))
    return acc


def main():
    a = parse_args()
    edges = np.array(BIN_EDGES)
    nb = len(edges) - 1

    simdir = os.path.join(a.data_dir, "sim", a.cov + "X", a.region, a.div, a.het, a.seed)
    truth = os.path.join(simdir, a.ref_str + "-truth.norm.vcf.gz")
    sites, snp_pos = load_truth(truth)

    # local divergence (div-SNP field) at every truth site's position, in one shot
    keys = list(sites.keys())
    pos0 = np.array([sites[k][1] for k in keys], dtype=np.int64)
    vtypes = np.array([sites[k][0] for k in keys])
    localn = dfield.windowed_density_at(pos0, snp_pos, a.window)   # density at every site
    binidx = np.clip(np.digitize(localn, edges) - 1, 0, nb - 1)

    rows = []
    for caller in a.callers:
        for it in a.iters:
            accf = os.path.join(a.data_dir, "iter", caller, a.cov + "X", a.region, a.div,
                                a.het, a.seed, "iter" + it, a.ref_str + "_recovered.acc.vcf.gz")
            if not os.path.isfile(accf):
                sys.stderr.write("WARN: missing acc VCF (skipped): %s\n" % accf)
                continue
            acc = load_acc(accf)
            recovered = np.fromiter((1 if k in acc else 0 for k in keys),
                                    dtype=np.int8, count=len(keys))
            for vt in ("snp", "indel"):
                sel = vtypes == vt
                if not sel.any():
                    continue
                for b in range(nb):
                    m = sel & (binidx == b)
                    n = int(m.sum())
                    if n == 0:
                        continue
                    nrec = int(recovered[m].sum())
                    rows.append((caller, it, vt,
                                 "%.4g" % edges[b], "%.4g" % edges[b + 1],
                                 "%.4g" % ((edges[b] + edges[b + 1]) / 2.0),
                                 n, nrec, nrec / n))

    with open(a.out, "w") as out:
        out.write("caller,iteration,vtype,locald_lo,locald_hi,locald_mid,"
                  "n_sites,n_recovered,recall\n")
        for r in rows:
            out.write("%s,%s,%s,%s,%s,%s,%d,%d,%.6g\n" % r)

    with open(a.out + ".dict.txt", "w") as d:
        d.write(provenance.header(a.out))
        d.write("s05 het-indels: SNP and indel recovery binned by local divergence.\n")
        d.write("One row per (caller, iteration, vtype, local-D bin). Window = +/-%d bp.\n\n" % a.window)
        d.write("Fields:\n")
        d.write("  caller, iteration : condition (mean divergence fixed at %s)\n" % a.div)
        d.write("  vtype             : snp | indel (both GT 1|1 divergence truth sites)\n")
        d.write("  locald_lo/hi/mid  : local-divergence bin edges + midpoint (div-SNP density)\n")
        d.write("  n_sites           : truth sites of this vtype in the bin\n")
        d.write("  n_recovered       : matched (POS,REF,ALT) in the accumulated recovered VCF\n")
        d.write("  recall            : n_recovered / n_sites\n\n")
        d.write("Overlay against s02 homogeneous recall-vs-mean-divergence: het recall at\n")
        d.write("local divergence x below s02 recall at mean divergence x => super-additive.\n")

    sys.stderr.write("# %d truth sites (%d snp, %d indel); wrote %d rows\n"
                     % (len(keys), int((vtypes == "snp").sum()),
                        int((vtypes == "indel").sum()), len(rows)))


if __name__ == "__main__":
    main()
