#!/usr/bin/env python3
"""
sv_confounder.py -- the s07 SV-confounder analysis: small-variant recovery as a
function of DISTANCE from the nearest SV breakpoint.

SVs are never recovered (our callers can't emit them), so they matter only by
disrupting read mapping in their neighborhood. This bins every SNP / small-indel
truth site by its distance to the nearest SV breakpoint and asks what fraction was
recovered, per caller x iteration. The prediction: recovery craters within a read/
fragment length of a breakpoint (reads clip / don't map) and returns to the genome
baseline far away -- and iteration lifts the neighborhood but not the breakpoint
itself (the SV stays unrepresented in the consensus).

COORDINATES (all reference):
  * SV breakpoints from the mason truth VCF (div_sv_*): a deletion contributes TWO
    junctions (its start POS and end POS+len(REF)-1); an insertion contributes one
    (its anchor POS). Reads clip at these junctions.
  * Divergent SNP/indel truth sites from truth.norm (GT 1|1), matched to the
    accumulated recovered VCF (acc) by (POS,REF,ALT) -- both bcftools-normalized.
  * Distance = min |site POS - any breakpoint|. Sites are disjoint from SV spans by
    construction (make_truth_vcf places SVs first), so distances are >= 1.

OUTPUT: tidy CSV (caller, iteration, vtype, dist_lo, dist_hi, dist_mid, n_sites,
n_recovered, recall) + dictionary. The far bin is the SV-free baseline.
"""

import os
import sys
import gzip
import argparse
import numpy as np

DIST_EDGES = [0, 50, 100, 150, 250, 400, 700, 1200, 2500, 5000, 10**9]


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data-dir", required=True, help="s07 data root (sim/ and iter/)")
    p.add_argument("--ref-str", default="mm39")
    p.add_argument("--region", default="19")
    p.add_argument("--cov", default="30")
    p.add_argument("--div", default="0.02")
    p.add_argument("--het", default="0.005")
    p.add_argument("--seed", default="42")
    p.add_argument("--callers", nargs="+", default=["gatk", "bcftools", "freebayes"])
    p.add_argument("--iters", nargs="+", default=["1", "2", "3"])
    p.add_argument("--out", required=True)
    return p.parse_args()


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def load_sv_breakpoints(mason_vcf):
    """Return (positions int64, types) of reference breakpoints from div_sv_* records,
    sorted by position. A deletion contributes two junctions (both 'del'); an insertion
    one ('ins')."""
    pos, typ = [], []
    with opener(mason_vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            vid = f[2]
            if not vid.startswith("div_sv"):
                continue
            p, ref = int(f[1]), f[3]
            if "_del_" in vid:
                pos.append(p); typ.append("del")
                pos.append(p + len(ref) - 1); typ.append("del")
            else:
                pos.append(p); typ.append("ins")
    pos = np.array(pos, dtype=np.int64)
    typ = np.array(typ)
    order = np.argsort(pos)
    return pos[order], typ[order]


def load_truth_div(truth_norm):
    """Return dict {(pos,ref,alt): (vtype, pos)} for GT 1|1 SNP/indel truth sites."""
    sites = {}
    with opener(truth_norm) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2].startswith("div_sv"):
                continue                           # SVs are the confounder, not a small-var site
            if f[9].split(":")[0] not in ("1|1", "1/1"):
                continue
            pos, ref, alt = int(f[1]), f[3].upper(), f[4].upper()
            vt = "snp" if (len(ref) == 1 and len(alt) == 1) else "indel"
            sites[(pos, ref, alt)] = (vt, pos)
    return sites


def load_acc(path):
    acc = set()
    with opener(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            acc.add((int(f[1]), f[3].upper(), f[4].upper()))
    return acc


def load_keyset(path, exclude_sv=False):
    """Set of 'POS:REF:ALT' for all variant records (SNP/indel/large), optionally
    excluding div_sv records. Used for the FP-near-breakpoint enrichment."""
    s = set()
    with opener(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if exclude_sv and f[2].startswith("div_sv"):
                continue
            s.add("%s:%s:%s" % (f[1], f[3].upper(), f[4].upper()))
    return s


def called_positions(path):
    """(sorted int64 positions, 'POS:REF:ALT' keys) for all called records."""
    pos, keys = [], []
    with opener(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            pos.append(int(f[1]))
            keys.append("%s:%s:%s" % (f[1], f[3].upper(), f[4].upper()))
    return np.array(pos, dtype=np.int64), np.array(keys)


def contig_length(data_dir, ref_str, region):
    fai = os.path.join(data_dir, "ref-regions", "%s-%s.fa.fai" % (ref_str, region))
    with open(fai) as fh:
        return int(fh.readline().split("\t")[1])


def nearest_breakpoint(pos, bp_pos, bp_typ):
    """Return (distance, nearest-breakpoint-type) from each pos to the closest
    breakpoint in sorted bp_pos."""
    n = len(bp_pos)
    idx = np.searchsorted(bp_pos, pos)
    ld = np.where(idx > 0, pos - bp_pos[np.clip(idx - 1, 0, n - 1)], np.inf)
    rd = np.where(idx < n, bp_pos[np.clip(idx, 0, n - 1)] - pos, np.inf)
    dist = np.minimum(ld, rd)
    near_idx = np.where(ld <= rd, np.clip(idx - 1, 0, n - 1), np.clip(idx, 0, n - 1))
    return dist, bp_typ[near_idx]


def main():
    a = parse_args()
    sim = os.path.join(a.data_dir, "sim", a.cov + "X", a.region, a.div, a.het, a.seed)
    bp_pos, bp_typ = load_sv_breakpoints(os.path.join(sim, a.ref_str + "-truth.vcf.gz"))
    if bp_pos.size == 0:
        sys.exit("ERROR: no div_sv breakpoints found -- is this an SV run?")
    sites = load_truth_div(os.path.join(sim, a.ref_str + "-truth.norm.vcf.gz"))

    keys = list(sites.keys())
    pos = np.array([sites[k][1] for k in keys], dtype=np.int64)
    vtypes = np.array([sites[k][0] for k in keys])
    dist, near_typ = nearest_breakpoint(pos, bp_pos, bp_typ)
    edges = np.array(DIST_EDGES, dtype=float)
    binidx = np.clip(np.digitize(dist, edges) - 1, 0, len(edges) - 2)

    def emit(rows, caller, it, recovered, sel, bptype):
        for b in range(len(edges) - 1):
            m = sel & (binidx == b)
            n = int(m.sum())
            if n == 0:
                continue
            nrec = int(recovered[m].sum())
            hi = edges[b + 1] if edges[b + 1] < 10**8 else -1     # -1 = far/baseline
            mid = edges[b] + (edges[b + 1] - edges[b]) / 2.0 if hi > 0 else edges[b]
            rows.append((caller, it, vt, bptype, "%g" % edges[b], "%g" % hi,
                         "%g" % mid, n, nrec, nrec / n))

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
                if not (vtypes == vt).any():
                    continue
                emit(rows, caller, it, recovered, vtypes == vt, "any")             # all breakpoints
                emit(rows, caller, it, recovered, (vtypes == vt) & (near_typ == "del"), "del")
                emit(rows, caller, it, recovered, (vtypes == vt) & (near_typ == "ins"), "ins")

    with open(a.out, "w") as out:
        out.write("caller,iteration,vtype,bptype,dist_lo,dist_hi,dist_mid,"
                  "n_sites,n_recovered,recall\n")
        for r in rows:
            out.write("%s,%s,%s,%s,%s,%s,%s,%d,%d,%.6g\n" % r)

    # --- FP enrichment near breakpoints (single-shot iter1 raw calls) -----------
    # Do SVs also CREATE false calls near breakpoints? Compare the fraction of a
    # caller's FPs within W bp of a breakpoint to the genomic fraction that close.
    W = 150
    L = contig_length(a.data_dir, a.ref_str, a.region)
    truth_all = np.array(sorted(load_keyset(
        os.path.join(sim, a.ref_str + "-truth.norm.vcf.gz"), exclude_sv=True)))
    genomic_frac = min(1.0, bp_pos.size * 2.0 * W / L)  # sparse breakpoints -> ~no overlap
    fp_rows = []
    for caller in a.callers:
        cvcf = os.path.join(a.data_dir, "iter", caller, a.cov + "X", a.region, a.div,
                            a.het, a.seed, "iter1", a.ref_str + "_calls.norm.vcf.gz")
        if not os.path.isfile(cvcf):
            continue
        cpos, ckeys = called_positions(cvcf)
        is_fp = ~np.isin(ckeys, truth_all)               # called but not a truth site
        fp_pos = cpos[is_fp]
        if fp_pos.size == 0:
            continue
        d, _ = nearest_breakpoint(fp_pos, bp_pos, bp_typ)
        near_frac = float(np.mean(d <= W))
        fp_rows.append((caller, int(fp_pos.size), near_frac, genomic_frac,
                        near_frac / genomic_frac if genomic_frac > 0 else float("nan")))

    fp_out = a.out.replace(".csv", "-fp.csv")
    with open(fp_out, "w") as out:
        out.write("caller,n_fp,frac_fp_within_%dbp,genomic_frac_within_%dbp,enrichment\n" % (W, W))
        for r in fp_rows:
            out.write("%s,%d,%.6g,%.6g,%.4g\n" % r)

    with open(a.out + ".dict.txt", "w") as d:
        d.write("Data dictionary for %s\n\n" % a.out)
        d.write("s07 SV-confounder: small-variant recovery vs distance from nearest SV\n")
        d.write("breakpoint. One row per (caller, iteration, vtype, distance bin).\n\n")
        d.write("Fields:\n")
        d.write("  caller, iteration : condition (SNP div background = %s; SVs uncalled)\n" % a.div)
        d.write("  vtype     : snp | indel (GT 1|1 divergence truth sites)\n")
        d.write("  bptype    : type of the NEAREST breakpoint -- any | del | ins (does a\n")
        d.write("              deletion's coverage gap hurt more than an insertion's anchor?)\n")
        d.write("  dist_lo/hi: distance-to-nearest-SV-breakpoint bin edges (bp); hi=-1 is the\n")
        d.write("              far/SV-free baseline bin (>= %g bp)\n" % DIST_EDGES[-2])
        d.write("  dist_mid  : bin midpoint (plotting x); = dist_lo for the baseline bin\n")
        d.write("  n_sites, n_recovered, recall : truth sites / recovered / recall in the bin\n\n")
        d.write("SV breakpoints: deletions contribute both ends, insertions their anchor.\n")
        d.write("Recovery should crater near 0 and rise to the baseline far away.\n")

    sys.stderr.write("# %d breakpoints, %d truth sites (%d snp, %d indel); %d rows\n"
                     % (bps.size, len(keys), int((vtypes == "snp").sum()),
                        int((vtypes == "indel").sum()), len(rows)))


if __name__ == "__main__":
    main()
