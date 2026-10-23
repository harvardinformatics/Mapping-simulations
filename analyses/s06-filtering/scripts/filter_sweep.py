#!/usr/bin/env python3
"""
filter_sweep.py -- can you filter your way out of reference bias? (s06)

Sweeps standard hard-filter thresholds (QUAL, per-sample depth DP, mapping quality
MQ) over the RAW single-shot caller VCFs (s01 iter1 `calls.norm`, SNP-only), scoring
recall of true DIVERGENCE and precision at each threshold. No new mapping/calling --
the caller output already exists; each VCF is read once and all thresholds are swept
in memory.

The point: a filter only ever REMOVES calls, so divergence recall can only go DOWN as
you filter -- its ceiling is the unfiltered value, which is already < 1 because of
reference bias (the lost variants are absent from the VCF, not low-quality entries in
it). Filtering buys PRECISION (removes false positives, which matters most for
freebayes at high divergence) but cannot recover the reference-bias (false-negative)
gap. That is the mechanism from s05's mapping-vs-calling split: filtering acts on
emitted calls; the losses are un-emitted.

DEFINITIONS (per filtered call set):
  * truth_div      = truth 1|1 (fixed-difference) SNP sites (POS,REF,ALT).
  * truth_all      = all truth SNP sites (1|1 and 0|1) -- used to flag false positives.
  * div_tp         = calls that are HOM-ALT and match a truth_div site.
  * recall         = div_tp / |truth_div|      (the reference-bias / divergence-recovery axis)
  * n_called       = SNP calls passing the filter
  * tp_any         = calls matching ANY truth SNP site (correct locus+allele)
  * precision      = tp_any / n_called ;  fp = n_called - tp_any

A call's (POS,REF,ALT) is compared after identical bcftools normalization on both
sides, so exact key match is valid.
"""

import os
import sys
import argparse
import subprocess
import numpy as np

QUAL_T = [0, 10, 20, 30, 50, 100]
DP_T   = [0, 2, 4, 6, 8, 10]
MQ_T   = [0, 20, 30, 40, 50, 60]
MQ_TAG = {"gatk": "MQ", "bcftools": "MQ", "freebayes": "MQM"}
BASE = {"A": 0, "C": 1, "G": 2, "T": 3}


def snp_key(pos, ref, alt):
    """Encode a biallelic-SNP (POS,REF,ALT) as one int64 for fast np.isin membership;
    None if REF/ALT are not single A/C/G/T bases."""
    r, a = BASE.get(ref), BASE.get(alt)
    if r is None or a is None:
        return None
    return (int(pos) << 4) | (r << 2) | a


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data-dir", required=True, help="s01 data root (has sim/ and iter/)")
    p.add_argument("--ref-str", default="mm39")
    p.add_argument("--region", default="19")
    p.add_argument("--cov", default="30")
    p.add_argument("--het", default="0.005")
    p.add_argument("--seed", default="42")
    p.add_argument("--iter", default="1", help="which iteration's raw calls (default 1 = single-shot)")
    p.add_argument("--callers", nargs="+", default=["gatk", "bcftools", "freebayes"])
    p.add_argument("--divs", nargs="+", required=True)
    p.add_argument("--out", required=True)
    return p.parse_args()


def bcf_query(path, fmt):
    """Stream `bcftools query -f fmt` over a VCF, yielding split lines."""
    p = subprocess.Popen(["bcftools", "query", "-f", fmt, path],
                         stdout=subprocess.PIPE, text=True)
    for line in p.stdout:
        yield line.rstrip("\n").split("\t")
    p.wait()
    if p.returncode != 0:
        sys.exit("ERROR: bcftools query failed on %s" % path)


def load_truth(path):
    """Return {'snp': (tdiv, tall) int64 keys, 'indel': (tdiv, tall) string keys}.
    SNP keys are int64-encoded (fast np.isin over millions); indels are ~8x fewer, so
    'POS:REF:ALT' string keys are fine."""
    snp_div, snp_all, ind_div, ind_all = [], [], [], []
    for f in bcf_query(path, "%POS\t%REF\t%ALT\t[%GT]\n"):
        pos, ref, alt, gt = f
        ref, alt = ref.upper(), alt.upper()
        div = gt in ("1|1", "1/1")
        if len(ref) == 1 and len(alt) == 1:
            k = snp_key(pos, ref, alt)
            if k is None:
                continue
            snp_all.append(k)
            if div:
                snp_div.append(k)
        else:
            k = "%s:%s:%s" % (pos, ref, alt)
            ind_all.append(k)
            if div:
                ind_div.append(k)
    return {
        "snp": (np.unique(np.array(snp_div, dtype=np.int64)),
                np.unique(np.array(snp_all, dtype=np.int64))),
        "indel": (np.unique(np.array(ind_div)), np.unique(np.array(ind_all))),
    }


def load_calls(path, mq_tag):
    """Return {'snp': (...), 'indel': (...)} where each value is
    (keys, gts, quals, dps, mqs) -- SNP keys int64, indel keys 'POS:REF:ALT' strings."""
    buf = {"snp": ([], [], [], [], []), "indel": ([], [], [], [], [])}
    fmt = "%%POS\t%%REF\t%%ALT\t[%%GT]\t%%QUAL\t[%%DP]\t%%INFO/%s\n" % mq_tag
    for f in bcf_query(path, fmt):
        pos, ref, alt, gt, qual, dp, mq = f
        ref, alt = ref.upper(), alt.upper()
        q = float(qual) if qual not in (".", "") else 0.0
        d = int(dp) if dp not in (".", "") else 0
        m = float(mq) if mq not in (".", "") else 0.0
        if len(ref) == 1 and len(alt) == 1:
            k = snp_key(pos, ref, alt)
            if k is None:
                continue
            vt = "snp"
        else:
            k = "%s:%s:%s" % (pos, ref, alt)
            vt = "indel"
        b = buf[vt]
        b[0].append(k); b[1].append(gt); b[2].append(q); b[3].append(d); b[4].append(m)
    out = {}
    out["snp"] = (np.array(buf["snp"][0], dtype=np.int64), np.array(buf["snp"][1]),
                  np.array(buf["snp"][2]), np.array(buf["snp"][3], dtype=float),
                  np.array(buf["snp"][4]))
    out["indel"] = (np.array(buf["indel"][0]), np.array(buf["indel"][1]),
                    np.array(buf["indel"][2]), np.array(buf["indel"][3], dtype=float),
                    np.array(buf["indel"][4]))
    return out


def main():
    a = parse_args()
    rows = []
    for caller in a.callers:
        mq_tag = MQ_TAG.get(caller, "MQ")
        for div in a.divs:
            base = os.path.join(a.data_dir)
            truth = os.path.join(base, "sim", a.cov + "X", a.region, div, a.het, a.seed,
                                 a.ref_str + "-truth.norm.vcf.gz")
            calls = os.path.join(base, "iter", caller, a.cov + "X", a.region, div, a.het,
                                 a.seed, "iter" + a.iter, a.ref_str + "_calls.norm.vcf.gz")
            if not (os.path.isfile(truth) and os.path.isfile(calls)):
                sys.stderr.write("WARN: missing files for %s div=%s (skipped)\n" % (caller, div))
                continue

            truth_sets = load_truth(truth)
            call_sets = load_calls(calls, mq_tag)

            for vt in ("snp", "indel"):
                tdiv, tall = truth_sets[vt]
                keys, gts, quals, dps, mqs = call_sets[vt]
                n_truth_div = tdiv.size
                if n_truth_div == 0 or keys.size == 0:
                    continue                       # e.g. no indels in an s01 (SNP-only) run

                is_hom = np.isin(gts, ("1|1", "1/1"))
                is_div_tp = np.isin(keys, tdiv, assume_unique=False) & is_hom
                is_any_tp = np.isin(keys, tall, assume_unique=False)

                for axis, field, thresholds in (("QUAL", quals, QUAL_T),
                                                ("DP", dps, DP_T),
                                                ("MQ", mqs, MQ_T)):
                    for t in thresholds:
                        mask = field >= t
                        n_called = int(mask.sum())
                        div_tp = int((is_div_tp & mask).sum())
                        tp_any = int((is_any_tp & mask).sum())
                        fp = n_called - tp_any
                        recall = div_tp / n_truth_div if n_truth_div else 0.0
                        precision = tp_any / n_called if n_called else 0.0
                        rows.append((caller, div, vt, axis, t, n_truth_div, div_tp,
                                     n_called, tp_any, fp, recall, precision))

                sys.stderr.write("# %s div=%s %s: n_truth_div=%d n_calls=%d\n"
                                 % (caller, div, vt, n_truth_div, keys.size))

    with open(a.out, "w") as out:
        out.write("caller,divergence,vtype,axis,threshold,n_truth_div,div_tp,"
                  "n_called,tp_any,fp,recall,precision\n")
        for r in rows:
            out.write("%s,%s,%s,%s,%g,%d,%d,%d,%d,%d,%.6g,%.6g\n" % r)

    with open(a.out + ".dict.txt", "w") as d:
        d.write("Data dictionary for %s\n\n" % a.out)
        d.write("s06 filtering sweep: recall of true divergence + precision vs hard-filter\n")
        d.write("threshold, over the RAW single-shot (iter %s) SNP caller VCFs.\n\n" % a.iter)
        d.write("Fields:\n")
        d.write("  caller, divergence : condition (homogeneous; s01 = SNP-only, s02 = SNP+indel)\n")
        d.write("  vtype     : snp | indel (indel rows only present when the run has indels, i.e. s02)\n")
        d.write("  axis      : filter axis -- QUAL | DP (per-sample depth) | MQ (mapping quality;\n")
        d.write("              INFO/MQ for gatk & bcftools, INFO/MQM for freebayes)\n")
        d.write("  threshold : keep calls with axis >= threshold\n")
        d.write("  n_truth_div : truth 1|1 SNP sites\n")
        d.write("  div_tp    : hom-alt calls matching a truth 1|1 site (passing the filter)\n")
        d.write("  n_called  : SNP calls passing the filter\n")
        d.write("  tp_any    : calls matching ANY truth SNP site (correct locus+allele)\n")
        d.write("  fp        : n_called - tp_any (false positives passing the filter)\n")
        d.write("  recall    : div_tp / n_truth_div  (divergence recovery = the reference-bias axis)\n")
        d.write("  precision : tp_any / n_called\n\n")
        d.write("threshold=0 (or the axis minimum) = UNFILTERED. recall at threshold 0 is the\n")
        d.write("ceiling; any filter can only lower it. Reference bias is the gap 1 - recall(0).\n")

    sys.stderr.write("# wrote %d rows to %s\n" % (len(rows), a.out))


if __name__ == "__main__":
    main()
