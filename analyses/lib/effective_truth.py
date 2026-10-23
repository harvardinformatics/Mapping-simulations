#!/usr/bin/env python3
"""
effective_truth.py -- separate SIMULATED truth variants from PHANTOM ones.

A "phantom" truth site is a variant that is in the truth VCF but was never materialized into the
reads by mason (typically because it sits in a dense variant cluster mason drops silently). At such a
site the mason `-oa` truth BAM -- which places every read at its TRUE origin, so mapping cannot be
blamed -- carries the REFERENCE allele instead of the alt. Phantoms are UNRECOVERABLE by any caller or
iteration (the alt is simply not in the data), so leaving them in the recall denominator understates
recall and mislabels them as reference-bias loss. This tool scans the truth BAM once, measures the
alt-allele fraction among truth reads at every truth SNP, and flags those with the alt essentially
absent. It writes: the phantom list, an EFFECTIVE truth VCF (truth minus phantoms) for re-scoring, and
a one-row summary. SNPs only (biallelic); indels and multiallelics are passed through to the effective
VCF unchanged (their read-level presence needs CIGAR inspection -- not handled here).

Usage:
  python effective_truth.py --truth-bam mm39-truth.oa.bam --truth-vcf mm39-truth.norm.vcf.gz \
      --region 19 --divergence 0.10 --seed 42 --out-prefix results/eff/eff-0.10
"""
import argparse
import bisect
import gzip
import subprocess

import pysam


def load_snps(vcf_path):
    """Return sorted [(pos0, ref, alt, is_het)] for biallelic SNP truth sites (both 1|1 and 0|1)."""
    out = []
    op = gzip.open(vcf_path, "rt") if vcf_path.endswith(".gz") else open(vcf_path)
    for line in op:
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        ref, alt, gt = c[3], c[4], c[9].split(":")[0]
        if len(ref) == 1 and len(alt) == 1:
            is_het = gt in ("0|1", "1|0", "0/1", "1/0")
            out.append((int(c[1]) - 1, ref.upper(), alt.upper(), is_het))
    out.sort()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth-bam", required=True, help="mason -oa truth BAM (reads at TRUE positions)")
    ap.add_argument("--truth-vcf", required=True, help="normalized truth VCF (bgzipped)")
    ap.add_argument("--region", required=True)
    ap.add_argument("--divergence", required=True)
    ap.add_argument("--seed", required=True)
    ap.add_argument("--min-depth", type=int, default=5,
                    help="truth-read depth below which a site is left UNCLASSIFIED (kept as truth)")
    ap.add_argument("--alt-frac", type=float, default=0.10,
                    help="alt-allele fraction below which the alt is 'absent' -> phantom (default 0.10)")
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()

    snps = load_snps(a.truth_vcf)
    poslist = [p for p, _, _, _ in snps]
    altc = [0] * len(snps)
    refc = [0] * len(snps)
    othc = [0] * len(snps)

    # single linear pass over the truth BAM, tallying the base each read presents at each covered SNP
    with pysam.AlignmentFile(a.truth_bam, "rb") as tf:
        for r in tf:
            if r.is_secondary or r.is_supplementary or r.is_unmapped:
                continue
            lo = bisect.bisect_left(poslist, r.reference_start)
            hi = bisect.bisect_left(poslist, r.reference_end)
            if hi <= lo:
                continue
            aln = {rp: qp for qp, rp in r.get_aligned_pairs(matches_only=True)}
            seq = r.query_sequence
            for j in range(lo, hi):
                qp = aln.get(poslist[j])
                if qp is None:
                    continue
                b = seq[qp].upper()
                if b == snps[j][2]:
                    altc[j] += 1
                elif b == snps[j][1]:
                    refc[j] += 1
                else:
                    othc[j] += 1

    # classify: phantom = enough depth to judge AND alt essentially absent
    phantom = set()                                   # 0-based positions dropped from effective truth
    n_phantom_div = n_phantom_het = n_lowcov = 0
    alt_frac_sum = n_scored = 0.0
    with open(a.out_prefix + "-phantom.tsv", "w") as pf:
        pf.write("pos\tref\talt\tzyg\tdepth\talt_count\talt_frac\n")
        for j, (p, ref, alt, is_het) in enumerate(snps):
            depth = altc[j] + refc[j] + othc[j]
            if depth < a.min_depth:
                n_lowcov += 1
                continue
            af = altc[j] / depth
            alt_frac_sum += af
            n_scored += 1
            if af < a.alt_frac:
                phantom.add(p)
                zyg = "het" if is_het else "div"
                if is_het:
                    n_phantom_het += 1
                else:
                    n_phantom_div += 1
                pf.write(f"{p + 1}\t{ref}\t{alt}\t{zyg}\t{depth}\t{altc[j]}\t{af:.3f}\n")

    # effective truth VCF = truth minus phantom SNP positions (non-SNP lines passed through)
    eff_plain = a.out_prefix + "-effective.vcf"
    op = gzip.open(a.truth_vcf, "rt") if a.truth_vcf.endswith(".gz") else open(a.truth_vcf)
    with open(eff_plain, "w") as out:
        for line in op:
            if line.startswith("#"):
                out.write(line)
                continue
            c = line.split("\t")
            pos0 = int(c[1]) - 1
            if len(c[3]) == 1 and len(c[4]) == 1 and pos0 in phantom:
                continue                              # drop the phantom SNP
            out.write(line)
    subprocess.run(["bgzip", "-f", eff_plain], check=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", eff_plain + ".gz"], check=True)

    n_div = sum(1 for _, _, _, h in snps if not h)
    n_het = sum(1 for _, _, _, h in snps if h)
    with open(a.out_prefix + "-summary.csv", "w") as s:
        s.write("region,divergence,seed,n_truth_snp,n_truth_div,n_truth_het,"
                "n_phantom_div,n_phantom_het,phantom_div_rate,phantom_het_rate,"
                "n_lowcov,mean_alt_frac\n")
        s.write(f"{a.region},{a.divergence},{a.seed},{len(snps)},{n_div},{n_het},"
                f"{n_phantom_div},{n_phantom_het},"
                f"{n_phantom_div / n_div if n_div else 0:.5f},"
                f"{n_phantom_het / n_het if n_het else 0:.5f},"
                f"{n_lowcov},{alt_frac_sum / n_scored if n_scored else 0:.4f}\n")
    print(f"div phantom {n_phantom_div}/{n_div} ({100 * n_phantom_div / max(1, n_div):.2f}%), "
          f"het phantom {n_phantom_het}/{n_het} ({100 * n_phantom_het / max(1, n_het):.2f}%); "
          f"wrote {a.out_prefix}-effective.vcf.gz")


if __name__ == "__main__":
    main()
