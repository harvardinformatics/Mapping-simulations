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

import provenance                       # shared provenance block for the data dictionaries


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


# Written into every dictionary of this tool's outputs (and of tables derived from them).
ARTIFACT_NOTE = ("KNOWN ARTIFACT (found 2026-10-22): most sites flagged phantom here are NOT phantoms. "
                 "mason -oa writes dense SNP clusters in the truth BAM as gapped 1D..1I CIGAR pairs, so "
                 "reading bases through those CIGARs puts them at shifted positions and real alt sites "
                 "look all-reference; bwa-mapped reads carry the alt at most of them. Do not use these "
                 "counts or the effective truth until the scan compares read SEQUENCES to the truth "
                 "haplotype instead.")


def write_dictionaries(a, notes):
    """Dictionaries for <prefix>-phantom.tsv, -effective.vcf.gz and -summary.csv. With --dicts-only
    they document existing outputs (backfill) under the ORIGINAL command."""
    def hdr(table):
        notes_all = [ARTIFACT_NOTE] + list(notes)
        if a.dicts_only:
            return provenance.header(table, command=a.dicts_only, backfill=True, notes=notes_all)
        return provenance.header(table, notes=notes_all)

    method = ("Method: one pass over the mason -oa truth BAM (%s); for every biallelic truth SNP (1|1\n"
              "and 0|1) the base each primary read presents at the SNP is tallied as alt / ref / other.\n"
              "A SNP with truth-read depth >= %d and alt fraction < %.2f is flagged phantom; depth < %d\n"
              "leaves it unclassified (kept as truth). Truth VCF: %s.\n\n"
              % (a.truth_bam, a.min_depth, a.alt_frac, a.min_depth, a.truth_vcf))
    ph = a.out_prefix + "-phantom.tsv"
    with open(ph + ".dict.txt", "w") as d:
        d.write(hdr(ph) + method)
        d.write("One row per flagged phantom SNP (tab-separated).\n\nFields:\n"
                "  pos       : 1-based reference position\n"
                "  ref, alt  : truth alleles\n"
                "  zyg       : div (hom-div 1|1) | het (0|1)\n"
                "  depth     : truth-BAM reads presenting a base at pos\n"
                "  alt_count : of those, presenting the alt base\n"
                "  alt_frac  : alt_count / depth\n")
    eff = a.out_prefix + "-effective.vcf.gz"
    with open(eff + ".dict.txt", "w") as d:
        d.write(hdr(eff) + method)
        d.write("The truth VCF with the flagged phantom SNP records removed; every other record (incl.\n"
                "indels and multiallelics) passed through unchanged. bgzipped + tabix-indexed.\n")
    sm = a.out_prefix + "-summary.csv"
    with open(sm + ".dict.txt", "w") as d:
        d.write(hdr(sm) + method)
        d.write("One row (this region / divergence / seed).\n\nFields:\n"
                "  region, divergence, seed : condition labels\n"
                "  n_truth_snp       : biallelic truth SNPs (div + het)\n"
                "  n_truth_div, n_truth_het : of those, hom-div (1|1) / het (0|1)\n"
                "  n_phantom_div, n_phantom_het : flagged phantom, by zygosity\n"
                "  phantom_div_rate  : n_phantom_div / n_truth_div\n"
                "  phantom_het_rate  : n_phantom_het / n_truth_het\n"
                "  n_lowcov          : SNPs with truth-read depth < min depth (unclassified)\n"
                "  mean_alt_frac     : mean alt fraction over the classified SNPs\n")


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
    ap.add_argument("--dicts-only", metavar="ORIGINAL_COMMAND",
                    help="write ONLY the data dictionaries for existing outputs (same arguments as the "
                         "run that made them), recording ORIGINAL_COMMAND as provenance")
    ap.add_argument("--dict-note", action="append", default=[],
                    help="extra provenance note for the dictionaries (repeatable)")
    a = ap.parse_args()

    # --dicts-only: document existing outputs without recomputing them
    if a.dicts_only:
        write_dictionaries(a, a.dict_note)
        return

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
    write_dictionaries(a, a.dict_note)
    print(f"div phantom {n_phantom_div}/{n_div} ({100 * n_phantom_div / max(1, n_div):.2f}%), "
          f"het phantom {n_phantom_het}/{n_het} ({100 * n_phantom_het / max(1, n_het):.2f}%); "
          f"wrote {a.out_prefix}-effective.vcf.gz")


if __name__ == "__main__":
    main()
