#!/usr/bin/env python3
"""
compare_vcfs.py -- variant-calling accuracy for the v2 mason pipeline (SNPs + indels).

Compares the GATK-called VCF against the controlled truth VCF from make_truth_vcf.py.
Both are expected to be NORMALIZED (bcftools norm -f ref, multiallelics split) so that
indel representations are left-aligned and comparable; the workflow does this before
calling this script. Both are read with pysam.VariantFile.

Truth sites are homozygous-alt (GT 1|1, "divergence") or heterozygous (0|1). Each is a
SNP (single-base REF and ALT) or a small indel. We classify every truth site by how it
was called, and report SNPs and indels separately.

SNP divergence sites (1|1) -> exactly one of:
  div_tp             called hom-alt with MATCHING allele (recovered)
  div_het            called heterozygous
  div_ref_or_absent  called hom-ref OR not called            <-- reference bias
  div_other          hom-alt wrong allele / multiallelic edge case
SNP heterozygous sites (0|1) -> het_tp / het_hom_alt / het_ref_or_absent / het_other.

Indel divergence sites (1|1) -> indel_div_tp / indel_div_ref_or_absent / indel_div_other.
Indel heterozygous sites (0|1) -> indel_het_tp / indel_het_ref_or_absent / indel_het_other.
(A called indel matches a truth indel when they share (POS, REF, ALT) after normalization.)

False positives: called variants not matching any truth site -> fp (SNP) / fp_indel.

The SNP output column names are unchanged from the SNP-only version so existing
downstream analyses keep working; indel and seed columns are appended.
"""

import sys
import argparse
import pysam


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--truth-vcf", required=True, help="normalized truth VCF (bgzipped)")
    p.add_argument("--called-vcf", required=True, help="normalized GATK VCF (bgzipped)")
    p.add_argument("--region", required=True)
    p.add_argument("--coverage", required=True)
    p.add_argument("--divergence", required=True)
    p.add_argument("--heterozygosity", required=True)
    p.add_argument("--seed", required=True)
    p.add_argument("--caller", default="gatk", help="variant caller label for the summary row")
    p.add_argument("--iteration", default="1", help="mapping iteration label (1 = single-shot)")
    p.add_argument("--extra", nargs="*", default=[], metavar="NAME=VALUE",
                   help="optional extra identifier columns appended to the summary row "
                        "(e.g. --extra field=L20000-S1.0). Omitting it leaves the header "
                        "unchanged, so other sections' output is unaffected.")
    p.add_argument("--out-prefix", required=True)
    return p.parse_args()


def genotype_category(rec):
    """Return (category, alt_allele) for a record: hom_ref / het / hom_alt / other.
    alt_allele is the relevant ALT string (any length) for het/hom_alt, else None."""
    gt = rec.samples[0].get("GT", (None, None))
    if gt is None or all(a is None for a in gt):
        return "hom_ref", None
    idx = [a for a in gt if a is not None]
    if all(a == 0 for a in idx):
        return "hom_ref", None
    if len(set(idx)) == 1 and idx[0] != 0:            # e.g. 1/1
        return "hom_alt", rec.alleles[idx[0]]
    nonzero = [a for a in idx if a != 0]
    if 0 in idx and len(set(nonzero)) == 1:           # e.g. 0/1
        return "het", rec.alleles[nonzero[0]]
    return "other", None


def is_snp(ref, alt):
    return len(ref) == 1 and len(alt) == 1


def load_truth(path):
    """(pos, ref, alt) -> (is_indel, is_div). Keyed by the full allele so that two
    variants left-aligned onto the same POS by normalization are both kept."""
    truth = {}
    vf = pysam.VariantFile(path)
    for rec in vf:
        alt = rec.alts[0]
        gt = rec.samples[0].get("GT", (None, None))
        is_div = tuple(gt) == (1, 1)                  # 1|1 divergence vs 0|1 het
        truth[(rec.pos, rec.ref, alt)] = (not is_snp(rec.ref, alt), is_div)
    vf.close()
    return truth


def load_called(path):
    """pos -> list of (ref, alt, category). Normalized/split, so usually one per pos."""
    called = {}
    vf = pysam.VariantFile(path)
    for rec in vf:
        cat, alt = genotype_category(rec)
        for a in rec.alts:                            # after -m- there is one ALT
            called.setdefault(rec.pos, []).append((rec.ref, a, cat))
    vf.close()
    return called


def matched_status(called, pos, ref, alt):
    """How was a truth allele called? One of:
       'hom_alt'/'het'/'hom_ref' (exact allele match at this POS with that genotype),
       'other_call' (a DIFFERENT variant allele was called here),
       'hom_ref'    (only a hom-ref record is present), or None (nothing called here)."""
    recs = called.get(pos)
    if not recs:
        return None
    for r_ref, r_alt, cat in recs:
        if r_ref == ref and r_alt == alt:
            return cat
    if any(cat in ("het", "hom_alt", "other") for _, _, cat in recs):
        return "other_call"
    return "hom_ref"


def main():
    args = parse_args()
    truth = load_truth(args.truth_vcf)
    called = load_called(args.called_vcf)

    c = {k: 0 for k in [
        "div_tp", "div_het", "div_ref_or_absent", "div_other",
        "het_tp", "het_hom_alt", "het_ref_or_absent", "het_other",
        "indel_div_tp", "indel_div_ref_or_absent", "indel_div_other",
        "indel_het_tp", "indel_het_ref_or_absent", "indel_het_other"]}
    n = {"snp_div": 0, "snp_het": 0, "indel_div": 0, "indel_het": 0}

    for (pos, t_ref, t_alt), (is_indel, is_div) in truth.items():
        st = matched_status(called, pos, t_ref, t_alt)
        absent_or_ref = st is None or st == "hom_ref"

        if not is_indel and is_div:                    # SNP divergence
            n["snp_div"] += 1
            if st == "hom_alt":
                c["div_tp"] += 1
            elif st == "het":
                c["div_het"] += 1
            elif absent_or_ref:
                c["div_ref_or_absent"] += 1
            else:
                c["div_other"] += 1
        elif not is_indel and not is_div:              # SNP heterozygous
            n["snp_het"] += 1
            if st == "het":
                c["het_tp"] += 1
            elif st == "hom_alt":
                c["het_hom_alt"] += 1
            elif absent_or_ref:
                c["het_ref_or_absent"] += 1
            else:
                c["het_other"] += 1
        elif is_indel and is_div:                      # indel divergence
            n["indel_div"] += 1
            if st == "hom_alt":
                c["indel_div_tp"] += 1
            elif absent_or_ref:
                c["indel_div_ref_or_absent"] += 1
            else:
                c["indel_div_other"] += 1
        else:                                          # indel heterozygous
            n["indel_het"] += 1
            if st == "het":
                c["indel_het_tp"] += 1
            elif absent_or_ref:
                c["indel_het_ref_or_absent"] += 1
            else:
                c["indel_het_other"] += 1

    # false positives + called counts: called variant alleles not matching any truth site
    truth_alleles = set(truth.keys())
    fp = fp_indel = n_called_snp = n_called_indel = 0
    for pos, recs in called.items():
        for r_ref, r_alt, cat in recs:
            if cat in ("het", "hom_alt", "other"):
                if is_snp(r_ref, r_alt):
                    n_called_snp += 1
                else:
                    n_called_indel += 1
                if (pos, r_ref, r_alt) not in truth_alleles:
                    if is_snp(r_ref, r_alt):
                        fp += 1
                    else:
                        fp_indel += 1

    # --- write summary CSV (SNP columns unchanged; indel/seed appended) ----
    cols = ["region", "coverage", "divergence", "heterozygosity", "seed", "caller", "iteration",
            "n_truth_div", "div_tp", "div_het", "div_ref_or_absent", "div_other",
            "n_truth_het", "het_tp", "het_hom_alt", "het_ref_or_absent", "het_other",
            "fp", "n_called_snp",
            "n_truth_indel_div", "indel_div_tp", "indel_div_ref_or_absent", "indel_div_other",
            "n_truth_indel_het", "indel_het_tp", "indel_het_ref_or_absent", "indel_het_other",
            "fp_indel", "n_called_indel"]
    vals = [args.region, args.coverage, args.divergence, args.heterozygosity, args.seed, args.caller, args.iteration,
            n["snp_div"], c["div_tp"], c["div_het"], c["div_ref_or_absent"], c["div_other"],
            n["snp_het"], c["het_tp"], c["het_hom_alt"], c["het_ref_or_absent"], c["het_other"],
            fp, n_called_snp,
            n["indel_div"], c["indel_div_tp"], c["indel_div_ref_or_absent"], c["indel_div_other"],
            n["indel_het"], c["indel_het_tp"], c["indel_het_ref_or_absent"], c["indel_het_other"],
            fp_indel, n_called_indel]

    # optional extra identifier columns (e.g. s05 field tag); appended so the base
    # header is unchanged when --extra is not given.
    for kv in args.extra:
        if "=" not in kv:
            sys.exit("ERROR: --extra expects NAME=VALUE (got %r)" % kv)
        name, value = kv.split("=", 1)
        cols.append(name)
        vals.append(value)

    csv_path = args.out_prefix + "-vcf-summary.csv"
    with open(csv_path, "w") as out:
        out.write(",".join(cols) + "\n")
        out.write(",".join(str(x) for x in vals) + "\n")

    with open(args.out_prefix + "-vcf-summary.dict.txt", "w") as d:
        d.write("Data dictionary for %s\n" % csv_path)
        d.write("Variant-calling accuracy vs controlled truth (SNPs + indels), normalized.\n")
        d.write("Divergence sites are hom-alt (1|1); het sites are 0|1. A called variant\n")
        d.write("matches a truth site by (POS, REF, ALT) after bcftools normalization.\n\n")
        d.write("region,coverage,divergence,heterozygosity,seed : simulation condition\n")
        d.write("caller           : variant caller (gatk | bcftools | freebayes)\n")
        d.write("iteration        : mapping iteration (1 = single-shot; 2..N = iterative)\n")
        d.write("n_truth_div       : number of truth SNP divergence (1|1) sites\n")
        d.write("div_tp            : SNP divergence called hom-alt, matching allele (recovered)\n")
        d.write("div_het           : SNP divergence called heterozygous\n")
        d.write("div_ref_or_absent : SNP divergence called hom-ref OR not called  <-- reference bias\n")
        d.write("div_other         : SNP divergence hom-alt wrong allele / edge case\n")
        d.write("n_truth_het       : number of truth SNP heterozygous (0|1) sites\n")
        d.write("het_tp/het_hom_alt/het_ref_or_absent/het_other : SNP het site call outcome\n")
        d.write("fp                : called SNP variants not matching any truth site (false positives)\n")
        d.write("n_called_snp      : total confident SNP variant calls\n")
        d.write("n_truth_indel_div : number of truth indel divergence (1|1) sites\n")
        d.write("indel_div_tp      : indel divergence recovered (hom-alt, matching REF/ALT)\n")
        d.write("indel_div_ref_or_absent : indel divergence called hom-ref OR not called\n")
        d.write("indel_div_other   : indel divergence other genotype (e.g. called het)\n")
        d.write("n_truth_indel_het : number of truth indel heterozygous (0|1) sites\n")
        d.write("indel_het_tp/indel_het_ref_or_absent/indel_het_other : indel het site outcome\n")
        d.write("fp_indel          : called indel variants not matching any truth site\n")
        d.write("n_called_indel    : total confident indel variant calls\n")

    print("wrote %s" % csv_path)
    for k in cols[5:]:
        print("  %-24s %s" % (k, dict(zip(cols, vals))[k]))


if __name__ == "__main__":
    main()
