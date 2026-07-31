#!/usr/bin/env python3
"""
compare_callers.py -- does a pileup-based caller (bcftools) recover the divergence
sites that GATK HaplotypeCaller lost, broken down by the site-linkage outcome class?

Recovery = the caller reports homozygous-alt with the matching truth allele.

NOTE on the truth VCF: make_truth_vcf.py writes mason format (FORMAT="."), whose
genotype is NOT readable by pysam/bcftools. We therefore parse the truth as plain
text (as read_site_linkage.py does), not with pysam.
"""

import sys
import gzip
import csv
import pysam
from collections import defaultdict


def load_truth_div(path):
    """pos -> alt base for homozygous-alt (1|1) SNP truth sites (mason plain-text)."""
    truth = {}
    op = gzip.open(path, "rt") if path.endswith(".gz") else open(path)
    with op as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[9] == "1|1" and len(f[3]) == 1 and len(f[4]) == 1:
                truth[int(f[1])] = f[4]
    return truth


def load_homalt_snps(path):
    """pos -> alt base for homozygous-alt SNP calls in a standard VCF (via pysam)."""
    calls = {}
    vf = pysam.VariantFile(path)
    for rec in vf:
        if len(rec.ref) != 1:
            continue
        gt = rec.samples[0].get("GT", (None, None))
        idx = [a for a in gt if a is not None]
        if not idx or not all(i != 0 for i in idx) or len(set(idx)) != 1:
            continue
        alt = rec.alleles[idx[0]]
        if len(alt) == 1:
            calls[rec.pos] = alt
    vf.close()
    return calls


def main():
    bcf_path, truth_path, link_path, divergence, out_path = sys.argv[1:6]
    truth = load_truth_div(truth_path)
    bcf = load_homalt_snps(bcf_path)
    sys.stderr.write("# truth div SNPs=%d  bcftools hom-alt SNPs=%d\n" % (len(truth), len(bcf)))

    groups = defaultdict(list)
    for r in csv.DictReader(open(link_path)):
        if r["divergence"] != divergence:
            continue
        vaf = float(r["pileup_alt_frac_hq"])
        cls = ("recovered by GATK (tp)" if r["outcome"] == "tp"
               else "GATK-lost, evidence present" if vaf >= 0.5
               else "GATK-lost, evidence degraded")
        groups[cls].append(int(r["pos"]))

    lines = ["Caller comparison at div=%s (SNP-only condition)." % divergence,
             "bcftools (pileup) vs GATK HaplotypeCaller (assembly).",
             "recovered = caller emits hom-alt with the matching truth allele.\n"]
    for cls in ["recovered by GATK (tp)", "GATK-lost, evidence degraded",
                "GATK-lost, evidence present"]:
        ps = groups[cls]
        if not ps:
            continue
        rec = sum(1 for p in ps if bcf.get(p) == truth.get(p))
        lines.append("%-32s n=%6d  bcftools recovers %6d  (%.1f%%)"
                     % (cls, len(ps), rec, 100 * rec / len(ps)))
    open(out_path, "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
