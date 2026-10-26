#!/usr/bin/env python3
"""
annotate_local_divergence.py -- attach a MEASURED local divergence to every
fixed-divergence (1|1) truth SNP, for the s05 heterogeneous-divergence analysis.

WHY THIS EXISTS
---------------
The s05 truth VCF places divergent sites with an autocorrelated field, so each
site sits in a local neighborhood that is more or less diverged than the genome
mean. Every downstream per-site analysis (recall, read-fate decomposition,
phantom rate, assemblable neighborhood) is binned against this LOCAL divergence
rather than the single global D. This script computes, for each divergent truth
site, the local divergent-SNP density in a window around it -- the divergence a
read/fragment covering the site actually experiences.

It reads the divergent (1|1) SNP positions straight from the truth VCF (records
with ID "div_snp_*"), so it works for both homogeneous (sigma=0) and field
(sigma>0) truth sets; a homogeneous set simply yields a near-flat local D.

DEFINITION
----------
  local_d_<W> = (number of divergent truth SNPs within +/- W bp of the site,
                 including the site itself) / (2 * W)
  Two scales are emitted by default: W=150 (read footprint) and W=500 (fragment
  footprint). This is divergent-site density per bp over the window; reads are
  drawn from A/C/G/T sequence, so the small chance a window overlaps an N gap is
  ignored (documented rather than silently corrected).

OUTPUT
------
  A TSV keyed by position (CHROM, POS 1-based) with one local_d_<W> column per
  window, plus a sidecar data dictionary. Downstream tables join on POS.
"""

import sys
import argparse
import gzip
import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0])  # allow import of sibling module
import divergence_field as dfield
import provenance                       # shared provenance block for the data dictionaries


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--truth-vcf", required=True,
                   help="Truth VCF (mason format, .vcf or .vcf.gz) from make_truth_vcf.py.")
    p.add_argument("--windows", default="150,500",
                   help="Comma-separated half-window sizes in bp (default 150,500).")
    p.add_argument("--out", required=True, help="Output TSV path.")
    return p.parse_args()


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def load_div_sites(vcf_path):
    """Return (chrom, sorted 0-based div-SNP positions, sorted 0-based div-indel
    positions). Indels are the div_indel_* records (empty array if none)."""
    chrom = None
    snp, ind = [], []
    with opener(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            vid = f[2]
            if vid.startswith("div_snp"):
                chrom = f[0]; snp.append(int(f[1]) - 1)      # 1-based -> 0-based
            elif vid.startswith("div_indel"):
                chrom = f[0]; ind.append(int(f[1]) - 1)
    if chrom is None:
        sys.exit("ERROR: no div_snp_*/div_indel_* records found in %s" % vcf_path)
    return (chrom,
            np.sort(np.array(snp, dtype=np.int64)),
            np.sort(np.array(ind, dtype=np.int64)))


def main():
    args = parse_args()
    windows = [int(w) for w in args.windows.split(",")]
    chrom, snp_pos, ind_pos = load_div_sites(args.truth_vcf)

    # Local divergence = divergent-SNP density (the field) at each window scale.
    # SNP sites: density around each div-SNP (self-inclusive). Indel sites (if any):
    # density of div-SNPs around the indel position (co-clustering shows up as an
    # elevated local SNP density). Same axis as the SNP-only s05 analysis.
    snp_local = {w: dfield.windowed_local_rate(snp_pos, w) for w in windows}
    ind_local = {w: dfield.windowed_density_at(ind_pos, snp_pos, w) for w in windows} \
        if ind_pos.size else {w: np.array([]) for w in windows}

    with open(args.out, "w") as out:
        cols = ["CHROM", "POS", "TYPE"] + ["local_d_%d" % w for w in windows]
        out.write("\t".join(cols) + "\n")
        for i in range(snp_pos.size):
            row = [chrom, str(snp_pos[i] + 1), "snp"] + ["%.8g" % snp_local[w][i] for w in windows]
            out.write("\t".join(row) + "\n")
        for i in range(ind_pos.size):
            row = [chrom, str(ind_pos[i] + 1), "indel"] + ["%.8g" % ind_local[w][i] for w in windows]
            out.write("\t".join(row) + "\n")

    with open(args.out + ".dict.txt", "w") as d:
        d.write(provenance.header(args.out))
        d.write("One row per FIXED-divergence (1|1) truth SNP, with its measured local\n")
        d.write("divergence at each window scale. Downstream s05 tables join on POS.\n\n")
        d.write("Fields:\n")
        d.write("  CHROM      : contig name (matches the truth VCF / mapping reference)\n")
        d.write("  POS        : 1-based reference position of the divergent site\n")
        d.write("  TYPE       : snp (div_snp_*, 1|1) or indel (div_indel_*, 1|1)\n")
        for w in windows:
            d.write("  local_d_%-4d: divergent-SNP density (sites/bp) within +/- %d bp of POS\n"
                    % (w, w))
            d.write("               (the field). SNP rows are self-inclusive; indel rows\n")
            d.write("               measure the div-SNP density around the indel = count/(2*%d)\n" % w)
        d.write("\nWindows chosen: read footprint (150 bp) and fragment footprint (500 bp).\n")

    # per-scale summary for the log
    for w in windows:
        for tag, lr in (("snp", snp_local[w]), ("indel", ind_local[w])):
            if lr.size == 0:
                continue
            sys.stderr.write("# %s local_d_%d: n=%d mean=%.6g p05=%.6g median=%.6g p95=%.6g max=%.6g\n"
                             % (tag, w, lr.size, lr.mean(),
                                np.quantile(lr, 0.05), np.median(lr),
                                np.quantile(lr, 0.95), lr.max()))


if __name__ == "__main__":
    main()
