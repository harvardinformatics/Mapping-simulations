#!/usr/bin/env python3
"""
make_truth_vcf.py -- generate a controlled "golden"/truth VCF for the mason2-based
reference-bias simulations (v2 pipeline).

WHY THIS EXISTS
---------------
mason_simulator applies a mason-format VCF to a reference to build the diploid sample
it sequences. mason_variator can produce such a VCF but assigns genotypes we cannot
control (in practice ~99% of sites differ on both haplotypes and only ~1% are clean
0|1 heterozygotes), so it gives no independent knob for heterozygosity. This script
writes a truth VCF with EXACTLY the genotype structure we want, reproducing the v1
model so v1 and v2 results are directly comparable.

VARIANT CLASSES WRITTEN
-----------------------
Four independently-controlled classes of SNPs and small indels:

  SNP divergence   (fraction D):  homozygous-alt SNP,   GT 1|1  (fixed difference)
  SNP heterozygous (fraction H):  heterozygous  SNP,    GT 0|1  (within-sample polymorphism)
  indel divergence (fraction Di): homozygous-alt indel, GT 1|1
  indel heterozygous (fraction Hi): heterozygous indel, GT 0|1

Small indels are insertions or deletions (50/50) of a size drawn uniformly from
[--min-indel-size, --max-indel-size]. Setting Di = Hi = 0 reproduces the SNP-only
behaviour of the earlier version.

DEFINITIONS (documented so results are interpretable and replicable)
--------------------------------------------------------------------
  * Contig length L = length of the single sequence in --reference.
  * Site counts: n_snp_div = round(D*L), n_snp_het = round(H*L),
    n_indel_div = round(Di*L), n_indel_het = round(Hi*L).
  * Divergence D = fraction of sites fixed for a non-reference SNP allele; the indel
    divergence Di is the analogous fraction for fixed indels.
  * All variant sites are DISJOINT: a chosen SNP occupies its single base; a deletion
    occupies its REF span; an insertion occupies its anchor base. New variants are
    rejected if they would overlap an already-placed variant or an N base.
  * Only A/C/G/T reference bases are eligible; N bases are skipped.

ALLELE MODEL
------------
  * SNP: REF = reference base; ALT = one of the three other bases, uniform (no ts/tv bias).
  * Deletion (size Lb) at anchor p: REF = ref[p .. p+Lb], ALT = ref[p] (Lb bases removed).
  * Insertion (size Lb) at anchor p: REF = ref[p], ALT = ref[p] + Lb random bases.

OUTPUT
------
  * A plain (uncompressed) VCF written to --out; the workflow bgzips + tabixes it. The
    header mirrors mason_variator's format so mason_simulator -iv accepts it, and all
    parameters + realized counts are recorded as ## header comments for provenance.
  * A sidecar data dictionary --out + ".dict.txt".
  * Realized counts/rates are printed to stderr for logging / validation.
"""

import sys
import argparse
import numpy as np
import pysam

BASES = "ACGT"
ALT_CHOICES = {"A": "CGT", "C": "AGT", "G": "ACT", "T": "ACG"}  # SNP alternates


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--reference", required=True,
                   help="FASTA with a SINGLE sequence (e.g. the extracted chr19).")
    p.add_argument("--divergence", type=float, required=True,
                   help="SNP divergence D: fraction of sites made homozygous-alt (1|1).")
    p.add_argument("--heterozygosity", type=float, required=True,
                   help="SNP heterozygosity H: fraction of sites made het (0|1).")
    p.add_argument("--indel-divergence", type=float, default=0.0,
                   help="Indel divergence Di: fraction of sites made homozygous-alt indels.")
    p.add_argument("--indel-heterozygosity", type=float, default=0.0,
                   help="Indel heterozygosity Hi: fraction of sites made het indels.")
    p.add_argument("--min-indel-size", type=int, default=1)
    p.add_argument("--max-indel-size", type=int, default=6)
    p.add_argument("--seed", type=int, required=True,
                   help="Random seed (makes the truth set fully reproducible).")
    p.add_argument("--out", required=True,
                   help="Output VCF path (plain text; workflow handles bgzip/tabix).")
    p.add_argument("--sample-name", default="simulated")
    return p.parse_args()


def load_single_contig(ref_path):
    """Return (contig_name, length, uppercase sequence string, byte array)."""
    fa = pysam.FastaFile(ref_path)
    if len(fa.references) != 1:
        sys.exit("ERROR: --reference must contain exactly one sequence; found %d (%s)"
                 % (len(fa.references), ",".join(fa.references)))
    name = fa.references[0]
    seq = fa.fetch(name).upper()
    arr = np.frombuffer(seq.encode("ascii"), dtype="S1")
    return name, len(seq), seq, arr


def choose_snps(eligible, n, rng, occupied):
    """Pick n single-base positions from `eligible` (non-N) not already occupied,
    mark them occupied, and return the array of 0-based positions."""
    # eligible are non-N; at the point SNPs are chosen only N is occupied, so a plain
    # sample suffices. We still guard against occupied for safety.
    free = eligible[~occupied[eligible]]
    if n > free.size:
        sys.exit("ERROR: requested %d SNP sites but only %d free bases." % (n, free.size))
    pos = rng.choice(free, size=n, replace=False)
    occupied[pos] = True
    return pos


def place_indels(n, rng, seq, length, eligible, occupied, min_sz, max_sz):
    """Place n non-overlapping small indels by rejection sampling. Returns a list of
    (pos0, ref, alt) tuples and marks occupied spans. Anchor positions are drawn from
    `eligible` (non-N); candidates overlapping occupied bases or N are rejected."""
    out = []
    attempts = 0
    max_attempts = max(10000, n * 500)
    while len(out) < n and attempts < max_attempts:
        attempts += 1
        p0 = int(eligible[rng.integers(0, eligible.size)])
        size = int(rng.integers(min_sz, max_sz + 1))
        if rng.integers(0, 2) == 0:
            # insertion: anchor base + inserted sequence; occupies only the anchor
            if occupied[p0] or seq[p0] not in BASES:
                continue
            ins = "".join(BASES[i] for i in rng.integers(0, 4, size=size))
            ref, alt = seq[p0], seq[p0] + ins
            occupied[p0] = True
        else:
            # deletion: anchor + `size` deleted bases; occupies the whole REF span
            if p0 + size + 1 > length:
                continue
            span = seq[p0:p0 + size + 1]
            if occupied[p0:p0 + size + 1].any() or any(c not in BASES for c in span):
                continue
            ref, alt = span, seq[p0]
            occupied[p0:p0 + size + 1] = True
        out.append((p0, ref, alt))
    if len(out) < n:
        sys.exit("ERROR: could only place %d of %d indels after %d attempts "
                 "(sequence too crowded?)." % (len(out), n, attempts))
    return out


def main():
    args = parse_args()
    rng = np.random.default_rng(args.seed)

    contig, length, seq, arr = load_single_contig(args.reference)

    # eligible = non-N positions (0-based); occupied tracks bases already used
    valid = (arr == b"A") | (arr == b"C") | (arr == b"G") | (arr == b"T")
    eligible = np.where(valid)[0]
    occupied = ~valid  # N bases start occupied so nothing is placed on them

    n_snp_div = int(round(args.divergence * length))
    n_snp_het = int(round(args.heterozygosity * length))
    n_ind_div = int(round(args.indel_divergence * length))
    n_ind_het = int(round(args.indel_heterozygosity * length))

    # --- choose positions -------------------------------------------------
    # SNPs first (bulk), then indels by rejection so they avoid SNPs and each other.
    snp_pos = choose_snps(eligible, n_snp_div + n_snp_het, rng, occupied)
    snp_div_pos = set(int(x) for x in snp_pos[:n_snp_div])          # 1|1
    # remainder are het (0|1); membership test tells the two apart when writing

    indel_div = place_indels(n_ind_div, rng, seq, length, eligible, occupied,
                             args.min_indel_size, args.max_indel_size)
    indel_het = place_indels(n_ind_het, rng, seq, length, eligible, occupied,
                             args.min_indel_size, args.max_indel_size)

    # --- assemble records: (pos0, id, ref, alt, gt) -----------------------
    records = []
    n_div_c = n_het_c = 0
    for p0 in snp_pos:
        p0 = int(p0)
        ref = seq[p0]
        alt = ALT_CHOICES[ref][int(rng.integers(0, 3))]
        if p0 in snp_div_pos:
            records.append((p0, "div_snp_%d" % n_div_c, ref, alt, "1|1")); n_div_c += 1
        else:
            gt = "0|1" if rng.integers(0, 2) == 0 else "1|0"
            records.append((p0, "het_snp_%d" % n_het_c, ref, alt, gt)); n_het_c += 1

    for i, (p0, ref, alt) in enumerate(indel_div):
        records.append((p0, "div_indel_%d" % i, ref, alt, "1|1"))
    for i, (p0, ref, alt) in enumerate(indel_het):
        gt = "0|1" if rng.integers(0, 2) == 0 else "1|0"
        records.append((p0, "het_indel_%d" % i, ref, alt, gt))

    records.sort(key=lambda r: r[0])  # by position

    # --- write VCF --------------------------------------------------------
    with open(args.out, "w") as out:
        out.write("##fileformat=VCFv4.1\n")
        out.write("##source=make_truth_vcf.py (v2 mason pipeline)\n")
        out.write("##reference=%s\n" % args.reference)
        out.write("##contig=<ID=%s,length=%d>\n" % (contig, length))
        out.write("##make_truth_vcf_param=snp_divergence(D)=%g\n" % args.divergence)
        out.write("##make_truth_vcf_param=snp_heterozygosity(H)=%g\n" % args.heterozygosity)
        out.write("##make_truth_vcf_param=indel_divergence(Di)=%g\n" % args.indel_divergence)
        out.write("##make_truth_vcf_param=indel_heterozygosity(Hi)=%g\n" % args.indel_heterozygosity)
        out.write("##make_truth_vcf_param=indel_size_range=%d-%d\n"
                  % (args.min_indel_size, args.max_indel_size))
        out.write("##make_truth_vcf_param=seed=%d\n" % args.seed)
        out.write("##make_truth_vcf_realized=contig_length=%d\n" % length)
        out.write("##make_truth_vcf_realized=n_snp_div(1|1)=%d\n" % n_snp_div)
        out.write("##make_truth_vcf_realized=n_snp_het(0|1)=%d\n" % n_snp_het)
        out.write("##make_truth_vcf_realized=n_indel_div(1|1)=%d\n" % n_ind_div)
        out.write("##make_truth_vcf_realized=n_indel_het(0|1)=%d\n" % n_ind_het)
        out.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n"
                  % args.sample_name)
        for p0, vid, ref, alt, gt in records:
            out.write("%s\t%d\t%s\t%s\t%s\t.\tPASS\t.\t.\t%s\n"
                      % (contig, p0 + 1, vid, ref, alt, gt))

    # --- data dictionary --------------------------------------------------
    with open(args.out + ".dict.txt", "w") as d:
        d.write("Data dictionary for %s\n" % args.out)
        d.write("Controlled truth (golden) VCF for the v2 mason reference-bias pipeline.\n\n")
        d.write("Generation parameters:\n")
        d.write("  reference          : %s (contig %s, length %d)\n" % (args.reference, contig, length))
        d.write("  SNP divergence D   : %g  -> %d homozygous-alt (1|1) SNP sites\n" % (args.divergence, n_snp_div))
        d.write("  SNP heterozygosity : %g  -> %d heterozygous (0|1) SNP sites\n" % (args.heterozygosity, n_snp_het))
        d.write("  indel divergence Di: %g  -> %d homozygous-alt (1|1) indel sites\n" % (args.indel_divergence, n_ind_div))
        d.write("  indel het Hi       : %g  -> %d heterozygous (0|1) indel sites\n" % (args.indel_heterozygosity, n_ind_het))
        d.write("  indel size range   : %d-%d bp (insertions and deletions, 50/50)\n" % (args.min_indel_size, args.max_indel_size))
        d.write("  seed               : %d\n\n" % args.seed)
        d.write("Fields:\n")
        d.write("  CHROM   : contig name (matches the mapping reference)\n")
        d.write("  POS     : 1-based reference position (anchor base for indels)\n")
        d.write("  ID      : {div,het}_{snp,indel}_* -- variant class and index\n")
        d.write("  REF/ALT : SNP = single base swap; deletion = REF spans deleted bases, ALT = anchor;\n")
        d.write("            insertion = REF = anchor, ALT = anchor + inserted bases\n")
        d.write("  GT      : 1|1 = fixed divergence (both haplotypes); 0|1/1|0 = heterozygous\n\n")
        d.write("Definitions:\n")
        d.write("  Divergence (D, Di) = fraction of sites fixed for a non-reference allele.\n")
        d.write("  Heterozygosity (H, Hi) = fraction of sites polymorphic within the sample.\n")
        d.write("  All variant sites are disjoint (no overlapping SNPs/indels).\n")

    sys.stderr.write(
        "# contig=%s length=%d  snp_div=%d snp_het=%d indel_div=%d indel_het=%d seed=%d\n"
        % (contig, length, n_snp_div, n_snp_het, n_ind_div, n_ind_het, args.seed))


if __name__ == "__main__":
    main()
