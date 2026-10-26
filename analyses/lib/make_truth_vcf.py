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

sys.path.insert(0, __file__.rsplit("/", 1)[0])  # allow import of sibling module
import divergence_field as dfield
import provenance                       # shared provenance block for the data dictionaries

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
    # --- structural variants (large del/ins) as a CONFOUNDER (s07) --------------
    # SVs are rare and large, so they are specified as a COUNT (not a rate) with their
    # own size range. Homozygous (1|1) fixed differences, placed disjoint from SNPs and
    # small indels. Our small-variant callers cannot call them; they enter the truth only
    # to disrupt read mapping in their neighborhood (deletion = coverage gap, insertion =
    # unmappable novel sequence, both = clipped breakpoint reads).
    p.add_argument("--sv-count", type=int, default=0,
                   help="number of large SVs to place (~half deletions, half insertions).")
    p.add_argument("--sv-min-size", type=int, default=50)
    p.add_argument("--sv-max-size", type=int, default=5000,
                   help="SV size range in bp (uniform; insertions get random novel sequence).")
    p.add_argument("--seed", type=int, required=True,
                   help="Random seed (makes the truth set fully reproducible).")
    p.add_argument("--out", required=True,
                   help="Output VCF path (plain text; workflow handles bgzip/tabix).")
    p.add_argument("--sample-name", default="simulated")
    # --- s05 heterogeneous divergence (continuous autocorrelated field) -------
    # When --divergence-sigma > 0, the DIVERGENT (1|1) SNP positions are drawn with
    # probability proportional to an autocorrelated field instead of uniformly, so
    # divergence clusters/thins along the contig (see divergence_field.py). The mean
    # divergence is unchanged (still exactly round(D*L) sites). Het SNPs and indels
    # remain UNIFORM -- only fixed divergence is made heterogeneous. sigma == 0 (the
    # default) reproduces the uniform placement of s01-s04 byte-for-byte.
    p.add_argument("--divergence-sigma", type=float, default=0.0,
                   help="Field amplitude (0 = uniform/homogeneous; >0 = heterogeneous div).")
    p.add_argument("--divergence-corr-len", type=float, default=20000.0,
                   help="Field correlation length in bp (scale of heterogeneity). "
                        "Ignored when --divergence-sigma == 0.")
    p.add_argument("--divergence-field-seed", type=int, default=None,
                   help="Seed for the latent field. Defaults to --seed if unset. "
                        "Ignored when --divergence-sigma == 0.")
    p.add_argument("--divergence-cap", type=float, default=None,
                   help="Biological ceiling on LOCAL divergence (e.g. 0.20): bounds the "
                        "field's hot regions instead of letting the log-normal tail run "
                        "to saturation. Realized ceiling is validated via "
                        "make_divergence_field. Ignored when --divergence-sigma == 0.")
    p.add_argument("--het-on-field", action="store_true",
                   help="Place heterozygous SNPs on the SAME field as divergence "
                        "(diversity co-varies with divergence -- realistic). Default: "
                        "het uniform. Ignored when --divergence-sigma == 0.")
    p.add_argument("--indel-on-field", action="store_true",
                   help="Place indels on the SAME field as divergence (indel divergence "
                        "co-clusters with SNP divergence -- the s05 het-indels scenario). "
                        "Default: indels uniform. Ignored when --divergence-sigma == 0.")
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


def choose_snps_field(eligible, n_div, n_het, rng, occupied, weights, het_on_field):
    """s05 variant of choose_snps: place the n_div FIXED-divergence SNPs with
    probability proportional to `weights` (the autocorrelated field). Het SNPs are
    placed the same way (on the SAME field, so diversity co-varies with divergence
    -- the realistic default) when het_on_field is True, or uniformly otherwise.

    Returns (div_pos, het_pos) as 0-based position arrays and marks them occupied.
    Divergence and het draws are independent (different random keys) so the two
    site sets are disjoint but share the field's hot/cold structure.
    """
    free = eligible[~occupied[eligible]]
    if n_div + n_het > free.size:
        sys.exit("ERROR: requested %d SNP sites but only %d free bases."
                 % (n_div + n_het, free.size))
    # divergent sites: weighted, without replacement -> density follows the field
    div_local = dfield.weighted_sample_without_replacement(weights[free], n_div, rng)
    div_pos = free[div_local]
    occupied[div_pos] = True
    # het sites: on the same field (correlated) or uniform, from what remains
    if het_on_field:
        free2 = eligible[~occupied[eligible]]
        het_local = dfield.weighted_sample_without_replacement(weights[free2], n_het, rng)
        het_pos = free2[het_local]
        occupied[het_pos] = True
    else:
        het_pos = choose_snps(eligible, n_het, rng, occupied)
    return div_pos, het_pos


def _try_indel(p0, rng, seq, length, occupied, min_sz, max_sz):
    """Attempt one indel at anchor p0. Returns (ref, alt) and marks occupied on
    success, or None if it can't be placed (overlap / N / off the end). Draws the
    size and ins/del choice exactly as the original loop did, so the uniform path's
    RNG stream is unchanged."""
    size = int(rng.integers(min_sz, max_sz + 1))
    if rng.integers(0, 2) == 0:
        # insertion: anchor base + inserted sequence; occupies only the anchor
        if occupied[p0] or seq[p0] not in BASES:
            return None
        ins = "".join(BASES[i] for i in rng.integers(0, 4, size=size))
        occupied[p0] = True
        return seq[p0], seq[p0] + ins
    # deletion: anchor + `size` deleted bases; occupies the whole REF span
    if p0 + size + 1 > length:
        return None
    span = seq[p0:p0 + size + 1]
    if occupied[p0:p0 + size + 1].any() or any(c not in BASES for c in span):
        return None
    occupied[p0:p0 + size + 1] = True
    return span, seq[p0]


def place_indels(n, rng, seq, length, eligible, occupied, min_sz, max_sz, weights=None):
    """Place n non-overlapping small indels. Returns a list of (pos0, ref, alt)
    tuples and marks occupied spans.

    weights is None (default): uniform rejection sampling -- UNCHANGED from the
    original (preserves s01-s04 reproducibility). weights given (s05 --indel-on-
    field): anchors are visited in field-WEIGHTED random order (Efraimidis-Spirakis
    keys), so indel divergence co-clusters with SNP divergence."""
    out = []
    if weights is None:
        attempts = 0
        max_attempts = max(10000, n * 500)
        while len(out) < n and attempts < max_attempts:
            attempts += 1
            p0 = int(eligible[rng.integers(0, eligible.size)])
            ra = _try_indel(p0, rng, seq, length, occupied, min_sz, max_sz)
            if ra is not None:
                out.append((p0, ra[0], ra[1]))
        if len(out) < n:
            sys.exit("ERROR: could only place %d of %d indels after %d attempts "
                     "(sequence too crowded?)." % (len(out), n, attempts))
    else:
        # weighted-random visiting order over eligible anchors (rejection = skip)
        keys = rng.standard_exponential(eligible.size) / weights[eligible]
        order = eligible[np.argsort(keys, kind="stable")]
        for p0 in order:
            if len(out) >= n:
                break
            ra = _try_indel(int(p0), rng, seq, length, occupied, min_sz, max_sz)
            if ra is not None:
                out.append((int(p0), ra[0], ra[1]))
        if len(out) < n:
            sys.exit("ERROR: could only place %d of %d field-weighted indels "
                     "(sequence too crowded?)." % (len(out), n))
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
    # SVs FIRST (they are large, so reserve their spans before the dense SNP/indel
    # placement -- otherwise a big deletion's span always overlaps a SNP and is rejected).
    # Then SNPs (bulk), then small indels; each avoids everything already occupied.
    # place_indels(0, ...) consumes no RNG, so sv_count=0 leaves the SNP/indel stream
    # (and s01-s05 reproducibility) unchanged.
    sv = (place_indels(args.sv_count, rng, seq, length, eligible, occupied,
                       args.sv_min_size, args.sv_max_size)
          if args.sv_count > 0 else [])

    field_active = args.divergence_sigma > 0.0
    weights = None
    if field_active:
        # s05: divergent (1|1) sites follow the autocorrelated field; het on the
        # same field (correlated) or uniform per --het-on-field.
        field_seed = args.divergence_field_seed if args.divergence_field_seed is not None else args.seed
        cap_ratio = (args.divergence_cap / args.divergence) if args.divergence_cap else None
        weights = dfield.build_weights(length, args.divergence_corr_len,
                                       args.divergence_sigma, field_seed,
                                       cap_ratio=cap_ratio)
        div_pos, het_pos = choose_snps_field(eligible, n_snp_div, n_snp_het,
                                             rng, occupied, weights, args.het_on_field)
        snp_pos = np.concatenate([div_pos, het_pos])
        snp_div_pos = set(int(x) for x in div_pos)                 # 1|1
    else:
        # s01-s04: uniform placement (unchanged). First n_snp_div are 1|1.
        snp_pos = choose_snps(eligible, n_snp_div + n_snp_het, rng, occupied)
        snp_div_pos = set(int(x) for x in snp_pos[:n_snp_div])     # 1|1
        # remainder are het (0|1); membership test tells the two apart when writing

    # indels co-cluster with divergence when --indel-on-field (else uniform)
    ind_weights = weights if (field_active and args.indel_on_field) else None
    indel_div = place_indels(n_ind_div, rng, seq, length, eligible, occupied,
                             args.min_indel_size, args.max_indel_size, weights=ind_weights)
    indel_het = place_indels(n_ind_het, rng, seq, length, eligible, occupied,
                             args.min_indel_size, args.max_indel_size, weights=ind_weights)

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

    for i, (p0, ref, alt) in enumerate(sv):
        kind = "del" if len(ref) > len(alt) else "ins"   # del: REF spans; ins: ALT spans
        records.append((p0, "div_sv_%s_%d" % (kind, i), ref, alt, "1|1"))

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
        out.write("##make_truth_vcf_param=divergence_sigma=%g\n" % args.divergence_sigma)
        if field_active:
            out.write("##make_truth_vcf_param=divergence_corr_len_bp=%g\n" % args.divergence_corr_len)
            out.write("##make_truth_vcf_param=divergence_field_seed=%d\n"
                      % (args.divergence_field_seed if args.divergence_field_seed is not None else args.seed))
            out.write("##make_truth_vcf_param=divergence_cap=%s\n"
                      % ("%g" % args.divergence_cap if args.divergence_cap else "none"))
            out.write("##make_truth_vcf_param=het_on_field=%s\n" % bool(args.het_on_field))
            out.write("##make_truth_vcf_param=indel_on_field=%s\n" % bool(args.indel_on_field))
            out.write("##make_truth_vcf_note=fixed-divergence(1|1) SNP sites placed by "
                      "autocorrelated field; het sites %s; indels %s\n"
                      % ("on same field" if args.het_on_field else "uniform",
                         "on same field" if args.indel_on_field else "uniform"))
        out.write("##make_truth_vcf_realized=contig_length=%d\n" % length)
        out.write("##make_truth_vcf_realized=n_snp_div(1|1)=%d\n" % n_snp_div)
        out.write("##make_truth_vcf_realized=n_snp_het(0|1)=%d\n" % n_snp_het)
        out.write("##make_truth_vcf_realized=n_indel_div(1|1)=%d\n" % n_ind_div)
        out.write("##make_truth_vcf_realized=n_indel_het(0|1)=%d\n" % n_ind_het)
        if args.sv_count > 0:
            out.write("##make_truth_vcf_param=sv_count=%d\n" % args.sv_count)
            out.write("##make_truth_vcf_param=sv_size_range=%d-%d\n"
                      % (args.sv_min_size, args.sv_max_size))
            out.write("##make_truth_vcf_realized=n_sv(1|1)=%d\n" % len(sv))
            out.write("##make_truth_vcf_note=SVs (div_sv_{del,ins}_*) are a CONFOUNDER: "
                      "homozygous large del/ins, uncallable by small-variant callers\n")
        out.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n"
                  % args.sample_name)
        for p0, vid, ref, alt, gt in records:
            out.write("%s\t%d\t%s\t%s\t%s\t.\tPASS\t.\t.\t%s\n"
                      % (contig, p0 + 1, vid, ref, alt, gt))

    # --- data dictionary --------------------------------------------------
    with open(args.out + ".dict.txt", "w") as d:
        d.write(provenance.header(args.out))
        d.write("Controlled truth (golden) VCF for the v2 mason reference-bias pipeline.\n\n")
        d.write("Generation parameters:\n")
        d.write("  reference          : %s (contig %s, length %d)\n" % (args.reference, contig, length))
        d.write("  SNP divergence D   : %g  -> %d homozygous-alt (1|1) SNP sites\n" % (args.divergence, n_snp_div))
        d.write("  SNP heterozygosity : %g  -> %d heterozygous (0|1) SNP sites\n" % (args.heterozygosity, n_snp_het))
        d.write("  indel divergence Di: %g  -> %d homozygous-alt (1|1) indel sites\n" % (args.indel_divergence, n_ind_div))
        d.write("  indel het Hi       : %g  -> %d heterozygous (0|1) indel sites\n" % (args.indel_heterozygosity, n_ind_het))
        d.write("  indel size range   : %d-%d bp (insertions and deletions, 50/50)\n" % (args.min_indel_size, args.max_indel_size))
        d.write("  seed               : %d\n" % args.seed)
        if field_active:
            d.write("  divergence field   : autocorrelated (s05); sigma=%g corr_len=%g bp field_seed=%d\n"
                    % (args.divergence_sigma, args.divergence_corr_len,
                       args.divergence_field_seed if args.divergence_field_seed is not None else args.seed))
            d.write("                       cap=%s (local-D ceiling); het_on_field=%s; indel_on_field=%s\n"
                    % ("%g" % args.divergence_cap if args.divergence_cap else "none",
                       bool(args.het_on_field), bool(args.indel_on_field)))
            d.write("                       -> fixed-divergence (1|1) SNP sites clustered per the field\n")
            d.write("                          (%s); indels uniform; mean D unchanged.\n"
                    % ("het on the same field" if args.het_on_field else "het uniform"))
        else:
            d.write("  divergence field   : none (uniform placement; sigma=0)\n")
        d.write("\n")
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
        "# contig=%s length=%d  snp_div=%d snp_het=%d indel_div=%d indel_het=%d sv=%d seed=%d\n"
        % (contig, length, n_snp_div, n_snp_het, n_ind_div, n_ind_het, len(sv), args.seed))


if __name__ == "__main__":
    main()
