#!/usr/bin/env python3
"""
read_site_linkage.py -- link variant-call outcomes to read-level fate at each
divergence site (v2 mason pipeline). Answers: "when a divergence site is called as
reference instead of its true alternate, what happened to the reads that should have
supported the alternate allele?"

Every truth DIVERGENCE site is homozygous-alternate, so EVERY read that truly covers
it carries the alternate allele. For a sampled set of these sites we trace each
truly-covering read (identified from the mason truth BAM, which is in reference
coordinates) to its fate in the bwa-mapped BAM:

  retained_alt   read is aligned across the site and shows the ALT base (support kept)
  retained_ref   read is aligned across the site but shows the REF base (sequencing
                 error to reference; a small allelic-loss contribution)
  retained_other read is aligned across the site but shows a third base (seq error)
  clipped        read mapped to ~the right place (start within one read length of
                 truth) but the site base is NOT in its aligned span -> soft-clipped
                 over the site (alt evidence removed)
  mismapped      read mapped, but far (>1 read length) from its true start
  unmapped       read is unmapped
  missing        read absent from the mapped BAM entirely (expected 0)

Sites are stratified by their CALL OUTCOME (recovered = div_tp vs lost =
div_ref_or_absent) and analysed at the SAME divergence, so differences in read fate
between recovered and lost sites are not confounded by divergence itself.

We also record the caller's-eye view: the composition of ALL reads in the mapped
pileup at the site (including reads mismapped IN from elsewhere), and the resulting
alternate-allele fraction (VAF) the genotyper effectively sees.

SAMPLING: to keep this tractable over millions of sites, we randomly sample up to
--n-per-class sites from EACH outcome class (recovered / lost) per divergence level,
using --seed. This is a diagnostic estimate of mean read fate, not an exhaustive census.

Output: <out> (one row per sampled site) + <out>.dict.txt data dictionary.
"""

import sys
import argparse
import random
import pysam


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--truth-vcf", required=True)
    p.add_argument("--called-vcf", required=True)
    p.add_argument("--truth-bam", required=True)
    p.add_argument("--mapped-bam", required=True)
    p.add_argument("--region", required=True)
    p.add_argument("--coverage", required=True)
    p.add_argument("--divergence", required=True)
    p.add_argument("--heterozygosity", required=True)
    p.add_argument("--n-per-class", type=int, default=10000,
                   help="max sites sampled per outcome class per divergence")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--out", required=True)
    return p.parse_args()


TOL = 150       # bp; a read mapped within one read length of truth is "correctly placed"
MAPQ_MIN = 20   # GATK HaplotypeCaller default --minimum-mapping-quality; reads below
                # this are ignored by the caller even if they appear in a naive pileup


def read_key(aln):
    """(fragment number, mate) -> single int; mason names reads 'simulated.<N>'."""
    n = int(aln.query_name.rsplit(".", 1)[1])
    return n * 2 + (0 if aln.is_read1 else 1)


def load_truth_div_sites(path):
    """pos -> (ref_base, alt_base) for homozygous-alt (1|1) truth sites. Also the contig."""
    sites = {}
    contig = None
    opener = pysam.BGZFile(path, "rb") if path.endswith(".gz") else open(path, "rb")
    with opener as fh:
        for raw in fh:
            line = raw.decode("ascii") if isinstance(raw, bytes) else raw
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[9] == "1|1":                     # divergence (homozygous-alt) site
                sites[int(f[1])] = (f[3], f[4])
                contig = f[0]
    return sites, contig


def called_category(rec):
    """Return (category, alt_base) for a called record: hom_ref/het/hom_alt/other."""
    gt = rec.samples[0].get("GT", (None, None))
    if gt is None or all(a is None for a in gt):
        return "hom_ref", None
    idx = [a for a in gt if a is not None]
    if all(a == 0 for a in idx):
        return "hom_ref", None
    if len(set(idx)) == 1 and idx[0] != 0:
        allele = rec.alleles[idx[0]]
        return ("hom_alt", allele) if len(allele) == 1 else ("other", None)
    nonzero = [a for a in idx if a != 0]
    if 0 in idx and len(set(nonzero)) == 1:
        allele = rec.alleles[nonzero[0]]
        return ("het", allele) if len(allele) == 1 else ("other", None)
    return "other", None


def classify_outcome(truth_alt, call):
    """Outcome of a divergence (hom-alt) site given its call. 'lost' == ref or absent."""
    if call is None:
        return "lost"
    cat, alt = call
    if cat == "hom_alt":
        return "tp" if alt == truth_alt else "other"
    if cat == "hom_ref":
        return "lost"
    if cat == "het":
        return "het"
    return "other"


def pileup_bases(bam, contig, pos1):
    """read_key -> (aligned base upper, mapping quality) for all primary reads aligned
    across a 1-based position. Reads with the position soft-clipped or deleted do not
    appear (that is exactly the soft-clip signal we want to detect elsewhere)."""
    out = {}
    for col in bam.pileup(contig, pos1 - 1, pos1, truncate=True, min_base_quality=0):
        for pr in col.pileups:
            a = pr.alignment
            if a.is_secondary or a.is_supplementary:
                continue
            if pr.is_del or pr.is_refskip or pr.query_position is None:
                continue
            out[read_key(a)] = (a.query_sequence[pr.query_position].upper(),
                                a.mapping_quality)
    return out


def pileup_truth_cover(bam, contig, pos1):
    """read_key -> reference_start for all primary truth reads covering a position."""
    out = {}
    for col in bam.pileup(contig, pos1 - 1, pos1, truncate=True, min_base_quality=0):
        for pr in col.pileups:
            a = pr.alignment
            if a.is_secondary or a.is_supplementary:
                continue
            out[read_key(a)] = a.reference_start
    return out


def main():
    args = parse_args()
    rng = random.Random(args.seed)

    # --- truth divergence sites + their call outcome ----------------------
    truth_sites, contig = load_truth_div_sites(args.truth_vcf)
    if not truth_sites:
        sys.exit("ERROR: no divergence (1|1) sites in truth VCF (is divergence 0?).")

    called = {}
    vf = pysam.VariantFile(args.called_vcf)
    for rec in vf:
        if rec.ref is None or len(rec.ref) != 1:
            continue
        called[rec.pos] = called_category(rec)
    vf.close()

    by_outcome = {"tp": [], "lost": []}
    for pos, (ref_b, alt_b) in truth_sites.items():
        oc = classify_outcome(alt_b, called.get(pos))
        if oc in by_outcome:
            by_outcome[oc].append(pos)

    # --- stratified sample of sites per outcome ---------------------------
    sampled = []
    for oc, positions in by_outcome.items():
        rng.shuffle(positions)
        chosen = positions[:args.n_per_class]
        sampled.extend((pos, oc) for pos in chosen)
        sys.stderr.write("# outcome=%s available=%d sampled=%d\n"
                         % (oc, len(positions), len(chosen)))
    sampled.sort()  # by position, for cache-friendly BAM access

    # --- global mapped-read start positions (to explain non-covering reads) --
    # read_key -> reference_start (or -1 if unmapped). Reads absent from this dict
    # are 'missing' from the mapped BAM entirely.
    mapped_start = {}
    with pysam.AlignmentFile(args.mapped_bam, "rb") as mb:
        for a in mb:
            if a.is_secondary or a.is_supplementary:
                continue
            mapped_start[read_key(a)] = -1 if a.is_unmapped else a.reference_start
    sys.stderr.write("# mapped reads indexed: %d\n" % len(mapped_start))

    # --- per-site read-fate tracing --------------------------------------
    fate_cols = ["retained_alt", "retained_ref", "retained_other",
                 "clipped", "mismapped", "unmapped", "missing"]
    header = (["region", "coverage", "divergence", "heterozygosity",
               "pos", "outcome", "truth_depth"] + fate_cols +
              ["retention", "true_alt_frac",
               "pileup_depth", "pileup_alt", "pileup_ref", "pileup_alt_frac",
               "pileup_depth_hq", "pileup_alt_hq", "pileup_alt_frac_hq"])

    tbam = pysam.AlignmentFile(args.truth_bam, "rb")
    mbam = pysam.AlignmentFile(args.mapped_bam, "rb")

    with open(args.out, "w") as out:
        out.write(",".join(header) + "\n")
        for pos, outcome in sampled:
            ref_b, alt_b = truth_sites[pos]
            truth_cov = pileup_truth_cover(tbam, contig, pos)   # true covering reads
            mapped_base = pileup_bases(mbam, contig, pos)       # reads aligned across pos

            fate = {k: 0 for k in fate_cols}
            for k, tstart in truth_cov.items():
                if k in mapped_base:
                    b = mapped_base[k][0]
                    if b == alt_b:
                        fate["retained_alt"] += 1
                    elif b == ref_b:
                        fate["retained_ref"] += 1
                    else:
                        fate["retained_other"] += 1
                else:
                    ms = mapped_start.get(k, None)
                    if ms is None:
                        fate["missing"] += 1
                    elif ms < 0:
                        fate["unmapped"] += 1
                    elif abs(ms - tstart) <= TOL:
                        fate["clipped"] += 1     # placed right, site base soft-clipped
                    else:
                        fate["mismapped"] += 1

            truth_depth = len(truth_cov)
            retained = fate["retained_alt"] + fate["retained_ref"] + fate["retained_other"]
            retention = retained / truth_depth if truth_depth else 0.0
            true_alt_frac = fate["retained_alt"] / retained if retained else 0.0

            # caller's-eye pileup (all reads aligned at pos, incl. mismapped-in).
            # The "hq" variants keep only MAPQ >= MAPQ_MIN reads -- what the genotyper
            # actually uses -- so sites lost to low mapping quality are distinguishable
            # from sites lost to soft-clip/mismap dropout.
            pileup_depth = len(mapped_base)
            pileup_alt = sum(1 for b, q in mapped_base.values() if b == alt_b)
            pileup_ref = sum(1 for b, q in mapped_base.values() if b == ref_b)
            pileup_alt_frac = pileup_alt / pileup_depth if pileup_depth else 0.0

            hq = [(b, q) for b, q in mapped_base.values() if q >= MAPQ_MIN]
            pileup_depth_hq = len(hq)
            pileup_alt_hq = sum(1 for b, q in hq if b == alt_b)
            pileup_alt_frac_hq = pileup_alt_hq / pileup_depth_hq if pileup_depth_hq else 0.0

            row = ([args.region, args.coverage, args.divergence, args.heterozygosity,
                    pos, outcome, truth_depth] +
                   [fate[k] for k in fate_cols] +
                   ["%.4f" % retention, "%.4f" % true_alt_frac,
                    pileup_depth, pileup_alt, pileup_ref, "%.4f" % pileup_alt_frac,
                    pileup_depth_hq, pileup_alt_hq, "%.4f" % pileup_alt_frac_hq])
            out.write(",".join(str(x) for x in row) + "\n")

    tbam.close()
    mbam.close()

    # --- data dictionary --------------------------------------------------
    with open(args.out + ".dict.txt", "w") as d:
        d.write("Data dictionary for %s\n" % args.out)
        d.write("One row per SAMPLED truth divergence (1|1) site. Every such site is\n")
        d.write("homozygous-alt, so every truly-covering read should carry the ALT base.\n")
        d.write("Sites are stratified by call outcome and sampled up to %d per class.\n\n" % args.n_per_class)
        d.write("region,coverage,divergence,heterozygosity : simulation condition\n")
        d.write("pos            : 1-based reference position of the divergence site\n")
        d.write("outcome        : tp = called hom-alt (recovered); lost = called hom-ref or not called\n")
        d.write("truth_depth    : number of reads that TRULY cover the site (all carry ALT)\n")
        d.write("retained_alt   : truly-covering reads aligned across the site showing ALT (support kept)\n")
        d.write("retained_ref   : truly-covering reads aligned across the site showing REF (seq error)\n")
        d.write("retained_other : truly-covering reads aligned across the site showing a third base\n")
        d.write("clipped        : truly-covering reads placed ~correctly but with the site base soft-clipped\n")
        d.write("mismapped      : truly-covering reads mapped > %d bp from their true start\n" % TOL)
        d.write("unmapped       : truly-covering reads that are unmapped\n")
        d.write("missing        : truly-covering reads absent from the mapped BAM (expected 0)\n")
        d.write("retention      : (retained_alt+ref+other) / truth_depth (fraction whose base survived at the site)\n")
        d.write("true_alt_frac  : retained_alt / retained (ALT fraction among surviving true reads)\n")
        d.write("pileup_depth   : ALL reads aligned across the site in the mapped BAM (incl. mismapped-in)\n")
        d.write("pileup_alt     : those pileup reads showing ALT\n")
        d.write("pileup_ref     : those pileup reads showing REF\n")
        d.write("pileup_alt_frac: pileup_alt / pileup_depth (naive ALT fraction, any MAPQ)\n")
        d.write("pileup_depth_hq: pileup reads with MAPQ >= %d (what the caller actually uses)\n" % MAPQ_MIN)
        d.write("pileup_alt_hq  : high-MAPQ pileup reads showing ALT\n")
        d.write("pileup_alt_frac_hq: pileup_alt_hq / pileup_depth_hq (ALT fraction the genotyper effectively sees;\n")
        d.write("                 a site with high pileup_alt_frac but low pileup_alt_frac_hq was lost to low MAPQ)\n")

    print("wrote %s (%d sampled sites)" % (args.out, len(sampled)))


if __name__ == "__main__":
    main()
