#!/usr/bin/env python3
"""
variant_read_fate.py -- the EXACT per-SNP read ledger for iterative mapping.

A SNP is not one class: it is a *distribution over its covering reads*. In a simulation every
fact is deterministic, so we do not model or collapse -- we count, per (snp, iteration), how
many of the reads that TRULY cover it (by mason-truth position) fall in each cell of the full
cross, and attach the caller's exact call. Fractions are derived from the counts downstream;
nothing is reduced to a dominant label here.

Two independent, deterministic axes per read (+ the base it presents):
  placement (alignment start vs truth): exact (delta 0) | shifted (0<delta<=150) |
                                        mismapped (>150) | diff_chr | unmapped
  pos_at_p  (the SNP's own coordinate p in this read): aligned | clipped | deleted | absent
  base_at_p (only when aligned): alt (matches truth) | ref | err (neither)
read_clip (does the read carry any soft-clip) is kept for mismapped reads, since "mismapped
AND clipped" is a distinct, real cell (relevant to where false positives seed).

Cells (columns), all read COUNTS per (snp, iteration):
  ex_alt ex_ref ex_err ex_clip ex_del   -- exact placement, p aligned(alt/ref/err) / clipped / deleted
  sh_alt sh_ref sh_err sh_clip sh_del   -- shifted placement, same breakdown
  mm_clip mm_noclip                     -- mismapped, split by whether the read itself is soft-clipped
  dc un                                 -- diff_chr, unmapped
plus recovered (0/1) and the truth ref/alt. Identity coordinates only (s01 SNP-only; iter1).

Each read's TRUE span comes from its haplotype position (mason -oa tag oP), mapped back to the
reference through the haplotypes rebuilt from the mason-format truth VCF (truth_haplotypes.py),
not from the truth BAM's POS/CIGAR (1-3 bp off for ~1% of reads at dense SNP clusters). Outputs
made before 2026-10-27 used POS/CIGAR.

Usage:
  variant_read_fate.py --truth-bam T.bam --mason-vcf truth.vcf.gz --reference ref.fa \
      --truth-vcf truth.norm.vcf.gz --region 19 --divergence 0.10 --caller gatk --sample 8000 --out rows.csv \
      --iter 1:iter1/mapped.bam:iter1/recovered.vcf.gz  --iter 2:...  --iter 3:...
(--sample <=0 means ALL true SNPs, no sampling.)
"""
import argparse, bisect, sys
import pysam
import provenance                       # shared provenance block for the data dictionaries
import truth_haplotypes                 # each truth read's true span from its haplotype position (oP)

TOLERANCE = 150
CELLS = ["ex_alt", "ex_ref", "ex_err", "ex_clip", "ex_del",
         "sh_alt", "sh_ref", "sh_err", "sh_clip", "sh_del",
         "mm_clip", "mm_noclip", "dc", "un"]
CIX = {c: i for i, c in enumerate(CELLS)}


def read_key(aln):
    n = int(aln.query_name.rsplit(".", 1)[1])
    return n * 2 + (0 if aln.is_read1 else 1)


def load_hom_div_snps(vcf_path):
    """Sorted [(pos0, ref, alt, alt_haps)] for homozygous (1|1) biallelic SNP truth sites.
    alt_haps = the set of 1-based haplotypes (matching mason's oH tag) that carry the alt
    allele, from the PHASED genotype: 1|1 -> {1,2}, 1|0 -> {1}, 0|1 -> {2}. For hom sites it
    is {1,2}, so every covering read truly carries alt; kept general so het sites work too."""
    out = []
    op = pysam.BGZFile(vcf_path) if vcf_path.endswith(".gz") else open(vcf_path, "rb")
    for raw in op:
        line = raw.decode() if isinstance(raw, bytes) else raw
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        ref, alt, gt = c[3], c[4], c[9].split(":")[0]
        if len(ref) == 1 and len(alt) == 1 and gt in ("1|1", "1/1"):
            out.append((int(c[1]) - 1, ref.upper(), alt.upper(), frozenset({1, 2})))
    out.sort()
    return out


def pileup_at(mf, chrom, mapped_p, ref_p, refbase, altbase, truth):
    """POSITION-CENTRIC view of reference position p in the mapped BAM -- the caller's-eye pileup,
    i.e. EVERY read aligned over p regardless of where it truly came from. This is the complement of
    the truth-centric ledger: it counts reads MISMAPPED INTO p from elsewhere (which the ledger books
    at their own home loci), so it can show whether a SNP flips called<->not-called because the ref /
    contaminant reads at p, the depth, or the MAPQ changed while alt-support (f_alt) held steady.
    A read is MISMAPPED-IN when its true span (from the truth BAM) does not cover ref_p -- these are
    quantified directly (mmin_*), split by the base they present, so the "contamination" is measured
    not inferred. mapped_p = position in the mapped BAM (for the pileup fetch); ref_p = the same SNP in
    reference/truth coordinates (for the truth check; identical for s01 identity coords).
    Returns (depth, n_ref, n_alt, n_oth, mapq_mean, n_mapq20, mmin, mmin_ref, mmin_alt)."""
    depth = n_ref = n_alt = n_oth = n_mq20 = 0
    mmin = mmin_ref = mmin_alt = 0
    mq_sum = 0
    # samtools stepper applies mpileup-style read filtering (skips secondary/supp/dup/qcfail)
    for col in mf.pileup(chrom, mapped_p, mapped_p + 1, truncate=True, min_base_quality=0, stepper="samtools"):
        for pr in col.pileups:
            if pr.is_del or pr.is_refskip or pr.query_position is None:
                continue                                   # p deleted / skipped in this read -> no base
            aln = pr.alignment
            b = aln.query_sequence[pr.query_position].upper()
            is_alt = b == altbase
            is_ref = b == refbase
            depth += 1
            if is_alt:
                n_alt += 1
            elif is_ref:
                n_ref += 1
            else:
                n_oth += 1
            mq = aln.mapping_quality
            mq_sum += mq
            if mq >= 20:
                n_mq20 += 1
            # mismapped-IN: aligned over p but the read's TRUE span does not cover ref_p
            t = truth.get(read_key(aln))
            if t is None or t[0] != chrom or not (t[1] <= ref_p < t[2]):
                mmin += 1
                if is_ref:
                    mmin_ref += 1
                elif is_alt:
                    mmin_alt += 1
    return (depth, n_ref, n_alt, n_oth, (round(mq_sum / depth, 1) if depth else 0.0), n_mq20,
            mmin, mmin_ref, mmin_alt)


PILE_COLS = ["depth_p", "ref_p", "alt_p", "oth_p", "mapq_p", "mapq20_p",
             "mmin_p", "mmin_ref_p", "mmin_alt_p"]


def load_het_snps(vcf_path):
    """Sorted [(pos0, ref, alt, alt_haps)] for HETEROZYGOUS (0/1) biallelic SNP truth sites.
    mason applies PHASED genotypes, so the alt allele sits on exactly one haplotype: 1|0 -> {1},
    0|1 -> {2}. Thus only reads from that haplotype (mason oH tag) truly carry the alt, and
    n_true_alt = those ~half of covering reads. Unphased het (0/1, 1/0) can't be assigned a
    haplotype and is skipped (mason truth is always phased here)."""
    out = []
    op = pysam.BGZFile(vcf_path) if vcf_path.endswith(".gz") else open(vcf_path, "rb")
    for raw in op:
        line = raw.decode() if isinstance(raw, bytes) else raw
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        ref, alt, gt = c[3], c[4], c[9].split(":")[0]
        if len(ref) != 1 or len(alt) != 1:
            continue
        if gt == "1|0":
            haps = frozenset({1})
        elif gt == "0|1":
            haps = frozenset({2})
        else:
            continue                                   # only phased het sites
        out.append((int(c[1]) - 1, ref.upper(), alt.upper(), haps))
    out.sort()
    return out


def load_recovered(vcf_path, het=False):
    """Set of (pos0, ALT) the caller "recovered" at this iteration. For hom SNPs this is the
    accumulated/consensus recovered set (any GT). For HET (het=True) "recovered" means the caller
    called the site correctly AS A HET (GT 0/1) in its per-iteration VCF -- a het miscalled hom
    (0/0 missed, or 1/1 flipped to hom-alt by reference bias) does NOT count as recovered."""
    het_gts = {"0/1", "1/0", "0|1", "1|0"}
    s = set()
    op = pysam.BGZFile(vcf_path) if vcf_path.endswith(".gz") else open(vcf_path, "rb")
    for raw in op:
        line = raw.decode() if isinstance(raw, bytes) else raw
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 5:
            continue
        if len(c[3]) == 1 and len(c[4]) == 1:
            if het:
                gt = c[9].split(":")[0] if len(c) > 9 else ""
                if gt not in het_gts:
                    continue
            s.add((int(c[1]) - 1, c[4].upper()))
    return s


def make_lifter(chain_path):
    """Point-liftover from a single-contig, same-strand UCSC chain. lift(pos) maps a target
    coordinate to the query coordinate, or None if pos falls in a chain gap. We pass the
    per-iteration ref2cons.chain (target=reference, query=consensus), so lift() takes a
    REFERENCE coordinate to the mapped BAM's CONSENSUS coordinate -- needed for s02 iter>1,
    where reads are mapped to the length-shifted consensus rather than the reference."""
    starts, ends, offs = [], [], []
    t = q = 0
    with open(chain_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith("chain"):
                f = line.split()
                t, q = int(f[5]), int(f[10])          # tStart (reference), qStart (consensus)
                continue
            p = line.split()
            size = int(p[0])
            starts.append(t); ends.append(t + size); offs.append(q - t)
            if len(p) == 3:
                t += size + int(p[1]); q += size + int(p[2])

    def lift(pos):
        i = bisect.bisect_right(starts, pos) - 1
        if i >= 0 and pos < ends[i]:
            return pos + offs[i]
        return None
    return lift


def write_aggregate(a, n_iter, snp_traj, hcomp):
    """Write the two small summary tables consumed by figures G and H (aggregate mode).

    <out>-gtraj.csv : one row per distinct (call status, f_alt bin) trajectory across the iterations,
                      with n_snps = how many SNPs followed it  -> the joint-state flow (G).
    <out>-hcomp.csv : per (iteration, call status), the at-p read-fate categories summed over all
                      SNPs, plus total covering reads -> the composition (H).
    """
    from collections import Counter
    traj_counts = Counter(tuple(t) for t in snp_traj if len(t) == n_iter)
    with open(a.out + "-gtraj.csv", "w") as out:
        hdr = ["divergence", "caller"]
        for k in range(n_iter):
            hdr += ["s%d" % (k + 1), "fb%d" % (k + 1)]
        out.write(",".join(hdr + ["n_snps"]) + "\n")
        for traj, n in traj_counts.items():
            flat = [x for pair in traj for x in pair]          # (status, fb) pairs -> flat list
            out.write(",".join([a.divergence, a.caller] + flat + [str(n)]) + "\n")
    with open(a.out + "-hcomp.csv", "w") as out:
        out.write("divergence,caller,iteration,status,alt_at_p,ref_at_p,err_at_p,"
                  "clipped,deleted,mapped_elsewhere,unmapped,n_reads\n")
        for (it, status), h in sorted(hcomp.items()):
            out.write(",".join([a.divergence, a.caller, it, status] + [str(x) for x in h]) + "\n")


def write_dictionaries(a):
    """Write the data dictionary for every table this run produces (per-SNP mode: <out>.dict.txt;
    aggregate mode: <out>-gtraj.csv.dict.txt and <out>-hcomp.csv.dict.txt). Each starts with the
    shared provenance block (lib/provenance.py): this run's own command line, or with --dicts-only
    the ORIGINAL producing command, marked as a backfill (with --table-commit: as a rewrite for an
    output made by the current method at that commit)."""
    def hdr(table):
        if a.dicts_only:
            return provenance.header(table, command=a.dicts_only, backfill=True,
                                     table_commit=a.table_commit)
        return provenance.header(table)

    # how each read's true span was taken. --dicts-only alone documents outputs made before the
    # 2026-10-27 fix, which took it from the truth BAM's POS/CIGAR; with --table-commit, outputs made
    # by the current method.
    span_txt = ("mason -oa truth BAM POS/CIGAR (true start/end) + oH haplotype per read (the method\n"
                "             before 2026-10-27: POS/CIGAR is off for ~1% of reads SNP-only, ~30% in s02);"
                if a.dicts_only and not a.table_commit else
                "mason -oa truth BAM: oH haplotype + oP haplotype position per read; the true\n"
                "             reference span is [oP, oP + read length) on that haplotype, mapped back to\n"
                "             the reference (truth_haplotypes.py, haplotypes rebuilt from --mason-vcf);\n"
                "             each read is PLACED where an aligner would lay it unclipped: the start of\n"
                "             the best alignment of its leftmost 30 bases on its haplotype (bwa mem\n"
                "             scores; moves 1-3 bp off oP only for a sequencing indel among those bases,\n"
                "             ~0.02% of reads; since commit efa152b, 2026-10-27), minus any inserted\n"
                "             bases it starts on;")

    site = ("HETEROZYGOUS (phased 0|1 / 1|0) biallelic truth SNPs; true-alt reads are those from the "
            "alt-carrying haplotype (mason oH tag); recovered = called correctly AS A HET in the "
            "per-iteration caller VCF") if a.het else (
            "homozygous divergent (1|1) biallelic truth SNPs; every covering read truly carries alt; "
            "recovered = (position, alt) present in the iteration's recovered VCF")
    sampling = ("all truth SNPs on the region (--sample 0)" if not a.sample or a.sample <= 0 else
                "systematic sample of ~%d truth SNPs (every len/%d-th site along the region)"
                % (a.sample, a.sample))
    common = ("Sites      : %s.\n"
              "Sampling   : %s.\n"
              "Read truth : %s\n"
              "             reads are assigned to every SNP their TRUE span covers.\n"
              "Placement  : exact = mapped start == true start; shifted = 0 < |delta| <= %d bp;\n"
              "             mismapped = |delta| > %d bp; diff_chr = other contig (or true start in a\n"
              "             chain gap); unmapped. For s02 iter > 1 the true start is point-lifted to\n"
              "             the consensus with that iteration's ref2cons chain (4th --iter field).\n\n"
              % (site, sampling, span_txt, TOLERANCE, TOLERANCE))
    if a.aggregate:
        with open(a.out + "-gtraj.csv.dict.txt", "w") as d:
            d.write(hdr(a.out + "-gtraj.csv") + common)
            d.write("One row per distinct trajectory of (call status, f_alt bin) across the iterations;\n"
                    "the table is the aggregated input to the G (joint-state flow) plots.\n\n"
                    "Fields:\n"
                    "  divergence, caller : condition labels (as passed on the command line)\n"
                    "  s<k>               : call status at iteration k -- called | not_called\n"
                    "  fb<k>              : f_alt bin at iteration k, f_alt = (alt reads aligned at the\n"
                    "                       SNP, exact or shifted placement) / (reads truly carrying alt);\n"
                    "                       bins 0.0-0.2 ... 0.8-1.0 (upper edge inclusive in the top bin),\n"
                    "                       or 'nocov' when no read truly carries alt\n"
                    "  n_snps             : number of SNPs that followed this trajectory\n")
        with open(a.out + "-hcomp.csv.dict.txt", "w") as d:
            d.write(hdr(a.out + "-hcomp.csv") + common)
            d.write("Read-fate composition at the SNP position, SUMMED over all SNPs, per (iteration,\n"
                    "call status). Counts are read-SNP pairs (a read covering 2 SNPs counts twice).\n\n"
                    "Fields:\n"
                    "  divergence, caller : condition labels\n"
                    "  iteration          : mapping iteration (1 = single-shot)\n"
                    "  status             : called | not_called at that iteration\n"
                    "  alt_at_p           : reads at the locus (exact/shifted) presenting the alt base at p\n"
                    "  ref_at_p           : ... presenting the reference base at p\n"
                    "  err_at_p           : ... presenting another base (sequencing error)\n"
                    "  clipped            : at the locus but soft-clipped over p\n"
                    "  deleted            : at the locus but p falls in a deletion of the alignment\n"
                    "  mapped_elsewhere   : mismapped (> %d bp from the true start) or on another contig\n"
                    "  unmapped           : unmapped reads\n"
                    "  n_reads            : all reads truly covering the SNPs (sum of the above)\n"
                    % TOLERANCE)
        return
    pile_txt = ("  depth_p            : reads aligned over p in the mapped BAM (samtools-style pileup,\n"
                "                       secondary/supplementary/duplicate/QC-fail skipped), whatever\n"
                "                       their true origin -- the caller's-eye view\n"
                "  ref_p, alt_p, oth_p: of those, reads presenting the ref / alt / another base at p\n"
                "  mapq_p             : mean mapping quality of those reads (0.0 when depth_p = 0)\n"
                "  mapq20_p           : of those, reads with MAPQ >= 20\n"
                "  mmin_p             : of those, reads MISMAPPED INTO p (true span does not cover p)\n"
                "  mmin_ref_p, mmin_alt_p : mismapped-in reads presenting the ref / alt base\n"
                ) if a.pileup else ""
    with open(a.out + ".dict.txt", "w") as d:
        d.write(hdr(a.out) + common)
        d.write("One row per (SNP, iteration): the exact read ledger of the reads that TRULY cover the\n"
                "SNP. All cells are read counts.\n\n"
                "Fields:\n"
                "  divergence, caller : condition labels\n"
                "  iteration          : mapping iteration (1 = single-shot)\n"
                "  snp_pos            : 1-based reference position of the SNP\n"
                "  ref, alt           : truth alleles\n"
                "  recovered          : 1 if recovered at this iteration (see Sites), else 0\n"
                "  n_reads            : reads truly covering the SNP (sum of the 14 ledger cells)\n"
                "  n_true_alt         : of those, reads that truly carry the alt allele\n"
                "  ex_alt ex_ref ex_err ex_clip ex_del : exact placement; at p the read presents the\n"
                "                       alt / ref / another base, or is soft-clipped over / deletes p\n"
                "  sh_alt ... sh_del  : the same for shifted placement\n"
                "  mm_clip, mm_noclip : mismapped reads, split by whether the read has any soft-clip\n"
                "  dc, un             : reads on another contig / unmapped\n"
                + pile_txt +
                "\nDerived downstream: f_alt = (ex_alt + sh_alt) / n_true_alt; mismapped-out =\n"
                "mm_clip + mm_noclip + dc + un.\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth-bam")
    ap.add_argument("--mason-vcf",
                    help="mason-format truth VCF the truth BAM was simulated from (mm39-truth.vcf.gz); "
                         "rebuilds the haplotypes that give each read's true span from its oP tag")
    ap.add_argument("--reference", help="reference FASTA the simulation used (indexed)")
    ap.add_argument("--truth-vcf")
    ap.add_argument("--region")
    ap.add_argument("--divergence", required=True)
    ap.add_argument("--caller", required=True)
    ap.add_argument("--sample", type=int, default=8000, help="approx SNPs to sample; <=0 = ALL")
    ap.add_argument("--iter", action="append",
                    metavar="ITER:MAPPED_BAM:RECOVERED_VCF")
    ap.add_argument("--out", required=True)
    ap.add_argument("--dicts-only", metavar="ORIGINAL_COMMAND",
                    help="write ONLY the data dictionaries for an already-existing output (same --out, "
                         "--aggregate/--pileup/--het, --sample, --divergence, --caller as the run that made "
                         "it), recording ORIGINAL_COMMAND as its provenance; no BAM/VCF is read. Alone: an "
                         "output made before the 2026-10-27 true-span fix; with --table-commit: one made "
                         "by the current method at that commit.")
    ap.add_argument("--table-commit", metavar="COMMIT",
                    help="with --dicts-only: the commit whose variant_read_fate.py made the output (its "
                         "previous dictionary's Git commit); the dictionary then describes the current "
                         "method. For a fix to the dictionary text that leaves the computation unchanged.")
    ap.add_argument("--het", action="store_true",
                    help="score HETEROZYGOUS (0/1) truth sites instead of hom-div (1|1). n_true_alt is "
                         "then the alt-HAPLOTYPE reads only (via oH), f_alt measures how much of that "
                         "alt support survives (the het-dropout mechanism), and 'recovered' means the "
                         "caller called the site correctly as a het (GT 0/1). Pass the per-iteration "
                         "caller VCF (e.g. mm39_calls.norm.vcf.gz) as each --iter RECOVERED_VCF.")
    ap.add_argument("--pileup", action="store_true",
                    help="per-SNP mode only: also emit the POSITION-CENTRIC pileup at p from the mapped "
                         "BAM (depth_p, ref_p, alt_p, oth_p, mapq_p, mapq20_p) -- the caller's-eye view "
                         "including reads mismapped INTO p. Lets the companion decomposition ask what "
                         "changes (ref reads? depth? MAPQ?) when a high-f_alt SNP flips to called.")
    ap.add_argument("--aggregate", action="store_true",
                    help="do NOT write per-SNP rows; instead accumulate the small summary tables "
                         "for the joint-state flow (G) and the at-p composition (H) in memory and "
                         "write <out>-gtraj.csv and <out>-hcomp.csv. Lets --sample 0 (full chromosome) "
                         "run without emitting a multi-GB per-SNP ledger.")
    a = ap.parse_args()

    # --dicts-only: document an existing output without recomputing it
    if a.table_commit and not a.dicts_only:
        ap.error("--table-commit needs --dicts-only")
    if a.dicts_only:
        write_dictionaries(a)
        return
    missing = [f for f in ("truth_bam", "mason_vcf", "reference", "truth_vcf", "region", "iter")
               if not getattr(a, f)]
    if missing:
        ap.error("missing required arguments: " + ", ".join("--" + m.replace("_", "-") for m in missing))

    # truth reads: key -> (ref_name, true start, true end, oH-haplotype, placement start). oH is
    # mason's 1-based source-haplotype tag, so we know each read's TRUE allele independent of
    # sequencing error; the true span is the read's haplotype interval (from oP) mapped back to the
    # reference. The placement start is where an aligner would lay the read down unclipped (see
    # truth_haplotypes.load_truth_reads).
    truth, (n_tr, n_sd, n_ed, n_sh) = truth_haplotypes.load_truth_reads(a.truth_bam, a.mason_vcf,
                                                                        a.reference, read_key)
    sys.stderr.write("# %d truth reads; placement start differs from truth-BAM POS for %d (%.3f%%), "
                     "true end from its CIGAR end for %d (%.3f%%); placed off oP by the leftmost-bases "
                     "alignment: %d (%.3f%%)\n"
                     % (n_tr, n_sd, 100 * n_sd / n_tr, n_ed, 100 * n_ed / n_tr, n_sh, 100 * n_sh / n_tr))

    snps = load_het_snps(a.truth_vcf) if a.het else load_hom_div_snps(a.truth_vcf)
    if a.sample and a.sample > 0:
        snps = snps[::max(1, len(snps) // a.sample)]
    pos = [p for p, _, _, _ in snps]
    refb = [x for _, x, _, _ in snps]
    altb = [x for _, _, x, _ in snps]
    althap = [h for _, _, _, h in snps]      # haplotypes carrying alt at each SNP
    npos = len(pos)

    # --- output accumulators ----------------------------------------------
    # Per-SNP mode collects one row per (SNP, iteration). Aggregate mode keeps ONLY the small
    # summaries needed for figures G and H, so --sample 0 (full chromosome) never materializes a
    # multi-GB per-SNP table.
    BINLAB = ["0.0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0"]

    def falt_bin(alt_at_p, n_true):                    # f_alt = alt aligned at p / reads truly carrying alt
        if n_true <= 0:
            return "nocov"
        b = int((alt_at_p / n_true) / 0.2)
        return BINLAB[b if b < len(BINLAB) else len(BINLAB) - 1]

    rows = []                                          # per-SNP mode
    n_iter = len(a.iter)
    snp_traj = [[] for _ in range(npos)]               # aggregate: [(status, f_alt bin) per iter] -> G
    hcomp = {}                                         # aggregate: (iter, status) -> 8 running sums -> H

    for spec in a.iter:
        parts = spec.split(":")
        it, mapped_bam, recov_vcf = parts[0], parts[1], parts[2]
        chain = parts[3] if len(parts) > 3 and parts[3] else None    # ref2cons chain for s02 iter>1
        recovered = load_recovered(recov_vcf, het=a.het)
        # s02 iter>1: reads map to the length-shifted consensus. Lift reference coordinates INTO
        # consensus coordinates so the whole per-p classification runs in the mapped BAM's system.
        lift = make_lifter(chain) if chain else None
        pmap = [lift(pp) for pp in pos] if lift else pos             # each SNP's coord in the mapped BAM
        cnt = [[0] * len(CELLS) for _ in range(npos)]
        ntrue = [0] * npos                                # reads that TRULY carry alt (oH-based)
        with pysam.AlignmentFile(mapped_bam, "rb") as mf:
            for r in mf:
                if r.is_secondary or r.is_supplementary:
                    continue
                t = truth.get(read_key(r))
                if t is None:
                    continue
                lo = bisect.bisect_left(pos, t[1])         # reads assigned by TRUE span
                hi = bisect.bisect_left(pos, t[2])
                if hi <= lo:
                    continue
                oh = t[3]                                  # this read's true source haplotype
                for j in range(lo, hi):
                    if oh in althap[j]:                    # read truly carries alt at this SNP
                        ntrue[j] += 1
                # read's TRUE (placement) start in the mapped BAM's coordinate system (identity for
                # iter1; lifted reference->consensus for s02 iter>1)
                tstart = lift(t[4]) if lift else t[4]
                has_clip = r.cigartuples is not None and any(op == 4 for op, _ in r.cigartuples)
                # ---- placement (deterministic, from the alignment) ----
                if r.is_unmapped:
                    for j in range(lo, hi):
                        cnt[j][CIX["un"]] += 1
                    continue
                if r.reference_name != t[0] or tstart is None:   # off-contig, or true start in a chain gap
                    for j in range(lo, hi):
                        cnt[j][CIX["dc"]] += 1
                    continue
                delta = abs(r.reference_start - tstart)
                if delta > TOLERANCE:
                    cell = "mm_clip" if has_clip else "mm_noclip"
                    for j in range(lo, hi):
                        cnt[j][CIX[cell]] += 1
                    continue
                # ---- at the locus (exact or shifted): resolve pos_at_p + base per SNP ----
                pref = "ex" if delta == 0 else "sh"
                r_end = r.reference_end                      # exclusive (mapped-BAM coords)
                aln = {rp: qp for qp, rp in r.get_aligned_pairs(matches_only=True)}
                seq = r.query_sequence
                for j in range(lo, hi):
                    p = pmap[j]                              # SNP coord in the mapped BAM's system
                    if p is None:                           # SNP fell in a chain gap -> cannot classify
                        continue
                    qp = aln.get(p)
                    if qp is not None:                       # p aligned -> a base is presented
                        b = seq[qp].upper()
                        suf = "alt" if b == altb[j] else ("ref" if b == refb[j] else "err")
                    elif r.reference_start <= p < r_end:     # within the alignment span but not aligned
                        suf = "del"
                    else:                                    # outside aligned span -> soft-clipped over p
                        suf = "clip"
                    cnt[j][CIX[pref + "_" + suf]] += 1

        # position-centric pileup at p (separate handle -> clean random access via the BAM index)
        pile = None
        if a.pileup and not a.aggregate:
            empty = (0, 0, 0, 0, 0.0, 0, 0, 0, 0)
            with pysam.AlignmentFile(mapped_bam, "rb") as pf:
                pile = [pileup_at(pf, a.region, pmap[j], pos[j], refb[j], altb[j], truth)
                        if pmap[j] is not None else empty for j in range(npos)]

        for j in range(npos):
            rec = 1 if (pos[j], altb[j]) in recovered else 0
            c = cnt[j]
            if a.aggregate:
                # G: this SNP's (call status, f_alt bin) at this iteration
                alt_at_p = c[CIX["ex_alt"]] + c[CIX["sh_alt"]]
                status = "called" if rec else "not_called"
                snp_traj[j].append((status, falt_bin(alt_at_p, ntrue[j])))
                # H: sum the per-p read-fate categories over all SNPs, per (iteration, status)
                h = hcomp.setdefault((it, status), [0] * 8)
                h[0] += alt_at_p                                                # alt maps at p
                h[1] += c[CIX["ex_ref"]] + c[CIX["sh_ref"]]                     # ref maps at p
                h[2] += c[CIX["ex_err"]] + c[CIX["sh_err"]]                     # other base at p (error)
                h[3] += c[CIX["ex_clip"]] + c[CIX["sh_clip"]]                   # clipped over p
                h[4] += c[CIX["ex_del"]] + c[CIX["sh_del"]]                     # deleted at p
                h[5] += c[CIX["mm_clip"]] + c[CIX["mm_noclip"]] + c[CIX["dc"]]  # mapped elsewhere (absent)
                h[6] += c[CIX["un"]]                                            # unmapped
                h[7] += sum(c)                                                 # total covering reads
            else:
                extra = list(pile[j]) if pile is not None else []
                rows.append([a.divergence, a.caller, it, pos[j] + 1, refb[j], altb[j], rec,
                             sum(c), ntrue[j]] + c + extra)

    if a.aggregate:
        write_aggregate(a, n_iter, snp_traj, hcomp)
    else:
        with open(a.out, "w") as out:
            head = ["divergence", "caller", "iteration", "snp_pos", "ref", "alt",
                    "recovered", "n_reads", "n_true_alt"] + CELLS + (PILE_COLS if a.pileup else [])
            out.write(",".join(head) + "\n")
            for row in rows:
                out.write(",".join(str(x) for x in row) + "\n")
    write_dictionaries(a)


if __name__ == "__main__":
    main()
