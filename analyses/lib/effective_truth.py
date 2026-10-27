#!/usr/bin/env python3
"""
effective_truth.py -- check that every truth SNP was actually simulated into the reads.

A "phantom" truth SNP would be one that is in the truth VCF but that mason never put into the reads,
which no caller or iteration could ever recover. For every biallelic truth SNP this tool counts, among
the truth reads that come from a haplotype carrying the alt allele, how many present the alt, and flags
a phantom when the alt is essentially absent (depth >= MIN_DEPTH and alt fraction < ALT_FRAC).

HOW A READ'S BASE AT A SNP IS FOUND (no CIGAR, no POS). Both haplotypes are rebuilt from the
mason-format truth VCF (truth_haplotypes.py), the file mason itself applied. The mason -oa truth BAM
gives each read's source haplotype (tag oH) and 0-based leftmost position on it (tag oP), so the base
at a SNP whose haplotype position is h is SEQ[h - oP]. Reads are scored only when they are CONSISTENT
with their haplotype: their mismatches against it, counted away from that haplotype's truth SNPs, do
not exceed mason's own count of sequencing errors (tag XE). Mismatches AT truth SNPs are deliberately
not counted, so a phantom (the read carries the reference base at the SNP) cannot exclude the read that
reveals it. Inconsistent reads are tested for one 1-2 bp shift (a mason sequencing indel); the rest are
"unexplained", and the tool fails if they exceed MAX_UNEXPLAINED of the reads, or MAX_WINDOW_UNEXPLAINED
of any WINDOW-bp window, because that would mean the rebuilt haplotypes differ from mason's.

WHY (2026-10-26). The first version of this tool read bases through the -oa CIGAR. mason writes dense
SNP clusters there as fake 1D..1I gap pairs although the reads carry SNPs only, so bases landed shifted
and real alt sites looked all-reference: tens of thousands of SNPs were wrongly flagged. That method is
replayed here in the same pass (n_cigar_flag_* columns) so its counts stay reproducible, and so does
how many of its flags are phantoms under the corrected method.

Outputs (<out-prefix>-...): phantom.tsv (one row per phantom SNP) and summary.csv (one row), each with a
data dictionary. Usage (from the project root):
  python analyses/lib/effective_truth.py --mason-vcf .../mm39-truth.vcf.gz --reference .../mm39-19.fa \\
      --truth-bam .../truth.oa.bam --region 19 --divergence 0.10 --seed 42 --out-prefix .../eff-0.10
"""
import argparse
import os

import numpy as np
import pysam

import provenance                       # shared provenance block for the data dictionaries
import truth_haplotypes as th           # the two haplotypes, rebuilt as mason applied them

# --- thresholds ----------------------------------------------------------------------------------
MIN_DEPTH = 5                  # consistent alt-haplotype reads needed to judge a site (else low coverage)
ALT_FRAC = 0.10                # alt fraction below which a judged site is a phantom
CHUNK = 200_000                # reads processed together (numpy)
MAX_SHIFT = 2                  # sequencing-indel shifts tried for inconsistent reads (+-1..MAX_SHIFT bp)
MAX_UNEXPLAINED = 0.001        # fail if more than this fraction of reads is unexplained
WINDOW = 10_000                # haplotype window (bp, by oP) for the local unexplained check
MIN_WINDOW_READS = 100         # windows with fewer reads are not judged
MAX_WINDOW_UNEXPLAINED = 0.05  # fail if any judged window has more unexplained reads than this

# Written into every dictionary of this tool's outputs (and of tables combined from them).
CORRECTION_NOTE = ("Corrected method (2026-10-26). Earlier versions of this table came from reading "
                   "bases through mason's -oa CIGAR, which writes dense SNP clusters as fake 1D..1I gap "
                   "pairs; that wrongly flagged tens of thousands of SNPs as phantoms (still counted here "
                   "as n_cigar_flag_*). Bases are now read at each read's haplotype position (oP).")


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mason-vcf", required=True, help="mason-format truth VCF mason applied (mm39-truth.vcf.gz)")
    ap.add_argument("--reference", required=True, help="region FASTA the simulation used")
    ap.add_argument("--truth-bam", required=True, help="mason -oa truth BAM (oH, oP, XE tags)")
    ap.add_argument("--region", required=True)
    ap.add_argument("--divergence", required=True)
    ap.add_argument("--seed", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--dicts-only", metavar="ORIGINAL_COMMAND",
                    help="write ONLY the data dictionaries for existing outputs (same arguments as the run "
                         "that made them), recording ORIGINAL_COMMAND as provenance")
    ap.add_argument("--dict-note", action="append", default=[],
                    help="extra provenance note for the dictionaries (repeatable)")
    return ap.parse_args()


# --- data dictionaries ---------------------------------------------------------------------------

def write_dictionaries(a):
    def hdr(table):
        notes = [CORRECTION_NOTE] + list(a.dict_note)
        if a.dicts_only:
            return provenance.header(table, command=a.dicts_only, backfill=True, notes=notes)
        return provenance.header(table, notes=notes)

    method = ("Method: both haplotypes rebuilt from %s (the mason-format truth VCF); each truth read's base\n"
              "at a SNP is SEQ[h - oP] (h = the SNP's position on the read's haplotype oH; oP = the read's\n"
              "haplotype start, from %s). A site is judged from CONSISTENT reads of the haplotypes that\n"
              "carry its alt (both for 1|1, one for a het): reads whose mismatches against their haplotype,\n"
              "away from that haplotype's truth SNPs, do not exceed mason's error count XE. Phantom = depth\n"
              ">= %d and alt fraction < %.2f; depth < %d = low coverage (not judged).\n\n"
              % (a.mason_vcf, a.truth_bam, MIN_DEPTH, ALT_FRAC, MIN_DEPTH))
    ph = a.out_prefix + "-phantom.tsv"
    with open(ph + ".dict.txt", "w") as d:
        d.write(hdr(ph) + method)
        d.write("One row per phantom SNP (tab-separated; header only when there are none).\n\nFields:\n"
                "  pos         : 1-based reference position\n"
                "  ref, alt    : truth alleles\n"
                "  zyg         : div (1|1) | het (0|1 or 1|0)\n"
                "  depth       : consistent reads from the alt-carrying haplotype(s) presenting a base here\n"
                "  alt_count   : of those, presenting the alt\n"
                "  alt_frac    : alt_count / depth\n"
                "  ref_count, other_count : of those, presenting the reference / another base\n"
                "  alt_haps    : haplotypes carrying the alt (1, 2 or 1,2)\n"
                "  hap1_pos, hap2_pos : 1-based position of the SNP on each haplotype\n"
                "  depth_all, alt_count_all : the same counts over ALL reads of those haplotypes\n")
    sm = a.out_prefix + "-summary.csv"
    with open(sm + ".dict.txt", "w") as d:
        d.write(hdr(sm) + method)
        d.write("One row (this region / divergence / seed). The first 12 columns keep the meaning they\n"
                "had in earlier versions; rates are NA when their denominator is 0.\n\nFields:\n"
                "  region, divergence, seed : condition labels\n"
                "  n_truth_snp       : biallelic truth SNPs (div + het)\n"
                "  n_truth_div, n_truth_het : of those, 1|1 / heterozygous\n"
                "  n_phantom_div, n_phantom_het : phantom SNPs, by zygosity\n"
                "  phantom_div_rate  : n_phantom_div / n_truth_div\n"
                "  phantom_het_rate  : n_phantom_het / n_truth_het\n"
                "  n_lowcov          : SNPs with depth < %d (not judged)\n"
                "  mean_alt_frac     : mean alt fraction over the judged SNPs\n"
                "  n_lowcov_div, n_lowcov_het : n_lowcov by zygosity\n"
                "  min_alt_frac      : lowest alt fraction of any judged SNP\n"
                "  n_af_0.10_0.50, n_af_0.50_0.90, n_af_ge_0.90 : judged SNPs by alt-fraction bin\n"
                "  refhap_alt_frac   : at het SNPs, alt fraction among consistent reads of the haplotype\n"
                "                      WITHOUT the alt (should be about the sequencing-error rate)\n"
                "  n_reads           : truth reads scanned\n"
                "  frac_reads_exact  : reads whose mismatches against their haplotype equal XE exactly\n"
                "  frac_reads_consistent : reads used for scoring (see Method)\n"
                "  frac_reads_seq_indel  : inconsistent reads explained by one 1-%d bp shift (a mason\n"
                "                      sequencing indel); not scored\n"
                "  frac_reads_unexplained : the remaining inconsistent reads (tool fails above %.4f)\n"
                "  max_window_unexplained : highest unexplained fraction of any %d-bp haplotype window\n"
                "                      with >= %d reads (tool fails above %.2f)\n"
                "  n_cigar_flag_div, n_cigar_flag_het : SNPs the OLD method (bases read through the -oa\n"
                "                      CIGAR, all reads pooled) flags with the same thresholds\n"
                "  n_cigar_flag_now_phantom : of those, phantom under the corrected method\n"
                "  hap1_len, hap2_len : rebuilt haplotype lengths (checked against ref + inserted - deleted)\n"
                % (MIN_DEPTH, MAX_SHIFT, MAX_UNEXPLAINED, WINDOW, MIN_WINDOW_READS, MAX_WINDOW_UNEXPLAINED))


# --- scanning ------------------------------------------------------------------------------------

class Tally:
    """Per-SNP counters, all int64 arrays of length n_snp."""

    def __init__(self, n):
        z = lambda: np.zeros(n, dtype=np.int64)
        self.alt, self.ref, self.oth = z(), z(), z()            # consistent reads, alt-carrying haps
        self.all_alt, self.all_dep = z(), z()                   # all reads, alt-carrying haps
        self.rh_alt, self.rh_dep = z(), z()                     # consistent reads, the het's other hap
        self.old_alt, self.old_ref, self.old_oth = z(), z(), z()  # old CIGAR method, all reads pooled


def add(counter, sites, weights, n):
    counter += np.bincount(sites, weights=weights, minlength=n).astype(np.int64)


def score_hap_chunk(Q, OP, XE, hap, site, snp_alt, snp_ref, t, n_snp, cls, win):
    """Score one chunk of reads from one haplotype (Q: reads x L bases; OP: haplotype starts)."""
    L = Q.shape[1]
    idx = OP[:, None] + np.arange(L)
    if idx.max() >= len(hap.seq):
        raise ValueError("a read runs past the end of its rebuilt haplotype (oP + length > %d)" % len(hap.seq))
    H = hap.seq[idx]
    S = site[idx]
    mm = Q != H
    off = mm & (S < 0)                                   # mismatches away from this hap's truth SNPs
    mm_off = off.sum(1)
    exact = mm.sum(1) == XE
    consistent = mm_off <= XE
    # inconsistent reads: one sequencing indel (prefix unshifted + suffix shifted) explains them?
    bad = np.where(~consistent)[0]
    explained = np.zeros(len(bad), dtype=bool)
    if len(bad):
        pre = np.concatenate([np.zeros((len(bad), 1), dtype=np.int32),
                              np.cumsum(off[bad], axis=1, dtype=np.int32)[:, :-1]], axis=1)
        for s in [d for k in range(1, MAX_SHIFT + 1) for d in (-k, k)]:
            j = np.clip(idx[bad] + s, 0, len(hap.seq) - 1)
            offs = (Q[bad] != hap.seq[j]) & (site[j] < 0)
            suf = np.cumsum(offs[:, ::-1], axis=1, dtype=np.int32)[:, ::-1]
            explained |= (pre + suf).min(axis=1) <= XE[bad]
    seq_indel = np.zeros(len(Q), dtype=bool); seq_indel[bad[explained]] = True
    unexpl = ~consistent & ~seq_indel
    cls["n"] += len(Q); cls["exact"] += int(exact.sum()); cls["consistent"] += int(consistent.sum())
    cls["seq_indel"] += int(seq_indel.sum()); cls["unexplained"] += int(unexpl.sum())
    w = OP // WINDOW
    np.add.at(win[0], w, 1); np.add.at(win[1], w[unexpl], 1)

    # alt-carrying sites: consistent reads, and all reads
    m = S >= 0
    for rows, counters in ((consistent, ("alt", "ref", "oth")), (None, ("all_alt", None, None))):
        mk = m if rows is None else (m & rows[:, None])
        s = S[mk]; b = Q[mk]
        is_alt = b == snp_alt[s]
        add(getattr(t, counters[0]), s, is_alt, n_snp)
        if rows is None:
            add(t.all_dep, s, None, n_snp)
        else:
            is_ref = b == snp_ref[s]
            add(t.ref, s, is_ref, n_snp); add(t.oth, s, ~is_alt & ~is_ref, n_snp)
    # het sites on the haplotype WITHOUT the alt (sanity: alt there should be ~ the error rate)
    mk = (S <= -2) & consistent[:, None]
    s = -(S[mk] + 2)
    add(t.rh_alt, s, Q[mk] == snp_alt[s], n_snp); add(t.rh_dep, s, None, n_snp)


def old_method_150m(Q, RS, refsite, snp_alt, snp_ref, t, n_snp):
    """The old CIGAR method for 150M reads (aligned base i sits at reference POS-1 + i)."""
    s = refsite[RS[:, None] + np.arange(Q.shape[1])]
    m = s >= 0
    s = s[m]; b = Q[m]
    is_alt = b == snp_alt[s]; is_ref = b == snp_ref[s]
    add(t.old_alt, s, is_alt, n_snp); add(t.old_ref, s, is_ref, n_snp)
    add(t.old_oth, s, ~is_alt & ~is_ref, n_snp)


def main():
    a = parse_args()
    if a.dicts_only:
        write_dictionaries(a)
        return

    # --- haplotypes and the truth SNPs -----------------------------------------------------------
    ref_seq = pysam.FastaFile(a.reference).fetch(a.region)
    v = th.load_mason_vcf(a.mason_vcf, a.region)
    haps = th.build(ref_seq, v)
    for k in (0, 1):
        exp = th.expected_length(len(ref_seq), v, k)
        if len(haps[k].seq) != exp:
            raise ValueError("haplotype %d is %d bp, expected %d" % (k + 1, len(haps[k].seq), exp))
    snp = np.where(v.is_snp)[0]                          # variant index of each truth SNP
    n_snp = len(snp)
    snp_ref = np.frombuffer("".join(v.ref[i] for i in snp).encode(), dtype=np.uint8)
    snp_alt = np.frombuffer("".join(v.alt[i] for i in snp).encode(), dtype=np.uint8)
    gt = v.gt[snp]                                       # (n_snp, 2) alleles on hap 1, hap 2
    is_div = (gt[:, 0] == 1) & (gt[:, 1] == 1)
    # per-haplotype site map: >= 0 = SNP whose alt is on this hap; <= -2 = het SNP with the ref here
    sites = []
    for k in (0, 1):
        site = np.full(len(haps[k].seq), -1, dtype=np.int32)
        hp = haps[k].var_pos[snp]
        on = gt[:, k] == 1
        site[hp[on]] = np.where(on)[0]
        site[hp[~on]] = -(np.where(~on)[0] + 2)
        sites.append(site)
    refsite = np.full(len(ref_seq), -1, dtype=np.int32)  # the old method works in reference coordinates
    refsite[v.pos[snp]] = np.arange(n_snp)

    # --- one pass over the truth BAM -------------------------------------------------------------
    t = Tally(n_snp)
    cls = {"n": 0, "exact": 0, "consistent": 0, "seq_indel": 0, "unexplained": 0}
    wins = [[np.zeros(len(h.seq) // WINDOW + 1, dtype=np.int64) for _ in range(2)] for h in haps]
    buf = {0: [], 1: []}                                  # per haplotype: (seq, oP, XE)
    old_buf = []                                          # 150M reads for the old method: (seq, POS-1)
    L = None

    def flush(k):
        if not buf[k]:
            return
        Q = np.frombuffer(b"".join(x[0] for x in buf[k]), dtype=np.uint8).reshape(len(buf[k]), L)
        OP = np.array([x[1] for x in buf[k]], dtype=np.int64)
        XE = np.array([x[2] for x in buf[k]], dtype=np.int64)
        score_hap_chunk(Q, OP, XE, haps[k], sites[k], snp_alt, snp_ref, t, n_snp, cls, wins[k])
        buf[k] = []

    def flush_old():
        if not old_buf:
            return
        Q = np.frombuffer(b"".join(x[0] for x in old_buf), dtype=np.uint8).reshape(len(old_buf), L)
        RS = np.array([x[1] for x in old_buf], dtype=np.int64)
        old_method_150m(Q, RS, refsite, snp_alt, snp_ref, t, n_snp)
        old_buf.clear()

    with pysam.AlignmentFile(a.truth_bam, "rb") as bam:
        for r in bam:
            if r.is_secondary or r.is_supplementary or r.is_unmapped:
                raise ValueError("truth BAM read %s is not a primary mapped record" % r.query_name)
            if r.reference_name != a.region:
                raise ValueError("truth BAM read %s is on %s, expected %s" % (r.query_name, r.reference_name, a.region))
            seq = r.query_sequence.upper().encode()
            if L is None:
                L = len(seq)
            elif len(seq) != L:
                raise ValueError("read %s is %d bp; this tool expects one read length (%d)" % (r.query_name, len(seq), L))
            k = r.get_tag("oH") - 1
            if k not in (0, 1):
                raise ValueError("read %s has oH %d" % (r.query_name, k + 1))
            buf[k].append((seq, r.get_tag("oP"), r.get_tag("XE")))
            if len(buf[k]) >= CHUNK:
                flush(k)
            # the OLD method: bases read through the -oa CIGAR at reference coordinates
            if r.cigarstring == "%dM" % L:
                old_buf.append((seq, r.reference_start))
                if len(old_buf) >= CHUNK:
                    flush_old()
            else:
                pr = np.array(r.get_aligned_pairs(matches_only=True), dtype=np.int64)
                s = refsite[pr[:, 1]]
                m = s >= 0
                if m.any():
                    q = np.frombuffer(seq, dtype=np.uint8)[pr[m, 0]]; s = s[m]
                    is_alt = q == snp_alt[s]; is_ref = q == snp_ref[s]
                    np.add.at(t.old_alt, s, is_alt); np.add.at(t.old_ref, s, is_ref)
                    np.add.at(t.old_oth, s, ~is_alt & ~is_ref)
    flush(0); flush(1); flush_old()

    # --- read-level checks (fail loudly if the haplotypes do not match mason's) -------------------
    n = cls["n"]
    if n == 0:
        raise ValueError("no reads in %s" % a.truth_bam)
    frac_unexpl = cls["unexplained"] / n
    win_reads = np.concatenate([w[0] for w in wins]); win_unex = np.concatenate([w[1] for w in wins])
    judged = win_reads >= MIN_WINDOW_READS
    max_win = float((win_unex[judged] / win_reads[judged]).max()) if judged.any() else 0.0
    if frac_unexpl > MAX_UNEXPLAINED:
        raise ValueError("%.4f of reads are unexplained (> %.4f): the rebuilt haplotypes do not match "
                         "mason's" % (frac_unexpl, MAX_UNEXPLAINED))
    if max_win > MAX_WINDOW_UNEXPLAINED:
        raise ValueError("a %d-bp window has %.3f unexplained reads (> %.2f)" % (WINDOW, max_win, MAX_WINDOW_UNEXPLAINED))

    # --- per-site calls ---------------------------------------------------------------------------
    dep = t.alt + t.ref + t.oth
    judged_site = dep >= MIN_DEPTH
    af = np.where(dep > 0, t.alt / np.maximum(dep, 1), np.nan)
    phantom = judged_site & (af < ALT_FRAC)
    old_dep = t.old_alt + t.old_ref + t.old_oth
    old_flag = (old_dep >= MIN_DEPTH) & (t.old_alt < ALT_FRAC * old_dep)

    with open(a.out_prefix + "-phantom.tsv", "w") as out:
        out.write("pos\tref\talt\tzyg\tdepth\talt_count\talt_frac\tref_count\tother_count\talt_haps\t"
                  "hap1_pos\thap2_pos\tdepth_all\talt_count_all\n")
        for j in np.where(phantom)[0]:
            i = snp[j]
            haps_on = ",".join(str(k + 1) for k in (0, 1) if gt[j, k] == 1)
            out.write("%d\t%s\t%s\t%s\t%d\t%d\t%.3f\t%d\t%d\t%s\t%d\t%d\t%d\t%d\n" % (
                v.pos[i] + 1, v.ref[i], v.alt[i], "div" if is_div[j] else "het", dep[j], t.alt[j],
                af[j], t.ref[j], t.oth[j], haps_on, haps[0].var_pos[i] + 1, haps[1].var_pos[i] + 1,
                t.all_dep[j], t.all_alt[j]))

    n_div, n_het = int(is_div.sum()), int((~is_div).sum())
    rate = lambda x, d: "%.6f" % (x / d) if d else "NA"
    afj = af[judged_site]
    het = ~is_div
    cols = [("region", a.region), ("divergence", a.divergence), ("seed", a.seed),
            ("n_truth_snp", n_snp), ("n_truth_div", n_div), ("n_truth_het", n_het),
            ("n_phantom_div", int((phantom & is_div).sum())), ("n_phantom_het", int((phantom & het).sum())),
            ("phantom_div_rate", rate(int((phantom & is_div).sum()), n_div)),
            ("phantom_het_rate", rate(int((phantom & het).sum()), n_het)),
            ("n_lowcov", int((~judged_site).sum())),
            ("mean_alt_frac", "%.4f" % afj.mean() if len(afj) else "NA"),
            ("n_lowcov_div", int((~judged_site & is_div).sum())), ("n_lowcov_het", int((~judged_site & het).sum())),
            ("min_alt_frac", "%.4f" % afj.min() if len(afj) else "NA"),
            ("n_af_0.10_0.50", int(((afj >= 0.10) & (afj < 0.50)).sum())),
            ("n_af_0.50_0.90", int(((afj >= 0.50) & (afj < 0.90)).sum())),
            ("n_af_ge_0.90", int((afj >= 0.90).sum())),
            ("refhap_alt_frac", "%.5f" % (t.rh_alt.sum() / t.rh_dep.sum()) if t.rh_dep.sum() else "NA"),
            ("n_reads", n), ("frac_reads_exact", "%.5f" % (cls["exact"] / n)),
            ("frac_reads_consistent", "%.5f" % (cls["consistent"] / n)),
            ("frac_reads_seq_indel", "%.5f" % (cls["seq_indel"] / n)),
            ("frac_reads_unexplained", "%.6f" % frac_unexpl), ("max_window_unexplained", "%.4f" % max_win),
            ("n_cigar_flag_div", int((old_flag & is_div).sum())), ("n_cigar_flag_het", int((old_flag & het).sum())),
            ("n_cigar_flag_now_phantom", int((old_flag & phantom).sum())),
            ("hap1_len", len(haps[0].seq)), ("hap2_len", len(haps[1].seq))]
    with open(a.out_prefix + "-summary.csv", "w") as out:
        out.write(",".join(c for c, _ in cols) + "\n")
        out.write(",".join(str(x) for _, x in cols) + "\n")
    write_dictionaries(a)
    print("D=%s: %d SNPs, phantom div %s / het %s, min alt fraction %s; reads %d (consistent %.4f, "
          "unexplained %.6f); old CIGAR flags div %s / het %s"
          % (a.divergence, n_snp, cols[6][1], cols[7][1], cols[14][1], n, cls["consistent"] / n,
             frac_unexpl, cols[25][1], cols[26][1]))


if __name__ == "__main__":
    main()
