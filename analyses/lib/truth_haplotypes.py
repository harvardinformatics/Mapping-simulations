#!/usr/bin/env python3
"""
truth_haplotypes.py -- the two simulated haplotypes, rebuilt exactly as mason_simulator applied them.

`mason_simulator -iv <truth.vcf>` builds two haplotypes from the reference and the phased genotypes of
the MASON-FORMAT truth VCF (mm39-truth.vcf.gz) and samples reads from them. Its -oa truth BAM records,
for every read, the source haplotype (tag oH; 1 = the LEFT allele of the GT) and the 0-based leftmost
position ON THAT HAPLOTYPE (tag oP). Rebuilding the same haplotypes lets a read's bases be compared with
its haplotype at oP directly, using neither the -oa CIGAR (mason writes dense SNP clusters there as fake
1D..1I gap pairs, so bases read through it land shifted) nor POS (shifted by 1-3 bp by leading gap pairs
for ~1% of reads).

Use the mason-format VCF, NOT mm39-truth.norm.vcf.gz: normalisation left-shifts indels across adjacent
SNPs, so the norm VCF has overlapping records and does not describe mason's haplotypes. The mason VCF
writes the genotype straight into the sample column with FORMAT '.', so it is parsed here rather than
with pysam/bcftools.

Library use:
    import truth_haplotypes as th
    v = th.load_mason_vcf(path, contig)          # sorted, disjoint, phased variants (numpy arrays)
    haps = th.build(ref_seq, v)                  # [Haplotype for oH = 1, Haplotype for oH = 2]
    haps[0].seq            # haplotype sequence, np.uint8 (ASCII, upper case)
    haps[0].var_pos[i]     # haplotype position of variant i's first base (its anchor for indels);
                           # for a variant NOT on this haplotype, the position of its first REF base
    haps[0].hap_to_ref(h)  # reference position of haplotype position h, or -1 (an inserted base)
    haps[0].ref_span(h0, n)  # reference span [start, end) of haplotype intervals [h0, h0 + n),
                             # plus the number of inserted bases before each interval's first
                             # reference base
    th.anchor_offsets(hap_seq, h0, bases)  # where an aligner would start reads (offset from h0)
    truth, diffs = th.load_truth_reads(bam, vcf, fasta, read_key)   # each truth read's TRUE span
                                                                     # and its placement start
Every check fails loudly (ValueError); nothing is skipped silently.
"""
import gzip

import numpy as np

GT_ALLELES = {"1|1": (1, 1), "0|1": (0, 1), "1|0": (1, 0)}   # phased genotypes the simulations write

# --- placement start: where an aligner would start the read (see anchor_offsets) ------------------
# A read with a sequencing indel among its leftmost bases is laid down by an aligner 1-3 bp off oP
# (the indel absorbed as a shift plus mismatches). Its placement start comes from aligning its
# leftmost ANCHOR_LEN bases to its own haplotype with bwa mem's default scores.
ANCHOR_LEN   = 30     # leftmost read bases aligned to the haplotype
ANCHOR_SHIFT = 3      # candidate starts on the haplotype: oP - 3 .. oP + 3
MATCH, MISMATCH, N_SCORE = 1, -4, -1   # bwa mem -A 1 -B 4; a base vs N scores -1 (bwa's ambiguous score)
GAP_OPEN, GAP_EXTEND = 6, 1            # bwa mem -O 6 -E 1: a gap of length l costs 6 + l
ANCHOR_CHUNK = 200_000                 # reads per vectorized alignment batch


class Variants:
    """The truth variants of one contig, as parallel numpy arrays in position order."""

    def __init__(self, pos, ref, alt, gt, ids):
        self.pos = np.asarray(pos, dtype=np.int64)    # 0-based position of REF's first base
        self.ref = ref                                 # list of str (upper case)
        self.alt = alt                                 # list of str (upper case)
        self.gt = np.asarray(gt, dtype=np.int8)        # shape (n, 2): allele on hap 1, hap 2
        self.ids = ids                                 # list of str (VCF ID column)
        self.is_snp = np.array([len(r) == 1 and len(a) == 1 for r, a in zip(ref, alt)])

    def __len__(self):
        return len(self.pos)


def load_mason_vcf(path, contig):
    """Read the mason-format truth VCF for one contig. Fails on anything the haplotype builder cannot
    represent exactly: unsorted or overlapping records, an unphased or unexpected genotype, a
    multiallelic record, or a variant that is neither a SNP nor an anchored insertion/deletion."""
    pos, ref, alt, gt, ids = [], [], [], [], []
    last_end = -1                                   # 0-based end (exclusive) of the previous REF span
    with gzip.open(path, "rt") if path.endswith(".gz") else open(path) as fh:
        for n, line in enumerate(fh, 1):
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if c[0] != contig:
                raise ValueError("%s line %d: contig %s, expected only %s" % (path, n, c[0], contig))
            p, r, a, g = int(c[1]) - 1, c[3].upper(), c[4].upper(), c[9].split(":")[0]
            if g not in GT_ALLELES:
                raise ValueError("%s line %d: genotype %r (expected one of %s)"
                                 % (path, n, g, ", ".join(GT_ALLELES)))
            if "," in a:
                raise ValueError("%s line %d: multiallelic record" % (path, n))
            snp = len(r) == 1 and len(a) == 1
            indel = (len(r) == 1) != (len(a) == 1) and r[0] == a[0]     # anchored insertion/deletion
            if not (snp or indel):
                raise ValueError("%s line %d: %s>%s is neither a SNP nor an anchored indel" % (path, n, r, a))
            if p < last_end:
                raise ValueError("%s line %d: position %d overlaps or precedes the previous record "
                                 "(which ends at %d)" % (path, n, p + 1, last_end))
            last_end = p + len(r)
            pos.append(p); ref.append(r); alt.append(a); gt.append(GT_ALLELES[g]); ids.append(c[2])
    return Variants(pos, ref, alt, gt, ids)


class Haplotype:
    """One rebuilt haplotype: its sequence, where every truth variant sits on it, and a map back to
    reference coordinates (blocks of haplotype positions that copy reference positions 1:1)."""

    def __init__(self, seq, var_pos, blk_hstart, blk_rstart, blk_len):
        self.seq = seq
        self.var_pos = var_pos
        self._hs, self._rs, self._len = blk_hstart, blk_rstart, blk_len

    def hap_to_ref(self, h):
        """Reference position of haplotype position(s) h; -1 for an inserted base (no reference
        position). Accepts a scalar or a numpy array."""
        h = np.asarray(h, dtype=np.int64)
        i = np.searchsorted(self._hs, h, side="right") - 1
        ok = (i >= 0) & (h < self._hs[np.maximum(i, 0)] + self._len[np.maximum(i, 0)])
        out = np.where(ok, self._rs[np.maximum(i, 0)] + (h - self._hs[np.maximum(i, 0)]), -1)
        return out if out.ndim else int(out)

    def ref_span(self, h0, n):
        """Reference span [start, end) of the haplotype interval [h0, h0 + n) (numpy arrays), and
        lead = how many inserted bases come before the interval's first reference base (0 unless the
        interval starts inside an insertion).
        start = reference position of the interval's first base that is not an inserted base;
        end = one past the reference position of its last such base. Reference bases DELETED on
        this haplotype inside the interval lie inside the span (a read crossing a deletion spans it
        on the reference). Fails if an interval holds no reference base at all."""
        hs, rs, ln = self._hs, self._rs, self._len
        h0 = np.asarray(h0, dtype=np.int64)
        h1 = h0 + np.asarray(n, dtype=np.int64) - 1                # last haplotype base of the interval
        # first base: inside block i -> its reference position; inside the insertion that follows
        # block i -> the first reference position of block i + 1
        i = np.searchsorted(hs, h0, side="right") - 1
        in_blk = h0 < hs[i] + ln[i]
        if np.any(~in_blk & (i + 1 >= len(hs))):
            raise ValueError("an interval starts in an insertion at the haplotype end (no reference base)")
        nxt = np.minimum(i + 1, len(hs) - 1)
        start = np.where(in_blk, rs[i] + (h0 - hs[i]), rs[nxt])
        lead = np.where(in_blk, 0, hs[nxt] - h0)
        # last base: inside block j -> one past its reference position; inside the insertion that
        # follows block j -> the end of block j
        j = np.searchsorted(hs, h1, side="right") - 1
        end = np.where(h1 < hs[j] + ln[j], rs[j] + (h1 - hs[j]) + 1, rs[j] + ln[j])
        if np.any(end <= start):
            raise ValueError("%d intervals hold no reference base (entirely inside an insertion)"
                             % int(np.sum(end <= start)))
        return start, end, lead


def build(ref_seq, v):
    """Both haplotypes from the reference sequence (str or bytes, one contig) and the variants.
    Haplotype k applies every variant whose allele k is 1, exactly once, in position order -- what
    mason does when the variants are disjoint (checked by load_mason_vcf)."""
    ref = np.frombuffer(ref_seq.upper().encode() if isinstance(ref_seq, str) else ref_seq.upper(),
                        dtype=np.uint8)
    for i in range(len(v)):                          # REF must match the reference it was written on
        p, r = v.pos[i], v.ref[i]
        if ref[p:p + len(r)].tobytes().decode() != r:
            raise ValueError("variant %s at %d: REF %s does not match the reference (%s)"
                             % (v.ids[i], p + 1, r, ref[p:p + len(r)].tobytes().decode()))
    haps = []
    for k in (0, 1):
        pieces, hs, rs, ln = [], [], [], []          # sequence pieces; 1:1 reference blocks
        var_pos = np.empty(len(v), dtype=np.int64)
        last, cur = 0, 0                             # next reference position to copy; haplotype length
        for i in range(len(v)):
            p, r, a = int(v.pos[i]), v.ref[i], v.alt[i]
            # reference bases up to this variant are copied unchanged
            if p > last:
                pieces.append(ref[last:p]); hs.append(cur); rs.append(last); ln.append(p - last)
                cur += p - last
                last = p
            var_pos[i] = cur                         # the variant's first base (its anchor for indels)
            if v.gt[i, k] == 1:
                pieces.append(np.frombuffer(a.encode(), dtype=np.uint8))
                # the first base of REF and ALT correspond (SNP: the substituted base; indel: anchor)
                hs.append(cur); rs.append(p); ln.append(1)
                cur += len(a)
                last = p + len(r)
        if last < len(ref):
            pieces.append(ref[last:]); hs.append(cur); rs.append(last); ln.append(len(ref) - last)
        seq = np.concatenate(pieces) if pieces else np.empty(0, dtype=np.uint8)
        haps.append(Haplotype(seq, var_pos, np.array(hs, dtype=np.int64),
                              np.array(rs, dtype=np.int64), np.array(ln, dtype=np.int64)))
    return haps


def expected_length(ref_len, v, k):
    """Haplotype k's length from the variants alone: reference + inserted - deleted bases."""
    on = v.gt[:, k] == 1
    delta = sum(len(a) - len(r) for a, r, o in zip(v.alt, v.ref, on) if o)
    return ref_len + delta


def anchor_offsets(hap_seq, h0, lead_bases):
    """Where an aligner would start each read on its haplotype, as an offset c (-ANCHOR_SHIFT ..
    ANCHOR_SHIFT) from the read's haplotype start h0 (its oP tag).

    The read's leftmost ANCHOR_LEN bases (lead_bases: uint8 ASCII, shape (n, ANCHOR_LEN)) are aligned
    to its haplotype (hap_seq) with bwa mem's default scores: MATCH / MISMATCH / N_SCORE per base, and
    GAP_OPEN + l * GAP_EXTEND for a gap of length l. The read's first base must be aligned (no clipping,
    no leading gap) and must sit within ANCHOR_SHIFT bases of h0; the alignment may end anywhere on
    the haplotype. c is the start of the best-scoring alignment. Ties go to the smaller |c|, then to
    c > 0 over c < 0, so c = 0 unless a shifted start scores strictly better than starting at h0.
    Example: a read whose first base is a sequencing insertion aligns from h0 - 1 (the inserted base
    over the haplotype base before h0, one mismatch at most) rather than with a gap at h0.

    A read whose leftmost bases all match the haplotype at h0 scores the maximum possible there, so
    only reads with at least one mismatch are aligned (an exact shortcut: c = 0 for all others)."""
    K, W = ANCHOR_LEN, ANCHOR_SHIFT
    L = K + 3 * W                       # window on the haplotype: [h0 - W, h0 + K + 2W)
    SC = 16                             # scores are scaled by SC; the low bits hold the tie-break rank
    if 2 * W + 1 >= SC:
        raise ValueError("ANCHOR_SHIFT %d too large for the tie-break encoding" % W)
    NEG = -(1 << 40)                    # "impossible" (far below any reachable score)
    h0 = np.asarray(h0, dtype=np.int64)
    if lead_bases.shape != (len(h0), K):
        raise ValueError("lead_bases shape %s, expected (%d, %d)" % (lead_bases.shape, len(h0), K))
    # pad the haplotype so every window exists; pad bytes (0) are no base: a read cannot START on one,
    # and aligned against one a base scores as a mismatch
    hp = np.concatenate([np.zeros(W, np.uint8), hap_seq, np.zeros(K + 2 * W, np.uint8)])
    col = np.arange(L)
    c_of_col = np.arange(2 * W + 1) - W                      # start column j -> offset c = j - W
    rank = 2 * np.abs(c_of_col) + (c_of_col < 0)            # tie-break rank: 0 for c = 0, then 1, -1, 2, ...
    off = np.zeros(len(h0), dtype=np.int64)

    def score(read_base, win):
        """per-base score of read bases (m, 1) against window bases (m, L), scaled by SC"""
        s = np.where(read_base == win, MATCH, MISMATCH)
        s = np.where((read_base == ord("N")) | (win == ord("N")), N_SCORE, s)
        return SC * s.astype(np.int64)

    def deletions(M):
        """best score ending in a deletion (haplotype bases skipped) after the row's aligned cells:
        Y[j] = max over k < j of M[k] - SC * (GAP_OPEN + (j - k) * GAP_EXTEND)"""
        A = np.maximum.accumulate(M + SC * GAP_EXTEND * col, axis=1)
        Y = np.full_like(M, NEG)
        Y[:, 1:] = A[:, :-1] - SC * (GAP_OPEN + GAP_EXTEND * col[1:])
        return Y

    for s0 in range(0, len(h0), ANCHOR_CHUNK):
        sl = slice(s0, s0 + ANCHOR_CHUNK)
        win = hp[h0[sl, None] + col]                         # hp index h0 + j = haplotype h0 - W + j
        rd = lead_bases[sl]
        need = np.nonzero(np.any(win[:, W:W + K] != rd, axis=1))[0]
        if len(need) == 0:
            continue
        win, rd = win[need], rd[need]
        # row 0: the read's first base is aligned at a start column (c = -W .. W) holding a real base
        M = np.full(win.shape, NEG, dtype=np.int64)
        first = score(rd[:, :1], win)[:, :2 * W + 1] - rank
        M[:, :2 * W + 1] = np.where(win[:, :2 * W + 1] != 0, first, NEG)
        X = np.full_like(M, NEG)                              # no leading insertion (that is a clip)
        Y = deletions(M)
        for i in range(1, K):
            prev = np.maximum(np.maximum(M, X), Y)
            Mn = np.full_like(M, NEG)
            Mn[:, 1:] = prev[:, :-1] + score(rd[:, i:i + 1], win)[:, 1:]
            X = np.maximum(M - SC * (GAP_OPEN + GAP_EXTEND), X - SC * GAP_EXTEND)   # read base i inserted
            M = Mn
            Y = deletions(M)
        best = np.max(np.maximum(M, X), axis=1)              # free end on the haplotype
        r = np.mod(-best, SC)                                 # the winning start's tie-break rank
        off[s0 + need] = (r // 2) * np.where(r % 2 == 1, -1, 1)
    return off


def load_truth_reads(truth_bam, mason_vcf, reference, read_key, chunk=1_000_000):
    """{read_key(read): (contig, true_start, true_end, oH, place_start)} for every primary read of a
    mason -oa truth BAM, with the TRUE reference span taken from the read's haplotype position: the
    haplotype interval [oP, oP + read length) mapped back to the reference (Haplotype.ref_span).
    The truth BAM's POS/CIGAR is not used: POS is off for ~1% of reads in SNP-only simulations
    (mason's fake 1D..1I pairs at dense SNP clusters) and for ~30% at s02's indel density.

    true_start/true_end: the reference positions the read's bases came from -- use these to decide
    which sites a read truly covers. place_start: the start to compare an alignment's start with --
    where an aligner would lay the read down unclipped. It differs from true_start in two cases:
      - a sequencing indel among the read's leftmost bases (~0.02% of reads): an aligner absorbs it
        as a 1-3 bp shift plus mismatches. The start comes from aligning the read's leftmost
        ANCHOR_LEN bases to its haplotype (anchor_offsets, bwa mem scores); added at commit efa152b
        (2026-10-27). Tables made at 8762eab placed these reads at oP, so they counted as 'shifted'
        when bwa agreed with mason's POS.
      - a read that begins inside a simulated insertion (no reference position for its first
        bases): an aligner can place it without clipping only by laying the inserted bases over the
        reference bases before the insertion, so place_start = its first reference base minus the
        leading inserted bases (mason's own POS does the same).

    Assumption: the read's footprint on its haplotype is its length (query_length). mason's
    sequencing indels (~1.5% of reads) make the true footprint 1-2 bp longer or shorter, which
    moves true_end by that much; only SNPs within 2 bp of such a read's end are affected.

    Also returns (n_reads, n_start_differs, n_end_differs, n_anchor_shifted): how many reads'
    place_start / true_end differ from the truth BAM's POS / CIGAR end, and how many reads the
    leftmost-bases alignment placed off oP, for the caller's log."""
    import array
    import pysam
    keys, oh, op, ln = array.array("q"), array.array("b"), array.array("q"), array.array("q")
    pos_s, pos_e, contigs = array.array("q"), array.array("q"), set()
    left = bytearray()                                 # each read's leftmost ANCHOR_LEN bases, concatenated
    with pysam.AlignmentFile(truth_bam, "rb") as tf:
        for r in tf:
            if r.is_secondary or r.is_supplementary:
                continue
            keys.append(read_key(r)); oh.append(r.get_tag("oH")); op.append(r.get_tag("oP"))
            ln.append(r.query_length); pos_s.append(r.reference_start); pos_e.append(r.reference_end)
            contigs.add(r.reference_name)
            if r.query_length < ANCHOR_LEN:
                raise ValueError("%s: read %s is %d bp, shorter than ANCHOR_LEN %d"
                                 % (truth_bam, r.query_name, r.query_length, ANCHOR_LEN))
            left += r.query_sequence[:ANCHOR_LEN].upper().encode()   # SEQ is in reference orientation
    if len(contigs) != 1:
        raise ValueError("%s: reads on %d contigs (%s); one contig expected"
                         % (truth_bam, len(contigs), ", ".join(sorted(contigs))))
    contig = contigs.pop()
    with pysam.FastaFile(reference) as fa:
        haps = build(fa.fetch(contig), load_mason_vcf(mason_vcf, contig))
    oh = np.frombuffer(oh, dtype=np.int8)
    if not np.isin(oh, (1, 2)).all():
        raise ValueError("%s: oH tag outside {1, 2}" % truth_bam)
    op, ln = np.frombuffer(op, dtype=np.int64), np.frombuffer(ln, dtype=np.int64)
    left = np.frombuffer(bytes(left), dtype=np.uint8).reshape(len(op), ANCHOR_LEN)
    start, end, place, shift = (np.empty(len(op), dtype=np.int64) for _ in range(4))
    for k in (1, 2):                                   # each read's span on its own haplotype
        m = oh == k
        start[m], end[m], _ = haps[k - 1].ref_span(op[m], ln[m])
        # placement: where the read's leftmost bases align on the haplotype (oP + shift), then the
        # reference position there, minus any inserted bases it starts on
        shift[m] = anchor_offsets(haps[k - 1].seq, op[m], left[m])
        p_start, _, p_lead = haps[k - 1].ref_span(op[m] + shift[m], ln[m] - shift[m])
        place[m] = p_start - p_lead                    # start to compare alignments with
    diffs = (len(op), int(np.sum(place != np.frombuffer(pos_s, dtype=np.int64))),
             int(np.sum(end != np.frombuffer(pos_e, dtype=np.int64))), int(np.sum(shift != 0)))
    # build the lookup in chunks, so the temporary Python ints stay small
    truth, keys = {}, np.frombuffer(keys, dtype=np.int64)
    for c in range(0, len(keys), chunk):
        s = slice(c, c + chunk)
        for kk, a, b, h, p in zip(keys[s].tolist(), start[s].tolist(), end[s].tolist(),
                                  oh[s].tolist(), place[s].tolist()):
            truth[kk] = (contig, a, b, h, p)
    if len(truth) != len(keys):
        raise ValueError("%s: %d duplicate read keys" % (truth_bam, len(keys) - len(truth)))
    return truth, diffs
