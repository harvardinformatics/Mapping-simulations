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
    truth, diffs = th.load_truth_reads(bam, vcf, fasta, read_key)   # each truth read's TRUE span
Every check fails loudly (ValueError); nothing is skipped silently.
"""
import gzip

import numpy as np

GT_ALLELES = {"1|1": (1, 1), "0|1": (0, 1), "1|0": (1, 0)}   # phased genotypes the simulations write


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


def load_truth_reads(truth_bam, mason_vcf, reference, read_key, chunk=1_000_000):
    """{read_key(read): (contig, true_start, true_end, oH, place_start)} for every primary read of a
    mason -oa truth BAM, with the TRUE reference span taken from the read's haplotype position: the
    haplotype interval [oP, oP + read length) mapped back to the reference (Haplotype.ref_span).
    The truth BAM's POS/CIGAR is not used: POS is off for ~1% of reads in SNP-only simulations
    (mason's fake 1D..1I pairs at dense SNP clusters) and for ~30% at s02's indel density.

    true_start/true_end: the reference positions the read's bases came from -- use these to decide
    which sites a read truly covers. place_start: the start to compare an alignment's start with.
    It equals true_start except for a read that begins inside an insertion (no reference position
    for its first bases): an aligner can place such a read without clipping only by laying the
    inserted bases over the reference bases before the insertion, so place_start = true_start minus
    the number of leading inserted bases (mason's own POS does the same).

    Assumption: the read's footprint on its haplotype is its length (query_length). mason's
    sequencing indels (~1.5% of reads) make the true footprint 1-2 bp longer or shorter, which
    moves true_end by that much; only SNPs within 2 bp of such a read's end are affected.

    Also returns (n_reads, n_start_differs, n_end_differs): how many reads' place_start / true_end
    differ from the truth BAM's POS / CIGAR end, for the caller's log."""
    import array
    import pysam
    keys, oh, op, ln = array.array("q"), array.array("b"), array.array("q"), array.array("q")
    pos_s, pos_e, contigs = array.array("q"), array.array("q"), set()
    with pysam.AlignmentFile(truth_bam, "rb") as tf:
        for r in tf:
            if r.is_secondary or r.is_supplementary:
                continue
            keys.append(read_key(r)); oh.append(r.get_tag("oH")); op.append(r.get_tag("oP"))
            ln.append(r.query_length); pos_s.append(r.reference_start); pos_e.append(r.reference_end)
            contigs.add(r.reference_name)
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
    start, end, lead = (np.empty(len(op), dtype=np.int64) for _ in range(3))
    for k in (1, 2):                                   # each read's span on its own haplotype
        m = oh == k
        start[m], end[m], lead[m] = haps[k - 1].ref_span(op[m], ln[m])
    place = start - lead                               # start to compare alignments with
    diffs = (len(op), int(np.sum(place != np.frombuffer(pos_s, dtype=np.int64))),
             int(np.sum(end != np.frombuffer(pos_e, dtype=np.int64))))
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
