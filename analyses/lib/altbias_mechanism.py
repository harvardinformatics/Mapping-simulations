#!/usr/bin/env python
"""Direct test of the s03 reference-bias mechanism: at each true het site, split the
overlapping primary reads by the base they carry (ref vs alt) and compare their counts,
MAPQ, and soft-clip rate. NO MAPQ filter, so low-MAPQ reads are visible.

The key output is the alt/ref *count* ratio among reads with an aligned base at the het
position: < 1 means alt-bearing reads are missing from the pileup (lost / mismapped /
clipped over the variant) -- i.e. reference bias is alt-read EXCLUSION. The MAPQ and clip
columns show whether the alt reads that *survive* are also degraded (they are not: the bias
is in the reads that vanish, not the ones we see).

Usage:
  altbias_mechanism.py BAM TRUTH.vcf[.gz] REGION [--stride N] [--label "div,iter,caller"] [--csv]
Het positions are reference coordinates; run on identity-coordinate BAMs (s01, or s01/s02
iter1). For length-shifted consensus BAMs (s02 iter>1) the het positions would need lifting.
"""
import argparse, pysam


def load_het(path):
    """[(pos0, ref, alt)] for heterozygous SNP truth sites (GT 0/1 or 1/0)."""
    out = []
    op = pysam.BGZFile(path) if path.endswith(".gz") else open(path, "rb")
    for raw in op:
        line = raw.decode() if isinstance(raw, bytes) else raw
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        ref, alt, gt = c[3], c[4], c[9].split(":")[0]
        if len(ref) == 1 and len(alt) == 1 and gt in ("0|1", "1|0", "0/1", "1/0"):
            out.append((int(c[1]) - 1, ref.upper(), alt.upper()))
    return out


def measure(bam_path, truth_path, region, stride):
    het = load_het(truth_path)[::stride]
    bam = pysam.AlignmentFile(bam_path, "rb")
    n = {"ref": 0, "alt": 0, "other": 0}
    mapq_sum = {"ref": 0, "alt": 0}
    lowmapq = {"ref": 0, "alt": 0}   # MAPQ < 20 (would be filtered by the RAF tool)
    clip = {"ref": 0, "alt": 0}      # read has any soft-clip in its CIGAR
    for pos0, refb, altb in het:
        for r in bam.fetch(region, pos0, pos0 + 1):
            if r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_unmapped:
                continue
            qpos = None
            for (qp, rp) in r.get_aligned_pairs(matches_only=True):
                if rp == pos0:
                    qpos = qp
                    break
            if qpos is None:
                continue  # het position falls in an indel/soft-clip of this read
            base = r.query_sequence[qpos].upper()
            cls = "ref" if base == refb else ("alt" if base == altb else "other")
            n[cls] += 1
            if cls in ("ref", "alt"):
                mapq_sum[cls] += r.mapping_quality
                if r.mapping_quality < 20:
                    lowmapq[cls] += 1
                if any(op == 4 for op, _ in r.cigartuples):
                    clip[cls] += 1
    return len(het), n, mapq_sum, lowmapq, clip


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("bam"); ap.add_argument("truth"); ap.add_argument("region")
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--label", default="", help="comma fields prepended to the --csv row")
    ap.add_argument("--csv", action="store_true", help="emit one CSV data row (no header)")
    a = ap.parse_args()

    sites, n, mq, low, clip = measure(a.bam, a.truth, a.region, a.stride)
    ratio = (n["alt"] / n["ref"]) if n["ref"] else float("nan")
    mq_ref = (mq["ref"] / n["ref"]) if n["ref"] else float("nan")
    mq_alt = (mq["alt"] / n["alt"]) if n["alt"] else float("nan")
    pct = lambda x, b: (100.0 * x / b) if b else float("nan")

    if a.csv:
        fields = ([f for f in a.label.split(",") if f != ""] +
                  [sites, n["ref"], n["alt"], "%.4f" % ratio,
                   "%.3f" % mq_ref, "%.3f" % mq_alt,
                   "%.3f" % pct(low["ref"], n["ref"]), "%.3f" % pct(low["alt"], n["alt"]),
                   "%.3f" % pct(clip["ref"], n["ref"]), "%.3f" % pct(clip["alt"], n["alt"])])
        print(",".join(str(x) for x in fields))
    else:
        print("het sites sampled: %d  (stride %d)" % (sites, a.stride))
        print("aligned reads at het positions:  ref=%d  alt=%d  other=%d" % (n["ref"], n["alt"], n["other"]))
        print("  alt/ref count ratio           : %.3f   (1.0 = no loss; <1 = alt reads missing)" % ratio)
        print("mean MAPQ                        : ref=%.2f   alt=%.2f   (delta=%.2f)" % (mq_ref, mq_alt, mq_ref - mq_alt))
        print("%% reads with MAPQ<20 (filtered)  : ref=%.2f%%  alt=%.2f%%" % (pct(low["ref"], n["ref"]), pct(low["alt"], n["alt"])))
        print("%% reads soft-clipped             : ref=%.2f%%  alt=%.2f%%" % (pct(clip["ref"], n["ref"]), pct(clip["alt"], n["alt"])))


if __name__ == "__main__":
    main()
