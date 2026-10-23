#!/usr/bin/env python3
"""
neighborhood_assembly.py -- what changes about the "assemblable neighborhood" when a gatk SNP flips
not-called -> called, even though its pileup at p barely moves?

Hypothesis: gatk's reference-anchored reassembler fails while the LOCAL sequence around p is too
divergent from the (consensus) reference -- too many still-unrecovered true SNPs in the assembly
window mean too few shared exact k-mers, so the graph can't close. As iteration bakes neighboring SNPs
into the consensus, the local unrecovered-SNP density falls; once it is low enough the haplotype
assembles and the SNP is called. This is OFF p, so the at-p pileup decomposition can't see it.

We measure, for every sampled SNP and each consecutive iteration pair (k -> k+1), the number of OTHER
true hom-div SNPs within +/-W bp of p that are NOT yet recovered at that iteration (= the local
mismatch load the assembler faces). We compare this BEFORE (iter k) vs AFTER (iter k+1) for flips
(not-called -> called) against the stayed-not-called and stayed-called controls.

Recovered set per iteration = iter{k}/mm39_recovered.vcf.gz (consensus-vs-ref diff, cumulative).
Run from project root:  python analyses/s01-baseline-iterative/scripts/neighborhood_assembly.py
"""
import bisect
import csv
import gzip
import os

BASE = "data/s01-baseline-iterative"
PILE = "analyses/s01-baseline-iterative/results/vrf-pileup"
OUT = "analyses/s01-baseline-iterative/results/neighborhood"
os.makedirs(OUT, exist_ok=True)
SEED = "42"
W = 200                                   # +/- window (bp) ~ HaplotypeCaller active region + padding
CALLERS = ["gatk", "bcftools", "freebayes"]
DIVS = ["0.06", "0.08", "0.10"]


def snp_positions(vcf_path, hom_div_only):
    """Sorted list of 0-based SNP positions. hom_div_only -> keep GT 1|1 hom-div biallelic SNPs."""
    out = []
    with gzip.open(vcf_path, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            c = line.rstrip("\n").split("\t")
            if len(c[3]) != 1 or len(c[4]) != 1:
                continue
            if hom_div_only and c[9].split(":")[0] not in ("1|1", "1/1"):
                continue
            out.append(int(c[1]) - 1)
    out.sort()
    return out


def count_between(sorted_pos, lo, hi, exclude):
    """# positions in [lo, hi] excluding `exclude`, via binary search on a sorted list."""
    n = bisect.bisect_right(sorted_pos, hi) - bisect.bisect_left(sorted_pos, lo)
    # subtract the excluded position if it is present in [lo,hi]
    if lo <= exclude <= hi:
        i = bisect.bisect_left(sorted_pos, exclude)
        if i < len(sorted_pos) and sorted_pos[i] == exclude:
            n -= 1
    return n


rows = []          # per (caller, div, fate) aggregate accumulator
persnp = []        # per (caller, div, SNP, iteration) -- joins into the pileup/boxplot data
for D in DIVS:
    truth = snp_positions(f"{BASE}/sim/30X/19/{D}/0.005/{SEED}/mm39-truth.norm.vcf.gz", hom_div_only=True)
    for C in CALLERS:
        # recovered positions per iteration (cumulative consensus-vs-ref SNP diff)
        rec = {}
        for k in (1, 2, 3):
            rec[k] = snp_positions(f"{BASE}/iter/{C}/30X/19/{D}/0.005/{SEED}/iter{k}/mm39_recovered.vcf.gz",
                                   hom_div_only=False)
        # per-SNP recovered flags from the pileup CSV (same sampled SNPs the decomposition uses)
        traj = {}
        with open(f"{PILE}/vrf-pileup-{C}-{D}.csv") as fh:
            for r in csv.DictReader(fh):
                traj.setdefault(int(r["snp_pos"]) - 1, {})[int(r["iteration"])] = int(r["recovered"])
        for p, tr in traj.items():
            lo, hi = p - W, p + W
            tw = count_between(truth, lo, hi, p)                          # neighbor true hom-div SNPs
            # (a) per-iteration unrecovered neighbor count -> join into pileup data by (caller,div,snp,iter)
            for k in (1, 2, 3):
                if k in tr:
                    persnp.append((C, D, p + 1, k, tw, tw - count_between(rec[k], lo, hi, p)))
            # (b) per consecutive pair: fate + before/after, for the aggregate summary
            for k in (1, 2):
                if k not in tr or (k + 1) not in tr:
                    continue
                a, b = tr[k], tr[k + 1]
                fate = ("flip" if a == 0 and b == 1 else "stayed_notcalled" if a == 0 and b == 0
                        else "stayed_called" if a == 1 and b == 1 else "lost")
                unrec_before = tw - count_between(rec[k], lo, hi, p)      # not yet in consensus
                unrec_after = tw - count_between(rec[k + 1], lo, hi, p)
                rows.append((C, D, fate, tw, unrec_before, unrec_after))

# aggregate: mean neighbor-SNP counts by caller x divergence x fate
agg = {}
for C, D, fate, tw, ub, ua in rows:
    key = (C, D, fate)
    s = agg.setdefault(key, [0, 0.0, 0.0, 0.0])
    s[0] += 1; s[1] += tw; s[2] += ub; s[3] += ua

with open(f"{OUT}/neighborhood-per-snp.csv", "w") as out:
    out.write("caller,divergence,snp_pos,iteration,n_neighbor_truth,n_unrecovered\n")
    for r in persnp:
        out.write(",".join(str(x) for x in r) + "\n")

with open(f"{OUT}/neighborhood-summary.csv", "w") as out:
    out.write("caller,divergence,fate,window_bp,n,mean_neighbor_truth,"
              "mean_unrecovered_before,mean_unrecovered_after\n")
    for (C, D, fate), (n, tw, ub, ua) in sorted(agg.items()):
        out.write(f"{C},{D},{fate},{W},{n},{tw/n:.2f},{ub/n:.2f},{ua/n:.2f}\n")

print(f"window +/-{W}bp. Mean # still-UNRECOVERED neighbor hom-div SNPs around p, before -> after:")
print(f"{'caller':>9} {'D':>5} {'fate':>17} {'n':>6} {'unrec_before':>13} {'unrec_after':>12}")
for (C, D, fate), (n, tw, ub, ua) in sorted(agg.items()):
    if fate in ("flip", "stayed_notcalled", "stayed_called"):
        print(f"{C:>9} {D:>5} {fate:>17} {n:>6} {ub/n:>13.2f} {ua/n:>12.2f}")
