#!/usr/bin/env python3
"""
neighborhood_assembly.py -- what changes about the "assemblable neighborhood" when a variant flips
not-called -> called, even though its pileup at p barely moves?

Hypothesis: gatk's reference-anchored reassembler fails while the LOCAL sequence around p is too
divergent from the (consensus) reference -- too many still-unrecovered true SNPs in the assembly
window mean too few shared exact k-mers, so the graph can't close. As iteration bakes neighboring SNPs
into the consensus, the local unrecovered-SNP density falls; once low enough the haplotype assembles
and the variant is called. This is OFF p, so the at-p pileup decomposition can't see it.

For every sampled query variant (from the --pile CSVs) and each consecutive iteration pair (k -> k+1),
we count OTHER true hom-div SNPs within +/-W bp of p that are NOT yet recovered at that iteration (the
local mismatch load), and compare BEFORE (iter k) vs AFTER (iter k+1) across fates (flip,
stayed-not-called, stayed-called). Path templates make it reusable across sections:
  s01: --pile s01/.../vrf-pileup-{C}-{D}.csv  --recovered s01/.../iter{k}/mm39_recovered.vcf.gz
  s02: recovered name = mm39_recovered.acc.vcf.gz (accumulated, reference coords)
  s03: query = het sites (its own --pile), but neighbor truth + recovered are the s01 hom-div sets
All coordinates are reference coordinates (truth, recovered, and the window), so no chain is needed.

Templates use {C} (caller), {D} (divergence), {k} (iteration). Run from project root.
"""
import argparse
import bisect
import csv
import gzip
import os

import provenance                       # shared provenance block for the data dictionaries


def snp_positions(vcf_path, hom_div_only):
    """Sorted list of 0-based SNP positions. hom_div_only -> keep GT 1|1 hom-div biallelic SNPs."""
    return variant_positions(vcf_path, hom_div_only)[0]


def variant_positions(vcf_path, hom_div_only):
    """Single pass over a VCF returning (snp_positions, indel_positions), both sorted 0-based.
    A record is a SNP when ref and alt are each length 1, otherwise an indel. hom_div_only keeps
    only GT 1|1 records (the fixed sample-vs-reference variants), applied to SNPs AND indels alike.
    Reading both classes in one pass avoids re-parsing the large accumulated recovered VCFs twice."""
    snps, indels = [], []
    op = gzip.open(vcf_path, "rt") if vcf_path.endswith(".gz") else open(vcf_path)
    for line in op:
        if line[0] == "#":
            continue
        c = line.rstrip("\n").split("\t")
        if hom_div_only and c[9].split(":")[0] not in ("1|1", "1/1"):
            continue
        pos0 = int(c[1]) - 1
        if len(c[3]) == 1 and len(c[4]) == 1:
            snps.append(pos0)
        else:
            indels.append(pos0)
    snps.sort()
    indels.sort()
    return snps, indels


def count_between(sorted_pos, lo, hi, exclude):
    """# positions in [lo, hi] excluding `exclude`, via binary search on a sorted list."""
    n = bisect.bisect_right(sorted_pos, hi) - bisect.bisect_left(sorted_pos, lo)
    if lo <= exclude <= hi:
        i = bisect.bisect_left(sorted_pos, exclude)
        if i < len(sorted_pos) and sorted_pos[i] == exclude:
            n -= 1
    return n


def _in_sorted(sorted_pos, x):
    i = bisect.bisect_left(sorted_pos, x)
    return i < len(sorted_pos) and sorted_pos[i] == x


def unrec_in_window(truth_sorted, rec_sorted, lo, hi, exclude):
    """# of TRUTH SNPs in [lo, hi] (excluding `exclude`) that are NOT members of the recovered set.
    Counts by set membership -- NOT truth_count - recovered_count -- so false-positive recovered SNPs
    at non-truth positions (freebayes over-calls heavily) cannot drive the count negative."""
    i = bisect.bisect_left(truth_sorted, lo)
    j = bisect.bisect_right(truth_sorted, hi)
    cnt = 0
    for x in truth_sorted[i:j]:
        if x == exclude:
            continue
        if not _in_sorted(rec_sorted, x):
            cnt += 1
    return cnt


def write_dictionaries(a):
    """Data dictionaries for neighborhood-per-snp.csv and neighborhood-summary.csv. With --dicts-only
    they document existing tables (backfill), recording the ORIGINAL command and any --dict-note."""
    common = ("Query sites: the SNPs sampled in the --pile per-SNP pileup CSVs (template %s);\n"
              "             their per-iteration 'recovered' flags come from that file.\n"
              "Neighbours : OTHER true hom-div (GT 1|1) variants within +/-%d bp of the query SNP p,\n"
              "             from the truth VCF (template %s). A record is a SNP when REF and ALT are\n"
              "             both 1 bp, otherwise an indel.\n"
              "Recovered  : per-iteration recovered VCF (template %s), i.e. the variants\n"
              "             already baked into the consensus at that iteration.\n"
              "Unrecovered: a neighbour truth variant whose position is NOT in the recovered set\n"
              "             (set membership, so false-positive recovered calls cannot lower the count).\n"
              "All coordinates are reference coordinates.\n\n"
              % (a.pile, a.window, a.truth, a.recovered))

    def hdr(table):
        if a.dicts_only:
            return provenance.header(table, command=a.dicts_only, backfill=True, notes=a.dict_note)
        return provenance.header(table, notes=a.dict_note)

    per_snp_text = ("One row per (caller, divergence, query SNP, iteration).\n\n"
                    "Fields:\n"
                    "  caller                 : variant caller (gatk | bcftools | freebayes)\n"
                    "  divergence             : simulated hom-div divergence D\n"
                    "  snp_pos                : 1-based reference position of the query SNP p\n"
                    "  iteration              : mapping iteration (1 = single-shot)\n"
                    "  n_neighbor_truth       : neighbour true hom-div SNPs in [p-W, p+W], excluding p\n"
                    "  n_unrecovered          : of those, NOT yet recovered at this iteration\n"
                    "  n_neighbor_truth_indel : neighbour true hom-div indels in the window\n"
                    "  n_unrecovered_indel    : of those, NOT yet recovered at this iteration. Indel\n"
                    "                           recovery is by POSITION match of the normalized truth vs\n"
                    "                           recovered indels, so it slightly understates recovery.\n")
    summary_text = ("For each query SNP and each consecutive iteration pair (k -> k+1) the SNP's fate over\n"
                    "the pair is classified, and the neighbour counts are averaged per (caller, divergence,\n"
                    "fate) group. Neighbour counts here are SNPs only.\n\n"
                    "Fields:\n"
                    "  caller                  : variant caller (gatk | bcftools | freebayes)\n"
                    "  divergence              : simulated hom-div divergence D\n"
                    "  fate                    : flip (not called at k -> called at k+1) | stayed_notcalled |\n"
                    "                            stayed_called | lost (called at k -> not called at k+1)\n"
                    "  window_bp               : half-window W (bp)\n"
                    "  n                       : number of (SNP, iteration pair) instances in the group\n"
                    "  mean_neighbor_truth     : mean neighbour true hom-div SNPs in the window\n"
                    "  mean_unrecovered_before : mean neighbours NOT yet recovered at iteration k\n"
                    "  mean_unrecovered_after  : mean neighbours NOT yet recovered at iteration k+1\n")
    tables = {"per-snp": ("neighborhood-per-snp.csv", per_snp_text),
              "summary": ("neighborhood-summary.csv", summary_text)}
    for name, (fname, text) in tables.items():
        if a.dicts_only_table not in ("both", name):
            continue                     # --dicts-only-table: leave the other dictionary alone
        table = os.path.join(a.out, fname)
        with open(table + ".dict.txt", "w") as d:
            d.write(hdr(table) + common + text)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", required=True, help="template for the sim truth VCF (hom-div NEIGHBOR SNPs); {D}")
    ap.add_argument("--recovered", required=True, help="template for the per-iter recovered VCF; {C},{D},{k}")
    ap.add_argument("--pile", required=True, help="template for the per-SNP pileup CSV (QUERY sites); {C},{D}")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--callers", default="gatk,bcftools,freebayes")
    ap.add_argument("--divs", default="0.06,0.08,0.10")
    ap.add_argument("--window", type=int, default=200, help="+/- window (bp) ~ assembly region")
    ap.add_argument("--dicts-only", metavar="ORIGINAL_COMMAND",
                    help="write ONLY the data dictionaries for existing outputs in --out (same other "
                         "arguments as the run that made them), recording ORIGINAL_COMMAND as provenance")
    ap.add_argument("--dict-note", action="append", default=[],
                    help="extra provenance note for the dictionaries (repeatable)")
    ap.add_argument("--dicts-only-table", choices=["both", "per-snp", "summary"], default="both",
                    help="with --dicts-only: which table's dictionary to write (default both)")
    a = ap.parse_args()
    if a.dicts_only_table != "both" and not a.dicts_only:
        ap.error("--dicts-only-table needs --dicts-only")

    # --dicts-only: document existing outputs without recomputing them
    if a.dicts_only:
        write_dictionaries(a)
        return
    os.makedirs(a.out, exist_ok=True)
    W = a.window
    CALLERS = a.callers.split(",")
    DIVS = a.divs.split(",")

    rows = []          # per (caller, div, fate) aggregate accumulator
    persnp = []        # per (caller, div, SNP, iteration) -- joins into the pileup/boxplot data
    for D in DIVS:
        # neighbour truth = the fixed (1|1) divergent SNPs and indels around each query SNP
        truth, truth_indel = variant_positions(a.truth.format(D=D), hom_div_only=True)
        for C in CALLERS:
            rec, rec_indel = {}, {}                       # recovered SNP / indel positions per iter
            for k in (1, 2, 3):
                rec[k], rec_indel[k] = variant_positions(
                    a.recovered.format(C=C, D=D, k=k), hom_div_only=False)
            traj = {}
            with open(a.pile.format(C=C, D=D)) as fh:
                for r in csv.DictReader(fh):
                    traj.setdefault(int(r["snp_pos"]) - 1, {})[int(r["iteration"])] = int(r["recovered"])
            for p, tr in traj.items():
                lo, hi = p - W, p + W
                tw = count_between(truth, lo, hi, p)                # neighbor true hom-div SNPs
                tw_i = count_between(truth_indel, lo, hi, p)        # neighbor true indels (p is a SNP)
                for k in (1, 2, 3):
                    if k in tr:
                        persnp.append((C, D, p + 1, k, tw,
                                       unrec_in_window(truth, rec[k], lo, hi, p),
                                       tw_i,
                                       unrec_in_window(truth_indel, rec_indel[k], lo, hi, p)))
                for k in (1, 2):
                    if k not in tr or (k + 1) not in tr:
                        continue
                    aa, bb = tr[k], tr[k + 1]
                    fate = ("flip" if aa == 0 and bb == 1 else "stayed_notcalled" if aa == 0 and bb == 0
                            else "stayed_called" if aa == 1 and bb == 1 else "lost")
                    ub = unrec_in_window(truth, rec[k], lo, hi, p)
                    ua = unrec_in_window(truth, rec[k + 1], lo, hi, p)
                    rows.append((C, D, fate, tw, ub, ua))

    with open(f"{a.out}/neighborhood-per-snp.csv", "w") as out:
        # n_neighbor_truth / n_unrecovered are SNPs; the *_indel columns mirror them for indels,
        # so the per-SNP file supports the neighbour SNP-vs-indel unrecovered-load comparison.
        out.write("caller,divergence,snp_pos,iteration,n_neighbor_truth,n_unrecovered,"
                  "n_neighbor_truth_indel,n_unrecovered_indel\n")
        for r in persnp:
            out.write(",".join(str(x) for x in r) + "\n")

    agg = {}
    for C, D, fate, tw, ub, ua in rows:
        s = agg.setdefault((C, D, fate), [0, 0.0, 0.0, 0.0])
        s[0] += 1; s[1] += tw; s[2] += ub; s[3] += ua
    with open(f"{a.out}/neighborhood-summary.csv", "w") as out:
        out.write("caller,divergence,fate,window_bp,n,mean_neighbor_truth,"
                  "mean_unrecovered_before,mean_unrecovered_after\n")
        for (C, D, fate), (n, tw, ub, ua) in sorted(agg.items()):
            out.write(f"{C},{D},{fate},{W},{n},{tw/n:.2f},{ub/n:.2f},{ua/n:.2f}\n")
    write_dictionaries(a)

    print(f"window +/-{W}bp. Mean # still-UNRECOVERED neighbor hom-div SNPs around p, before -> after:")
    print(f"{'caller':>9} {'D':>5} {'fate':>17} {'n':>6} {'unrec_before':>13} {'unrec_after':>12}")
    for (C, D, fate), (n, tw, ub, ua) in sorted(agg.items()):
        if fate in ("flip", "stayed_notcalled", "stayed_called"):
            print(f"{C:>9} {D:>5} {fate:>17} {n:>6} {ub/n:>13.2f} {ua/n:>12.2f}")


if __name__ == "__main__":
    main()
