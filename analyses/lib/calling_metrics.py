#!/usr/bin/env python3
"""
calling_metrics.py -- per-SNP CALLING-side covariates: local properties of the region
that determine whether gatk's local reassembler + genotyper can act, GIVEN the reads are
present. Complements the existing "unrecovered-neighbour SNP" metric with the other axes
of calling-limitation (see the discussion): assembly-graph complexity (variant load),
graph resolvability (low-complexity sequence), and genotyping complexity (heterozygosity).

Calling-limited is essentially a gatk-assembler phenomenon: bcftools/freebayes call what
is in the pileup. So these metrics should separate gatk's calling-limited sites (alt present
but uncalled) from its recovered sites -- while barely moving for the pileup callers.

Metrics, per divergent (1|1) truth SNP position p, over a +/- window (default 100 bp):
  * het_density_W      : heterozygous truth SNPs / bp near p         (genotyping complexity)
  * indel_density_W    : indel truth sites / bp near p               (assembly complexity;
                         indels stress the assembler far more than SNPs)
  * dist_to_indel      : bp to nearest indel truth site              (assembly complexity)
  * snp_density_W      : divergent truth SNPs / bp near p            (assembly complexity)
  * homopolymer_frac_W : fraction of the window in homopolymer runs >= 5 bp (low-complexity
                         sequence that breaks the de Bruijn assembly graph)

These are STATIC properties of the truth + reference (they do not change across iterations),
so downstream they are compared BY FATE (recovered / calling-limited / mapping-limited),
not before/after. All computed vectorized over every divergent SNP.
"""

import sys
import argparse
import gzip
import numpy as np
import pysam
import provenance                       # shared provenance block for the data dictionaries


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--truth-vcf", required=True, help="mason-format truth VCF (div_snp/div_indel/het_snp IDs)")
    p.add_argument("--reference", required=True, help="region FASTA (single contig)")
    p.add_argument("--divergence", required=True, help="divergence label to stamp in the output")
    p.add_argument("--window", type=int, default=100, help="half-window (bp); metrics use +/- this")
    p.add_argument("--hp-min", type=int, default=5, help="min homopolymer run length to flag low-complexity")
    p.add_argument("--out", required=True)
    return p.parse_args()


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def load_positions(vcf):
    """Sorted int64 arrays of 1-based POS for div-SNP, het-SNP, indel truth sites."""
    div, het, ind = [], [], []
    with opener(vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            vid, pos = f[2], int(f[1])
            if vid.startswith("div_snp"):
                div.append(pos)
            elif vid.startswith("het_snp"):
                het.append(pos)
            elif vid.startswith("div_indel") or vid.startswith("het_indel"):
                ind.append(pos)
    return (np.sort(np.array(div, dtype=np.int64)),
            np.sort(np.array(het, dtype=np.int64)),
            np.sort(np.array(ind, dtype=np.int64)))


def density_in_window(query, source_sorted, W):
    """Count of source positions within +/- W of each query, / (2W)."""
    lo = np.searchsorted(source_sorted, query - W, side="left")
    hi = np.searchsorted(source_sorted, query + W, side="right")
    return (hi - lo) / (2.0 * W)


def nearest_distance(query, source_sorted):
    if source_sorted.size == 0:
        return np.full(query.shape, np.inf)
    n = source_sorted.size
    idx = np.searchsorted(source_sorted, query)
    ld = np.where(idx > 0, query - source_sorted[np.clip(idx - 1, 0, n - 1)], np.inf)
    rd = np.where(idx < n, source_sorted[np.clip(idx, 0, n - 1)] - query, np.inf)
    return np.minimum(ld, rd)


def homopolymer_prefix(seq_bytes, hp_min):
    """Prefix sum of a per-position flag = 1 if the position is inside a homopolymer run
    of length >= hp_min. Lets any window's homopolymer content be read in O(1)."""
    n = seq_bytes.size
    same = seq_bytes[1:] == seq_bytes[:-1]              # same[i] : base i+1 == base i
    # run id: increments where the base changes
    run_id = np.concatenate([[0], np.cumsum(~same)])
    # length of each run, mapped back to positions
    _, inv, counts = np.unique(run_id, return_inverse=True, return_counts=True)
    run_len = counts[inv]
    flag = (run_len >= hp_min).astype(np.int64)
    return np.concatenate([[0], np.cumsum(flag)])        # prefix[i] = sum of flag[0:i]


def main():
    a = parse_args()
    W = a.window

    fa = pysam.FastaFile(a.reference)
    name = fa.references[0]
    seq = fa.fetch(name).upper()
    sb = np.frombuffer(seq.encode("ascii"), dtype="S1")
    hp_prefix = homopolymer_prefix(sb, a.hp_min)
    L = sb.size

    div, het, ind = load_positions(a.truth_vcf)
    if div.size == 0:
        sys.exit("ERROR: no div_snp positions in %s" % a.truth_vcf)

    het_d = density_in_window(div, het, W)
    ind_d = density_in_window(div, ind, W)
    snp_d = density_in_window(div, div, W)               # includes self
    d2ind = nearest_distance(div, ind)
    # homopolymer fraction in +/- W around each SNP (0-based window, clipped to contig)
    lo = np.clip(div - 1 - W, 0, L)
    hi = np.clip(div - 1 + W, 0, L)
    hp_frac = (hp_prefix[hi] - hp_prefix[lo]) / np.maximum(hi - lo, 1)

    with open(a.out, "w") as out:
        out.write("divergence,snp_pos,het_density,indel_density,dist_to_indel,"
                  "snp_density,homopolymer_frac\n")
        big = 10 ** 9
        for i in range(div.size):
            di = d2ind[i]
            out.write("%s,%d,%.6g,%.6g,%d,%.6g,%.6g\n" %
                      (a.divergence, div[i], het_d[i], ind_d[i],
                       (int(di) if np.isfinite(di) else big), snp_d[i], hp_frac[i]))

    with open(a.out + ".dict.txt", "w") as d:
        d.write(provenance.header(a.out))
        d.write("Per divergent (1|1) SNP: local CALLING-side covariates (window +/- %d bp).\n\n" % W)
        d.write("Fields:\n")
        d.write("  divergence      : condition label\n")
        d.write("  snp_pos         : 1-based reference position of the divergent SNP\n")
        d.write("  het_density     : heterozygous truth SNPs per bp near p (genotyping complexity)\n")
        d.write("  indel_density   : indel truth sites per bp near p (assembly complexity)\n")
        d.write("  dist_to_indel   : bp to nearest indel truth site (%d if none within the contig)\n" % (10**9))
        d.write("  snp_density     : divergent truth SNPs per bp near p, incl. self (assembly complexity)\n")
        d.write("  homopolymer_frac: fraction of the window in homopolymer runs >= %d bp\n" % a.hp_min)
        d.write("                    (low-complexity sequence that breaks de Bruijn assembly)\n\n")
        d.write("Static per site (do not change across iterations) -> compare BY FATE downstream.\n")

    sys.stderr.write("# %s: %d div SNPs; het=%d indel=%d; mean indel_density=%.4g mean hp_frac=%.4g\n"
                     % (a.divergence, div.size, het.size, ind.size, ind_d.mean(), hp_frac.mean()))


if __name__ == "__main__":
    main()
