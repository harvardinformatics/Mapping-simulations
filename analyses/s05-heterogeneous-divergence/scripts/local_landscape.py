#!/usr/bin/env python3
"""
local_landscape.py -- small per-config tables that DESCRIBE each s05 divergence
landscape (for the notebook figures), derived from the ACTUAL truth sites.

Two outputs, both tidy over all field configs:
  1. <out_prefix>-landscape.csv : coarse binned local divergence ALONG the contig
     (per eligible base, N-gap bins flagged) -- shows what each field looks like
     spatially.
  2. <out_prefix>-localdist.csv : histogram of per-site local divergence at the
     read (150) and fragment (500) window scales -- shows how heterogeneity
     spreads the divergence a read experiences (vs a spike at the mean for the
     homogeneous control).

Divergent-site positions come from each config's mm39-truth-locald.tsv (the same
sites the pipeline scored). The reference is read once to mask N gaps so the
landscape is divergence per MAPPABLE base, not per genomic bp.
"""

import os
import sys
import argparse
import numpy as np
import pysam

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "lib"))
import divergence_field as dfield
import provenance                       # shared provenance block for the data dictionaries


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data-dir", required=True)
    p.add_argument("--reference", required=True, help="region FASTA (for length + N mask)")
    p.add_argument("--ref-str", default="mm39")
    p.add_argument("--region", default="19")
    p.add_argument("--cov", default="30")
    p.add_argument("--div", default="0.06")
    p.add_argument("--het", default="0.005")
    p.add_argument("--seed", default="42")
    p.add_argument("--configs", nargs="+", required=True)
    p.add_argument("--bin-size", type=int, default=50000, help="landscape bin (bp)")
    p.add_argument("--hist-step", type=float, default=0.01, help="local-D histogram bin width")
    p.add_argument("--hist-max", type=float, default=0.60, help="local-D histogram top edge")
    p.add_argument("--windows", nargs="+", default=["150", "500"])
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--dicts-only", metavar="ORIGINAL_COMMAND",
                   help="write ONLY the data dictionaries for existing outputs (same arguments as the "
                        "run that made them), recording ORIGINAL_COMMAND as provenance")
    p.add_argument("--dict-note", action="append", default=[],
                   help="extra provenance note for the dictionaries (repeatable)")
    return p.parse_args()


def write_dictionaries(args):
    """One data dictionary per output table. With --dicts-only they document existing tables
    (backfill) under the ORIGINAL command."""
    def hdr(table):
        if args.dicts_only:
            return provenance.header(table, command=args.dicts_only, backfill=True, notes=args.dict_note)
        return provenance.header(table, notes=args.dict_note)

    sites = ("Sites: each config's FIXED-divergence (1|1) truth SNPs, with local divergence from\n"
             "its mm39-truth-locald.tsv (annotate_local_divergence.py). Configs: %s.\n\n"
             % ", ".join(args.configs))
    land = args.out_prefix + "-landscape.csv"
    with open(land + ".dict.txt", "w") as d:
        d.write(hdr(land) + sites)
        d.write("Local divergence along chr%s in %d-bp bins; one row per (config, bin).\n\n"
                % (args.region, args.bin_size))
        d.write("Fields:\n")
        d.write("  field/corr_len/sigma : divergence-field config tag L{corr_len}-S{sigma}[-C{cap}]\n")
        d.write("                         (mean D fixed at %s)\n" % args.div)
        d.write("  bin_center           : bp center of the bin\n")
        d.write("  local_divergence     : divergent sites per ELIGIBLE (non-N) base in the bin\n")
        d.write("  valid                : 1 if the bin is >=50% non-N (else density is meaningless)\n")
    dist = args.out_prefix + "-localdist.csv"
    with open(dist + ".dict.txt", "w") as d:
        d.write(hdr(dist) + sites)
        d.write("Histogram of per-site local divergence; one row per (config, window, bin).\n\n")
        d.write("Fields:\n")
        d.write("  field/corr_len/sigma : divergence-field config (as above)\n")
        d.write("  window     : local-D window half-size in bp (150 = read scale, 500 = fragment)\n")
        d.write("  locald_mid : histogram bin midpoint; bins are %.2g wide from 0 to %.2g, plus a\n"
                % (args.hist_step, args.hist_max))
        d.write("               '>%.2g' overflow row\n" % args.hist_max)
        d.write("  n_sites    : truth divergent SNP sites with local divergence in that bin\n\n")
        d.write("The homogeneous control (sigma=0) concentrates near the mean; larger sigma\n")
        d.write("spreads sites across a wide range of local divergence (the heterogeneity).\n")


def parse_field(tag):
    # 'L20000-S1.0' or 'L20000-S1.0-C0.20' -> (corr_len, sigma); cap (if any) ignored here
    parts = tag.split("-")
    return parts[0][1:], parts[1][1:]


def load_locald(path, windows):
    pos = []
    cols = {w: [] for w in windows}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ipos = header.index("POS")
        idx = {w: header.index("local_d_%s" % w) for w in windows}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            pos.append(int(f[ipos]))
            for w in windows:
                cols[w].append(float(f[idx[w]]))
    return (np.array(pos, dtype=np.int64) - 1,  # -> 0-based for binning
            {w: np.array(cols[w], dtype=np.float64) for w in windows})


def main():
    args = parse_args()
    # --dicts-only: document existing outputs without recomputing them
    if args.dicts_only:
        write_dictionaries(args)
        return

    # reference length + eligible (non-N) positions, once
    fa = pysam.FastaFile(args.reference)
    name = fa.references[0]
    seq = fa.fetch(name).upper()
    arr = np.frombuffer(seq.encode("ascii"), dtype="S1")
    valid = (arr == b"A") | (arr == b"C") | (arr == b"G") | (arr == b"T")
    eligible = np.where(valid)[0]
    length = len(seq)

    # histogram edges for per-site local divergence
    hedges = np.arange(0.0, args.hist_max + args.hist_step, args.hist_step)

    land_rows = []
    dist_rows = []
    for cfg in args.configs:
        corr, sigma = parse_field(cfg)
        tsv = os.path.join(args.data_dir, "sim", args.cov + "X", args.region,
                           args.div, args.het, cfg, args.seed,
                           args.ref_str + "-truth-locald.tsv")
        if not os.path.isfile(tsv):
            sys.exit("ERROR: missing %s" % tsv)
        pos0, local = load_locald(tsv, args.windows)

        # landscape: local divergence per eligible base in bins along the contig
        centers, rate, bvalid = dfield.binned_rate(pos0, length, args.bin_size, eligible)
        for c, r, v in zip(centers, rate, bvalid):
            land_rows.append((cfg, corr, sigma, int(c), r, int(v)))

        # per-site local-D histogram at each window scale
        for w in args.windows:
            counts, _ = np.histogram(local[w], bins=hedges)
            for k in range(len(hedges) - 1):
                mid = (hedges[k] + hedges[k + 1]) / 2.0
                dist_rows.append((cfg, corr, sigma, w, "%.4g" % mid, int(counts[k])))
            # overflow above hist_max
            over = int((local[w] > args.hist_max).sum())
            dist_rows.append((cfg, corr, sigma, w, ">%.2g" % args.hist_max, over))

    with open(args.out_prefix + "-landscape.csv", "w") as out:
        out.write("field,corr_len,sigma,bin_center,local_divergence,valid\n")
        for r in land_rows:
            out.write("%s,%s,%s,%d,%.8g,%d\n" % r)

    with open(args.out_prefix + "-localdist.csv", "w") as out:
        out.write("field,corr_len,sigma,window,locald_mid,n_sites\n")
        for r in dist_rows:
            out.write("%s,%s,%s,%s,%s,%d\n" % r)

    write_dictionaries(args)

    sys.stderr.write("# landscape rows=%d  localdist rows=%d\n" % (len(land_rows), len(dist_rows)))


if __name__ == "__main__":
    main()
