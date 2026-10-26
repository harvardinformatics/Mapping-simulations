#!/usr/bin/env python3
"""
make_divergence_field.py -- generate and VALIDATE an s05 continuous autocorrelated
divergence landscape, independent of read simulation.

This is the validation-gate tool: it exercises the exact placement path that
make_truth_vcf.py uses (via the shared divergence_field module), then reports
whether the realized landscape has the properties we asked for:

  * genome-wide mean divergence pinned to the target D (by construction),
  * realized correlation length ~ the requested corr_len,
  * amplitude (spread of local divergence) that grows with sigma, and
  * sigma = 0 reproducing a flat/uniform landscape.

It does NOT write a truth VCF or place any reads. It writes:
  * <out_prefix>-landscape.tsv : coarse binned local-divergence track (for plots),
  * <out_prefix>-localdist.tsv : distribution of per-site local divergence at each
                                 read/fragment window scale,
  * <out_prefix>-summary.txt   : realized statistics + the parameters used, and
  * <out_prefix>-summary.txt.dict.txt : data dictionary for the above.

Only A/C/G/T bases are eligible (N bases excluded), matching make_truth_vcf.py, so
the realized landscape reflects the same mask the real truth VCF will see.
"""

import sys
import argparse
import numpy as np
import pysam

sys.path.insert(0, __file__.rsplit("/", 1)[0])  # allow import of sibling module
import divergence_field as dfield
import provenance                       # shared provenance block for the data dictionaries


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--reference", required=True,
                   help="FASTA with a SINGLE sequence (same one make_truth_vcf uses).")
    p.add_argument("--divergence", type=float, required=True,
                   help="Target mean SNP divergence D (genome-wide fraction).")
    p.add_argument("--corr-len", type=float, required=True,
                   help="Correlation length of the field in bp (knob 1: scale).")
    p.add_argument("--sigma", type=float, required=True,
                   help="Log-amplitude of the field (knob 2: strength; 0 = uniform).")
    p.add_argument("--cap", type=float, default=None,
                   help="Biological ceiling on local divergence (e.g. 0.20). Bounds the "
                        "field's hot regions; realized ceiling reported below to tune it.")
    p.add_argument("--field-seed", type=int, required=True,
                   help="Seed for the latent field (independent of the site-sampling seed).")
    p.add_argument("--sample-seed", type=int, required=True,
                   help="Seed for drawing site positions given the weights.")
    p.add_argument("--bin-size", type=int, default=10000,
                   help="Bin size (bp) for the coarse landscape track (default 10 kb).")
    p.add_argument("--windows", default="150,500",
                   help="Comma-separated half-window sizes (bp) for per-site local "
                        "divergence (default read scale 150, fragment scale 500).")
    p.add_argument("--out-prefix", required=True,
                   help="Output path prefix for the landscape/dist/summary files.")
    return p.parse_args()


def load_eligible(ref_path):
    """Return (contig_name, length, 0-based indices of A/C/G/T bases)."""
    fa = pysam.FastaFile(ref_path)
    if len(fa.references) != 1:
        sys.exit("ERROR: --reference must contain exactly one sequence; found %d"
                 % len(fa.references))
    name = fa.references[0]
    seq = fa.fetch(name).upper()
    arr = np.frombuffer(seq.encode("ascii"), dtype="S1")
    valid = (arr == b"A") | (arr == b"C") | (arr == b"G") | (arr == b"T")
    return name, len(seq), np.where(valid)[0]


def main():
    args = parse_args()
    windows = [int(w) for w in args.windows.split(",")]

    contig, length, eligible = load_eligible(args.reference)
    n_sites = int(round(args.divergence * length))
    if n_sites > eligible.size:
        sys.exit("ERROR: requested %d sites but only %d eligible bases." % (n_sites, eligible.size))

    # --- build the field and sample divergent-site positions (same path as truth VCF)
    cap_ratio = (args.cap / args.divergence) if args.cap else None
    weights_full = dfield.build_weights(length, args.corr_len, args.sigma,
                                        args.field_seed, cap_ratio=cap_ratio)
    weights_elig = weights_full[eligible]
    rng = np.random.default_rng(args.sample_seed)
    picked_local = dfield.weighted_sample_without_replacement(weights_elig, n_sites, rng)
    pos = np.sort(eligible[picked_local])  # 0-based divergent-site positions

    # --- coarse binned landscape (for plotting), eligibility-aware -------------
    centers, rate, valid = dfield.binned_rate(pos, length, args.bin_size, eligible)
    realized_mean = pos.size / float(eligible.size)  # sites per eligible base
    corr_est = dfield.estimate_correlation_length(rate, valid, args.bin_size)
    # amplitude summaries computed over VALID (non-N) bins only
    vrate = rate[valid]
    cv = (vrate.std() / vrate.mean()) if vrate.mean() > 0 else float("nan")

    with open(args.out_prefix + "-landscape.tsv", "w") as f:
        f.write("bin_center\tlocal_divergence\tvalid\n")
        for c, r, v in zip(centers, rate, valid):
            f.write("%d\t%.8g\t%d\n" % (int(c), r, int(v)))

    # --- per-site local divergence at each requested window scale --------------
    with open(args.out_prefix + "-localdist.tsv", "w") as f:
        f.write("half_window\tmin\tp05\tp25\tmedian\tp75\tp95\tmax\tmean\n")
        for hw in windows:
            lr = dfield.windowed_local_rate(pos, hw)
            qs = np.quantile(lr, [0.0, 0.05, 0.25, 0.5, 0.75, 0.95, 1.0])
            f.write("%d\t%s\t%.6g\n" % (hw, "\t".join("%.6g" % q for q in qs), lr.mean()))

    # --- summary + validation --------------------------------------------------
    with open(args.out_prefix + "-summary.txt", "w") as f:
        f.write("# s05 divergence-field validation summary\n")
        f.write("contig\t%s\n" % contig)
        f.write("contig_length\t%d\n" % length)
        f.write("eligible_bases\t%d\n" % eligible.size)
        f.write("param_target_D\t%g\n" % args.divergence)
        f.write("param_corr_len_bp\t%g\n" % args.corr_len)
        f.write("param_sigma\t%g\n" % args.sigma)
        f.write("param_field_seed\t%d\n" % args.field_seed)
        f.write("param_sample_seed\t%d\n" % args.sample_seed)
        f.write("param_bin_size_bp\t%d\n" % args.bin_size)
        f.write("realized_n_div_sites\t%d\n" % pos.size)
        f.write("n_valid_bins\t%d\n" % int(valid.sum()))
        f.write("realized_mean_D\t%.8g\n" % realized_mean)
        f.write("realized_binned_D_min\t%.8g\n" % vrate.min())
        f.write("realized_binned_D_max\t%.8g\n" % vrate.max())
        f.write("realized_binned_D_cv\t%.6g\n" % cv)
        f.write("realized_corr_len_bp_1oe\t%s\n"
                % ("%.6g" % corr_est if corr_est == corr_est else "NA"))
        f.write("frac_bins_gt_2x_mean\t%.6g\n"
                % (np.mean(vrate > 2.0 * vrate.mean())))
        f.write("frac_bins_lt_half_mean\t%.6g\n"
                % (np.mean(vrate < 0.5 * vrate.mean())))

    with open(args.out_prefix + "-summary.txt.dict.txt", "w") as d:
        d.write(provenance.header(args.out_prefix + "-summary.txt"))
        d.write("Covers all three outputs: %s-{landscape.tsv,localdist.tsv,summary.txt}\n\n"
                % args.out_prefix)
        d.write("Purpose: validation of the s05 continuous autocorrelated divergence\n")
        d.write("landscape produced by divergence_field.build_weights (same code path\n")
        d.write("make_truth_vcf.py uses to place divergent SNP sites).\n\n")
        d.write("landscape.tsv:\n")
        d.write("  bin_center        : bp position of the center of a %d-bp bin\n" % args.bin_size)
        d.write("  local_divergence  : divergent-site density (sites/bp) in that bin\n\n")
        d.write("localdist.tsv (distribution of per-site local divergence):\n")
        d.write("  half_window       : +/- window (bp) over which local divergence is measured\n")
        d.write("  min..max,mean     : quantiles + mean of per-site local divergence\n\n")
        d.write("summary.txt (key = value):\n")
        d.write("  param_*           : field parameters used (D, corr_len, sigma, seeds, bin size)\n")
        d.write("  realized_mean_D   : n_div_sites / contig_length (pinned to target D by construction)\n")
        d.write("  realized_binned_D_cv : coefficient of variation of the binned local divergence\n")
        d.write("                         (amplitude of heterogeneity; ~0 when sigma=0)\n")
        d.write("  realized_corr_len_bp_1oe : lag at which the binned-track autocorrelation\n")
        d.write("                             first falls to 1/e (should track param_corr_len_bp)\n")
        d.write("  frac_bins_gt_2x_mean / lt_half_mean : fraction of bins far above/below the mean\n")

    # --- echo to stderr for the log -------------------------------------------
    sys.stderr.write(
        "# field: D=%g corr_len=%g sigma=%g -> n_sites=%d realized_mean_D=%.6g "
        "binned_CV=%.3g corr_len_est=%s\n"
        % (args.divergence, args.corr_len, args.sigma, pos.size, realized_mean, cv,
           ("%.4g" % corr_est if corr_est == corr_est else "NA")))


if __name__ == "__main__":
    main()
