#!/usr/bin/env python3
"""
divergence_field.py -- shared engine for the s05 continuous autocorrelated
divergence landscape.

WHY THIS EXISTS
---------------
The homogeneous simulations (s01-s04) place divergent SNP sites UNIFORMLY along
the contig: every eligible base is equally likely to be a fixed difference. Real
divergence is not uniform -- it varies along the genome at many scales (features,
genes, recombination/isochore structure). s05 replaces uniform placement with a
spatially AUTOCORRELATED density field so that divergent sites cluster and thin
out the way real substitutions do, while keeping the genome-wide mean divergence
pinned to the target D (so s05 is directly comparable to the homogeneous run at
the same D).

This module is imported by BOTH:
  * make_truth_vcf.py       -- to actually place the divergent SNP sites, and
  * make_divergence_field.py -- to generate/validate the landscape and emit the
                                per-window local-rate table for plots.
Sharing one implementation guarantees the validation landscape and the truth VCF
use identical field math.

THE MODEL (two knobs)
---------------------
  1. Latent field z(x): standard Gaussian white noise smoothed with a Gaussian
     kernel of standard deviation `corr_len` (in bp), then standardized to zero
     mean / unit variance. Smoothing is done in the frequency domain (FFT) with
     the analytic Gaussian transfer function, so no large kernel array is built.
     `corr_len` = CORRELATION LENGTH = the scale of heterogeneity (knob 1).
  2. Weights w(x) = exp(sigma * z(x)) -- a log-normal transform. `sigma` = the
     AMPLITUDE of heterogeneity (knob 2): sigma = 0 gives w == 1 everywhere
     (uniform placement, i.e. the homogeneous control); larger sigma produces
     stronger high-divergence spikes and low-divergence troughs.

  Divergent-site placement then draws EXACTLY round(D * L) positions WITHOUT
  replacement with probability proportional to w(x). Because the count is fixed,
  the genome-wide mean divergence equals the target D by construction (no
  re-centering needed); the field only controls WHERE the sites concentrate.

EDGE HANDLING
-------------
  The FFT smoothing is circular (periodic), so the two contig ends are treated
  as adjacent. This introduces a minor artifact only within ~4*corr_len of the
  single seam at the contig end; for chr19 (~61 Mb) and corr_len <= a few hundred
  kb this affects a negligible fraction of the sequence. It is documented here
  and reported by make_divergence_field.py rather than silently smoothed over.
"""

import numpy as np


def build_latent_field(length, corr_len, seed):
    """Return a length-`length` float64 array: standardized (mean 0, var 1)
    Gaussian field autocorrelated at scale `corr_len` (bp).

    White noise smoothed by a Gaussian kernel of SD = corr_len, implemented via
    the analytic frequency-domain transfer function of a Gaussian (avoids
    materializing a multi-Mb kernel). corr_len is in samples == bp here.
    """
    if corr_len <= 0:
        raise ValueError("corr_len must be > 0 (got %r)" % corr_len)
    rng = np.random.default_rng(seed)
    g = rng.standard_normal(length)

    # rfft frequencies in cycles/sample; Gaussian transfer fn for SD=corr_len.
    freqs = np.fft.rfftfreq(length)
    transfer = np.exp(-2.0 * (np.pi ** 2) * (corr_len ** 2) * (freqs ** 2))
    z = np.fft.irfft(np.fft.rfft(g) * transfer, n=length)

    # Smoothing shrinks the variance; restore to zero mean / unit variance so
    # `sigma` has a consistent meaning regardless of corr_len.
    z -= z.mean()
    sd = z.std()
    if sd > 0.0:
        z /= sd
    return z


def build_weights(length, corr_len, sigma, seed, cap_ratio=None):
    """Return per-position placement weights, normalized to mean 1.

    w(x) = exp(sigma * z(x)) / mean, so a position's weight is its local
    divergence relative to the genome mean: a site there has expected local
    divergence D_mean * w(x).

    cap_ratio bounds the heterogeneity to a biological ceiling. Because expected
    local divergence = D_mean * w, clipping w at `cap_ratio` caps expected local
    divergence at cap_ratio * D_mean (e.g. cap_ratio = 0.20 / 0.06 keeps a mean-
    0.06 field's hot regions near local divergence 0.20 instead of the log-normal
    tail running to 0.5+). Clipping shifts the mean slightly, so the REALIZED
    ceiling overshoots a little and is validated/tuned via make_divergence_field
    rather than assumed exact. Troughs are untouched (conserved regions -> ~0).

    sigma == 0 -> all weights 1.0 (uniform placement == homogeneous control).
    """
    if sigma == 0.0:
        return np.ones(length, dtype=np.float64)
    z = build_latent_field(length, corr_len, seed)
    w = np.exp(sigma * z)
    w /= w.mean()                       # weight ratio ~ local-D / mean-D
    if cap_ratio is not None:
        if cap_ratio <= 0:
            raise ValueError("cap_ratio must be > 0 (got %r)" % cap_ratio)
        w = np.minimum(w, cap_ratio)
    return w


def weighted_sample_without_replacement(weights, n, rng):
    """Draw `n` distinct indices from range(len(weights)) with probability
    proportional to `weights`, without replacement.

    Uses the Efraimidis-Spirakis exponential-key method: for each item draw
    key_i = Exp(1) / w_i and keep the n SMALLEST keys. This is exact weighted
    sampling without replacement and runs in O(len(weights)) via argpartition
    -- important because len(weights) is the whole eligible set (tens of
    millions of bases).
    """
    size = weights.shape[0]
    if n >= size:
        return np.arange(size)
    if np.any(weights <= 0.0):
        raise ValueError("all weights must be > 0 for weighted sampling")
    keys = rng.standard_exponential(size) / weights
    idx = np.argpartition(keys, n)[:n]
    return idx


def windowed_local_rate(pos_sorted, half_window):
    """Per-site local divergence = (number of divergent sites within +/-
    half_window bp of the site, including itself) / (2 * half_window).

    `pos_sorted` must be a sorted 1-D array of 0-based divergent-site positions.
    Returns a float array aligned to pos_sorted. This is the divergence a read
    covering the site actually experiences at the given window scale.
    """
    lo = np.searchsorted(pos_sorted, pos_sorted - half_window, side="left")
    hi = np.searchsorted(pos_sorted, pos_sorted + half_window, side="right")
    count = hi - lo  # inclusive of the site itself
    return count / (2.0 * half_window)


def windowed_density_at(query, source_sorted, half_window):
    """Density of `source_sorted` positions within +/- half_window bp of each
    position in `query` (sites/bp). Generalizes windowed_local_rate to evaluate at
    positions that are NOT themselves in the source set (e.g. indel positions
    scored against the divergent-SNP field). `source_sorted` must be sorted."""
    lo = np.searchsorted(source_sorted, query - half_window, side="left")
    hi = np.searchsorted(source_sorted, query + half_window, side="right")
    return (hi - lo) / (2.0 * half_window)


def binned_rate(positions, length, bin_size, eligible_pos):
    """Coarse landscape for plotting/validation: local divergence per bin,
    measured as divergent sites per ELIGIBLE (non-N) base -- NOT per genomic bp.

    Dividing by eligible bases (rather than bin_size) is essential: chr-scale N
    gaps otherwise read as megabase-wide zero-divergence blocks and swamp the
    field's real structure. Bins whose eligible fraction is below `min_frac` are
    marked invalid (they are mostly N and their density is meaningless).

    Returns (bin_center, rate, valid) where `valid` is a boolean mask.
    """
    min_frac = 0.5
    nbins = int(np.ceil(length / bin_size))
    edges = np.arange(nbins + 1, dtype=np.float64) * bin_size
    site_counts, _ = np.histogram(positions, bins=edges)
    elig_counts, _ = np.histogram(eligible_pos, bins=edges)
    centers = (edges[:-1] + edges[1:]) / 2.0
    valid = elig_counts >= (min_frac * bin_size)
    rate = np.zeros(nbins, dtype=np.float64)
    rate[valid] = site_counts[valid] / elig_counts[valid]
    return centers, rate, valid


def _longest_valid_run(valid):
    """Return (start, stop) of the longest contiguous True run in `valid`."""
    best_len = best_start = 0
    i = 0
    n = valid.shape[0]
    while i < n:
        if valid[i]:
            j = i
            while j < n and valid[j]:
                j += 1
            if j - i > best_len:
                best_len, best_start = j - i, i
            i = j
        else:
            i += 1
    return best_start, best_start + best_len


def estimate_correlation_length(rate, valid, bin_size):
    """Estimate the realized correlation length (bp) of the binned local-
    divergence track as the lag at which its autocorrelation first falls to 1/e.

    Computed on the LONGEST CONTIGUOUS run of valid (non-N) bins so that N gaps
    do not inject spurious long-range correlation. Returns np.nan if it never
    drops to 1/e within that run. Validation reporting only.
    """
    lo, hi = _longest_valid_run(valid)
    seg = rate[lo:hi]
    if seg.shape[0] < 10:
        return np.nan
    x = seg - seg.mean()
    var = np.dot(x, x)
    if var == 0.0:
        return np.nan
    max_lag = min(seg.shape[0] - 1, 2000)
    thresh = 1.0 / np.e
    for lag in range(1, max_lag + 1):
        ac = np.dot(x[:-lag], x[lag:]) / var
        if ac < thresh:
            return lag * bin_size
    return np.nan
