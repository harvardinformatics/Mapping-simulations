# ===================================================================================================
# plot_divergence_distribution.R -- QC figure for the s05 realistic capped config: the per-read and
# per-site divergence distributions that justify the mapping-vs-calling plot's x-axis.
#
# Usage (from project root, via run_divergence_distribution.sh which extracts the inputs):
#   Rscript plot_divergence_distribution.R <readdiv.txt> <locald150.txt> <out.png>
#     readdiv.txt   : one per-read divergence per line = NM / read_length  (grid = 1/readlen = 1/150)
#     locald150.txt : one per-site local_d_150 per line = SNPs in +/-150 bp / 300  (grid = 1/300)
#
# IMPORTANT — binning: both quantities are DISCRETE (integer / fixed denominator), so a free-running
# bin count aliases against the value grid and produces spurious gaps (per-read) / comb spikes
# (per-site). We bin EACH panel on its own value grid with a bin centred on every allowed value
# (per-read grid = 1/readlen = 1/150; per-site grid = 1/300), via binwidth = grid and center = 0.
# Then each bar is exactly one achievable value and the histogram is artifact-free.
# ===================================================================================================
source("analyses/lib/design.r")
suppressMessages({library(ggplot2); library(dplyr); library(patchwork)})

args <- commandArgs(trailingOnly = TRUE)
readdiv_path <- args[1]; locald_path <- args[2]; out_path <- args[3]

rd <- scan(readdiv_path, quiet = TRUE)      # per-read divergence, grid = 1/150
sd <- scan(locald_path,  quiet = TRUE)      # per-site local_d_150, grid = 1/300
mrd <- mean(rd); msd <- mean(sd)
BW_READ <- 1 / 150     # per-read grid: divergence = NM / 150
BW_SITE <- 1 / 300     # per-site grid: local_d_150 = SNP count / 300

# Panel A: per-read divergence; mark field mean (0.03) and the actual read mean
pA <- ggplot(data.frame(x = rd), aes(x)) +
  geom_histogram(binwidth = BW_READ, center = 0, fill = corecol(numcol = 1), color = NA) +
  geom_vline(xintercept = 0.03, linetype = "dashed", color = "grey35") +
  geom_vline(xintercept = mrd, linetype = "solid", color = "black", linewidth = 0.5) +
  annotate("text", x = 0.03, y = Inf, label = "field mean\n0.03", hjust = 1.1, vjust = 1.4,
           size = 3, color = "grey35", lineheight = 0.9) +
  annotate("text", x = mrd, y = Inf, label = sprintf("read mean %.3f", mrd), hjust = -0.06,
           vjust = 1.6, size = 3.1, color = "black") +
  scale_x_continuous(limits = c(-BW_READ, 0.35)) +
  labs(title = "A. Per-read divergence (truth BAM, NM / read length)",
       subtitle = sprintf("n=%s reads · mean=%.3f · p99=%.3f · max=%.3f  (1 bar = 1 NM value, grid 1/150)",
                          format(length(rd), big.mark = ","), mrd, quantile(rd, .99), max(rd)),
       x = "divergent bases per read (fraction)", y = "reads") + bartheme()

# Panel B: per-site local divergence; mark field mean, site mean, field cap
pB <- ggplot(data.frame(x = sd), aes(x)) +
  geom_histogram(binwidth = BW_SITE, center = 0, fill = corecol(numcol = 2)[2], color = NA) +
  geom_vline(xintercept = 0.03, linetype = "dashed", color = "grey35") +
  geom_vline(xintercept = msd, linetype = "solid",  color = "black", linewidth = 0.6) +
  geom_vline(xintercept = 0.20, linetype = "dashed", color = "firebrick") +
  annotate("text", x = 0.03, y = Inf, label = "field mean\n0.03", hjust = 1.1, vjust = 1.4,
           size = 3, color = "grey35", lineheight = 0.9) +
  annotate("text", x = msd, y = Inf, label = sprintf("site mean %.3f", msd), hjust = -0.06,
           vjust = 1.4, size = 3.2, color = "black") +
  annotate("text", x = 0.20, y = Inf, label = "field cap 0.20", hjust = -0.06, vjust = 1.4,
           size = 3, color = "firebrick") +
  scale_x_continuous(limits = c(-BW_SITE, 0.35)) +
  labs(title = "B. Per-site read-scale local divergence (local_d_150, at divergent SNPs)",
       subtitle = sprintf("n=%s sites · mean=%.3f (>2x field mean: size-biased toward dense regions) · p99=%.3f · p99.9=%.3f · max=%.3f",
                          format(length(sd), big.mark = ","), msd, quantile(sd, .99), quantile(sd, .999), max(sd)),
       x = "local divergence in ±150 bp window (fraction)", y = "divergent SNP sites") + bartheme()

ggsave(out_path, pA / pB, width = 10, height = 8, dpi = 135)
cat(sprintf("wrote %s  (per-read mean=%.4f, per-site mean=%.4f)\n", out_path, mrd, msd))
