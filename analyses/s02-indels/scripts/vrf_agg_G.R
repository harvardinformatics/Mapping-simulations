# ---------------------------------------------------------------------------------------------
# s02 (SNP + indel) driver for figure G -- the joint call-state x f_alt flow across iterations.
#
# Mirror of the s01 driver: all plotting logic lives in the shared builder
# `analyses/lib/vrf_G_figs.R`, and this file only pins the s02 input/output directories and the
# divergence sweep. It writes one standalone PNG per caller x divergence plus the combined
# 3 x 5 grid (G_grid.png) carrying a single shared f_alt legend.
#
# Before this script existed the s02 grid was produced by an uncommitted scratch script plus an
# unrecorded image-montage step, so the figure embedded in the synthesis notebook could not be
# regenerated from the repository. This driver closes that gap.
#
# Run from the project root:
#   Rscript analyses/s02-indels/scripts/vrf_agg_G.R
# ---------------------------------------------------------------------------------------------

IN_DIR  <- "analyses/s02-indels/results/vrf-agg"
OUT_DIR <- "analyses/s02-indels/results/vrf-viz-G"
DIVS    <- "0.02,0.04,0.06,0.08,0.10"   # same homogeneous-divergence sweep as s01, indels on
LABEL   <- "SNP + indel (s02)"

# hand off to the shared builder, which does the loading, panel construction and assembly
status <- system2("Rscript",
                  c("analyses/lib/vrf_G_figs.R",
                    "--in-dir",  shQuote(IN_DIR),
                    "--out-dir", shQuote(OUT_DIR),
                    "--divs",    shQuote(DIVS),
                    "--label",   shQuote(LABEL)))

# fail loudly: a non-zero exit from the builder must not look like a successful run
if (status != 0) stop(sprintf("analyses/lib/vrf_G_figs.R failed with exit status %d", status))
