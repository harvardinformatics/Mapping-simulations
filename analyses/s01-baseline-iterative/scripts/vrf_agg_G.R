# ---------------------------------------------------------------------------------------------
# s01 (SNP-only) driver for figure G -- the joint call-state x f_alt flow across iterations.
#
# All plotting logic lives in the shared builder `analyses/lib/vrf_G_figs.R`, which is also used
# by s02; this file only pins the s01 input/output directories and the divergence sweep. It
# writes one standalone PNG per caller x divergence plus the combined 3 x 5 grid (G_grid.png)
# carrying a single shared f_alt legend.
#
# Run from the project root:
#   Rscript analyses/s01-baseline-iterative/scripts/vrf_agg_G.R
# ---------------------------------------------------------------------------------------------

IN_DIR  <- "analyses/s01-baseline-iterative/results/vrf-agg"
OUT_DIR <- "analyses/s01-baseline-iterative/results/vrf-viz-full"
DIVS    <- "0.02,0.04,0.06,0.08,0.10"   # the s01 homogeneous-divergence sweep shown in the grid
LABEL   <- "SNP-only (s01)"

# hand off to the shared builder, which does the loading, panel construction and assembly
status <- system2("Rscript",
                  c("analyses/lib/vrf_G_figs.R",
                    "--in-dir",  shQuote(IN_DIR),
                    "--out-dir", shQuote(OUT_DIR),
                    "--divs",    shQuote(DIVS),
                    "--label",   shQuote(LABEL)))

# fail loudly: a non-zero exit from the builder must not look like a successful run
if (status != 0) stop(sprintf("analyses/lib/vrf_G_figs.R failed with exit status %d", status))
