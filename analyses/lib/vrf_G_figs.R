# ---------------------------------------------------------------------------------------------
# Figure G: joint call-state x f_alt flow across iterations, from the aggregate -gtraj tables.
#
# Each `vrf-agg-<caller>-<div>-gtraj.csv` already holds, per condition, the count of SNPs per
# (status, f_alt-bin) trajectory across the three iterations -- i.e. exactly ggalluvial's
# aggregated ("lodes") form, so no per-SNP data is re-read here.
#
# This script produces TWO kinds of output from the same panel builder:
#   1. one standalone PNG per caller x divergence (each keeps its own legend, since it is
#      meant to be read on its own), and
#   2. a single combined grid, callers as rows x divergence as columns, assembled with
#      patchwork and `guides = "collect"` so the 15 identical f_alt legends collapse to ONE
#      shared legend on the right. Collecting only works when every panel's fill scale is
#      identical, which is why the scale below pins `limits = BIN` and `drop = FALSE` -- a
#      panel missing an f_alt bin would otherwise emit a shorter, non-mergeable guide.
#
# It is shared by s01 (SNP-only) and s02 (SNP + indel); the section-specific input/output
# directories and divergence list are passed on the command line by each section's driver.
#
# Usage:
#   Rscript analyses/lib/vrf_G_figs.R --in-dir <vrf-agg dir> --out-dir <viz dir> \
#           [--divs 0.02,0.04,0.06,0.08,0.10] [--callers gatk,bcftools,freebayes] \
#           [--label "SNP-only (s01)"]
# ---------------------------------------------------------------------------------------------

suppressMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readr)
  library(ggalluvial); library(patchwork)
})

# --- Display / layout parameters (all figure knobs live here) --------------------------------
BIN       <- c("0.0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0")  # f_alt bins, low -> high
FB_COL    <- setNames(hcl.colors(length(BIN), "Viridis"), BIN)         # bin -> colour
GAP_FRAC  <- 0.03      # white gap between the called / not-called boxes, as a fraction of N
BOX_HALFW <- 0.22      # half-width of the called / not-called outline boxes (x units)
STRAT_W   <- 0.44      # width of the alluvial strata
FLOW_A    <- 0.4       # alpha of the connecting flows
PANEL_W   <- 8.5       # standalone panel size (inches)
PANEL_H   <- 6.0
CELL_W    <- 6.2       # per-cell size inside the combined grid (inches)
CELL_H    <- 4.8
LEGEND_W  <- 1.4       # width reserved for the single collected legend (inches)
DPI       <- 100

# stratum factor order: not-called bins, then the spacer, then called bins (bottom-up)
LEV <- c(paste("not called", BIN, sep = " · "), "gap", paste("called", BIN, sep = " · "))

# --- Command-line arguments ------------------------------------------------------------------
# Minimal named-argument parser (no optparse dependency): --key value pairs.
parse_args <- function(argv) {
  if (length(argv) %% 2 != 0) stop("arguments must be --key value pairs")
  keys <- sub("^--", "", argv[c(TRUE, FALSE)])
  vals <- argv[c(FALSE, TRUE)]
  setNames(as.list(vals), keys)
}
a <- parse_args(commandArgs(trailingOnly = TRUE))
for (req in c("in-dir", "out-dir")) {
  if (is.null(a[[req]])) stop(sprintf("missing required argument --%s", req))
}
IN_DIR  <- a[["in-dir"]]
OUT_DIR <- a[["out-dir"]]
DIVS    <- as.numeric(strsplit(if (is.null(a[["divs"]])) "0.02,0.04,0.06,0.08,0.10"
                               else a[["divs"]], ",")[[1]])
CALLERS <- strsplit(if (is.null(a[["callers"]])) "gatk,bcftools,freebayes"
                    else a[["callers"]], ",")[[1]]
LABEL   <- if (is.null(a[["label"]])) "" else a[["label"]]

if (!dir.exists(IN_DIR)) stop(sprintf("input dir does not exist: %s", IN_DIR))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# --- Load the aggregate trajectory tables ------------------------------------------------------
gtraj_files <- list.files(IN_DIR, "gtraj.csv$", full.names = TRUE)
if (length(gtraj_files) == 0) stop(sprintf("no *-gtraj.csv files found in %s", IN_DIR))
gt <- bind_rows(lapply(gtraj_files, read_csv, col_types = cols()))

# --- Panel builder -----------------------------------------------------------------------------
# Returns the ggplot for one caller x divergence cell, or NULL when that condition is absent.
build_G <- function(caller_name, div_value) {

  g <- gt %>% filter(caller == caller_name, divergence == div_value,
                     fb1 %in% BIN, fb2 %in% BIN, fb3 %in% BIN)
  if (nrow(g) == 0) return(NULL)

  # wide (one row per trajectory) -> long lodes form: one row per trajectory per iteration
  lodes <- g %>%
    mutate(id = row_number()) %>%
    pivot_longer(c(s1, fb1, s2, fb2, s3, fb3),
                 names_to = c(".value", "iter"), names_pattern = "([a-z]+)([0-9])") %>%
    transmute(id, n = n_snps, iter,
              status = ifelse(s == "not_called", "not called", "called"),
              fb     = factor(fb, levels = BIN),
              state  = factor(paste(status, fb, sep = " · "), levels = LEV),
              is_gap = FALSE)

  # a fully transparent spacer stratum separating the called and not-called boxes
  gap_n <- round(GAP_FRAC * sum(g$n_snps))
  lodes <- bind_rows(lodes,
                     tibble(id = max(lodes$id) + 1, n = gap_n, iter = c("1", "2", "3"),
                            status = "gap", fb = factor(NA, levels = BIN),
                            state = factor("gap", levels = LEV), is_gap = TRUE))

  # outline boxes: called sits at the bottom, not-called stacked above the gap
  box <- lodes %>% filter(!is_gap) %>%
    group_by(iter, status) %>% summarize(n = sum(n), .groups = "drop") %>%
    group_by(iter) %>% mutate(ncall = sum(n[status == "called"])) %>% ungroup() %>%
    mutate(pos  = as.integer(iter), xmin = pos - BOX_HALFW, xmax = pos + BOX_HALFW,
           ymin = ifelse(status == "called", 0, ncall + gap_n),
           ymax = ifelse(status == "called", ncall, ncall + gap_n + n),
           yc   = (ymin + ymax) / 2)

  ggplot(lodes, aes(x = iter, stratum = state, alluvium = id, y = n,
                    fill = fb, alpha = is_gap)) +
    geom_flow() +
    # white fill first, so the flows do not show through the box interiors
    geom_rect(data = box, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              inherit.aes = FALSE, fill = "white") +
    geom_stratum(color = NA, width = STRAT_W) +
    geom_rect(data = box, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              inherit.aes = FALSE, fill = NA, color = "grey15", linewidth = 0.9) +
    geom_text(data = box, aes(x = pos, y = yc, label = status), inherit.aes = FALSE,
              size = 4.2, fontface = "bold", colour = "grey10") +
    scale_alpha_manual(values = c(`FALSE` = FLOW_A, `TRUE` = 0), guide = "none") +
    # limits + drop=FALSE keep this guide byte-identical across panels so patchwork can
    # collect the 15 copies into one shared legend
    scale_fill_manual(values = FB_COL, name = "f_alt bin\n(alt at p / true-alt reads)",
                      limits = BIN, drop = FALSE, na.value = "white", na.translate = FALSE) +
    scale_x_discrete(labels = c("iter1", "iter2", "iter3"), expand = c(.12, .12)) +
    labs(x = NULL, y = "SNPs (whole chromosome)",
         title = sprintf("%s, D=%s (%s SNPs)", caller_name, format(div_value),
                         format(sum(g$n_snps), big.mark = ","))) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(), axis.text.y = element_blank())
}

# --- Build every panel, write standalone PNGs, keep them for the grid ---------------------------
# Row-major order (caller rows x divergence columns) is what patchwork's wrap_plots expects.
panels <- list()
for (cc in CALLERS) {
  for (dd in DIVS) {
    p <- build_G(cc, dd)
    if (is.null(p)) {
      warning(sprintf("no gtraj rows for caller=%s divergence=%s -- panel skipped", cc, dd))
      next
    }
    ggsave(file.path(OUT_DIR, sprintf("G_%s_D%s.png", cc, format(dd))),
           p, width = PANEL_W, height = PANEL_H, dpi = 120)
    panels[[length(panels) + 1]] <- p
  }
}
if (length(panels) == 0) stop("no panels were built -- nothing to assemble")

# --- Assemble the combined grid with ONE shared legend ------------------------------------------
# guides = "collect" merges the identical f_alt guides; the `&` theme then places that single
# collected legend on the right of the whole grid rather than beside each panel.
grid <- wrap_plots(panels, nrow = length(CALLERS), guides = "collect") &
  theme(legend.position = "right")

grid <- grid + plot_annotation(
  title = paste0("Figure G: joint state (called / not-called x f_alt) across iterations",
                 if (nzchar(LABEL)) paste0(" -- ", LABEL) else ""),
  subtitle = "rows = caller, columns = divergence; colour = f_alt bin (shared legend, right)",
  theme = theme(plot.title = element_text(face = "bold", size = 16),
                plot.subtitle = element_text(size = 12, colour = "grey30")))

grid_png <- file.path(OUT_DIR, "G_grid.png")
ggsave(grid_png, grid,
       width  = length(DIVS) * CELL_W + LEGEND_W,
       height = length(CALLERS) * CELL_H,
       dpi = DPI, limitsize = FALSE)

cat(sprintf("wrote %d standalone panels + %s (%d x %d cells, 1 shared legend)\n",
            length(panels), grid_png, length(CALLERS), length(DIVS)))
