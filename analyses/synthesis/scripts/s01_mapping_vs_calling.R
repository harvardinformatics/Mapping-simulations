#!/usr/bin/env Rscript
# s01_mapping_vs_calling.R -- classify every sampled s01 divergent SNP, per (caller, divergence,
# iteration), as recovered / calling-limited / mapping-limited, from the s01 per-site pileup ledger
# (analyses/s01-baseline-iterative/results/vrf-pileup/vrf-pileup-<caller>-<D>.csv, written by
# analyses/lib/variant_read_fate.py --pileup). Used by the synthesis notebook's mapping-vs-calling
# figure.
#
# Codifies the ad hoc command that first made the table (2026-08-13); the classification rule and
# output format are unchanged.
#
# Usage (from the project root):
#   conda run -n renv Rscript analyses/synthesis/scripts/s01_mapping_vs_calling.R [OUT]
#   OUT defaults to analyses/synthesis/s01-mapping-vs-calling.csv
#   --dict-only [OUT] : only (re)write the data dictionary of an existing table, as a backfill

suppressMessages({ library(dplyr); library(readr) })
source("analyses/lib/provenance.r")

# --- settings ------------------------------------------------------------------------------------
SCRIPT   <- "analyses/synthesis/scripts/s01_mapping_vs_calling.R"
PILE_DIR <- "analyses/s01-baseline-iterative/results/vrf-pileup"
MIN_ALT  <- 2        # an unrecovered site with >= MIN_ALT alt reads at p is calling-limited

args <- commandArgs(trailingOnly = TRUE)
dict_only <- length(args) > 0 && args[1] == "--dict-only"
if (dict_only) args <- args[-1]
OUT <- if (length(args) > 0) args[1] else "analyses/synthesis/s01-mapping-vs-calling.csv"

# --- data dictionary -----------------------------------------------------------------------------
dict_body <- c(
  "Per (caller, divergence, iteration): how the sampled s01 divergent (1|1) SNPs split into three",
  "outcomes. Input: every vrf-pileup-<caller>-<D>.csv in",
  sprintf("%s (30,000 sampled SNPs per caller x D, iterations 1-3).", PILE_DIR),
  "",
  "Classification of each (SNP, iteration):",
  "  recovered       : recovered == 1 (called at this iteration)",
  sprintf("  calling-limited : not recovered, but >= %d reads present the alt base at the SNP in the", MIN_ALT),
  "                    mapped BAM (alt_p): the evidence arrived and the caller did not emit it",
  sprintf("  mapping-limited : not recovered and alt_p < %d: the alt-carrying reads never arrived", MIN_ALT),
  "",
  "Fields:",
  "  caller     : variant caller (gatk | bcftools | freebayes)",
  "  divergence : simulated hom-div divergence D (s01 SNP-only, seed 42)",
  "  iteration  : mapping iteration (1 = single-shot)",
  "  class      : recovered | calling-limited | mapping-limited",
  "  n          : SNPs in that class",
  "  frac       : n / all sampled SNPs of that (caller, divergence, iteration); the three sum to 1")

if (dict_only) {
  prov_write_dictionary(OUT, dict_body, SCRIPT,
                        command = sprintf("conda run -n renv Rscript %s", SCRIPT), backfill = TRUE,
                        notes = c("first made 2026-08-13 by an equivalent ad hoc R command (same rule and output format); this script reproduces it byte for byte (checked 2026-10-26)"))
  quit(save = "no")
}

# --- classify ------------------------------------------------------------------------------------
files <- list.files(PILE_DIR, "vrf-pileup-.*\\.csv$", full.names = TRUE)
if (length(files) == 0) stop("no vrf-pileup-*.csv in ", PILE_DIR)
pl <- bind_rows(lapply(files, read_csv, col_types = cols()))

cls <- pl %>%
  mutate(class = case_when(recovered == 1 ~ "recovered",
                           alt_p >= MIN_ALT ~ "calling-limited",
                           TRUE ~ "mapping-limited")) %>%
  group_by(caller, divergence, iteration, class) %>% summarize(n = n(), .groups = "drop") %>%
  group_by(caller, divergence, iteration) %>% mutate(frac = n / sum(n)) %>% ungroup()

write_csv(cls, OUT)
prov_write_dictionary(OUT, dict_body, SCRIPT)
cat("wrote", nrow(cls), "rows to", OUT, "; divergences:", paste(sort(unique(cls$divergence)), collapse = ","), "\n")
