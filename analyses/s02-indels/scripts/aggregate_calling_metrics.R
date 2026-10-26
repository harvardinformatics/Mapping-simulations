# ===================================================================================================
# aggregate_calling_metrics.R -- turn the per-SNP static calling covariates (calling-metrics.csv,
# from run_calling_metrics.sh) into the two BY-FATE summary tables the synthesis figures read.
#
# Run from the project root:
#   conda run -n renv Rscript analyses/s02-indels/scripts/aggregate_calling_metrics.R
#
# Inputs:
#   analyses/s02-indels/results/calling-metrics/calling-metrics.csv   (per-SNP static covariates)
#   analyses/s02-indels/results/vrf-pileup/vrf-pileup-<caller>-<D>.csv (per-SNP call fate + alt_p)
#
# Outputs (mean + SE of each covariate, per caller x divergence x fate x metric):
#   calling-metrics-by-fate.csv      -- SINGLE-SHOT fate (iter1): recovered / calling-limited
#                                       (alt_p>=2 present but uncalled) / mapping-limited (alt absent)
#   calling-metrics-by-iterfate.csv  -- ITERATION fate (iter1->iter2): stayed not-called /
#                                       flip(->called) / stayed called
# The metrics are STATIC (identical across iterations); the two tables just group them two ways.
# ===================================================================================================
suppressMessages({library(dplyr); library(readr); library(tidyr)})

DIR   <- "analyses/s02-indels/results/calling-metrics"
PILE  <- "analyses/s02-indels/results/vrf-pileup"
METS  <- c("homopolymer_frac", "indel_density", "snp_density", "het_density")

cm <- read_csv(file.path(DIR, "calling-metrics.csv"), col_types = cols())
pile_files <- list.files(PILE, "vrf-pileup-.*\\.csv$", full.names = TRUE)
pl_all <- bind_rows(lapply(pile_files, read_csv, col_types = cols()))

# helper: mean + SE of each metric over a fate grouping
summarise_by <- function(df, fate_col) {
  df %>%
    pivot_longer(all_of(METS), names_to = "metric", values_to = "v") %>%
    group_by(caller, divergence, !!sym(fate_col), metric) %>%
    summarize(n = n(), mean = mean(v), se = sd(v) / sqrt(n()), .groups = "drop") %>%
    rename(fate_value = !!sym(fate_col))
}

# --- SINGLE-SHOT fate (iteration 1) -------------------------------------------------------------------
single <- pl_all %>%
  filter(iteration == 1) %>%
  mutate(outcome = case_when(recovered == 1 ~ "recovered",
                             alt_p >= 2    ~ "calling-limited",
                             TRUE          ~ "mapping-limited")) %>%
  select(caller, divergence, snp_pos, outcome) %>%
  inner_join(cm, by = c("divergence", "snp_pos"))
single_agg <- summarise_by(single, "outcome") %>% rename(outcome = fate_value)
write_csv(single_agg, file.path(DIR, "calling-metrics-by-fate.csv"))

# --- ITERATION transition fate (iter1 -> iter2) -------------------------------------------------------
iter <- pl_all %>%
  filter(iteration %in% c(1, 2)) %>%
  select(caller, divergence, snp_pos, iteration, recovered) %>%
  pivot_wider(names_from = iteration, values_from = recovered, names_prefix = "it") %>%
  filter(!is.na(it1), !is.na(it2)) %>%
  mutate(fate = case_when(it1 == 0 & it2 == 1 ~ "flip (→called)",
                          it1 == 1 & it2 == 1 ~ "stayed called",
                          it1 == 0 & it2 == 0 ~ "stayed not-called",
                          TRUE                ~ "lost")) %>%
  filter(fate != "lost") %>%                       # drop the rare called->uncalled transition
  inner_join(cm, by = c("divergence", "snp_pos"))
iter_agg <- summarise_by(iter, "fate") %>% rename(fate = fate_value)
write_csv(iter_agg, file.path(DIR, "calling-metrics-by-iterfate.csv"))

cat(sprintf("wrote calling-metrics-by-fate.csv (%d rows) and calling-metrics-by-iterfate.csv (%d rows)\n",
            nrow(single_agg), nrow(iter_agg)))
