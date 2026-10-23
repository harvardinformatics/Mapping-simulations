# ===================================================================================================
# decomp_figs.R -- shared companion-decomposition figures, so s01/s02/s03 render the SAME analysis.
#
# Source this from a notebook (knit from project root; it assumes design.r is already sourced for
# bartheme()/corecol()). Then:
#   pairs <- decomp_load(pile_dir, neigh_per_snp_csv)   # per-SNP pileup + neighborhood, joined
#   decomp_distfig(pairs, 0.10)                         # per-SNP DISTRIBUTIONS by fate (boxplots)
#   decomp_meantrend(pairs, "flip (→called)")           # mean trend across divergences, one fate
#   decomp_bybin(pairs)                                 # flips by starting f_alt bin (D=0.10)
#   decomp_neighfig(neigh_summary_csv)                  # assemblable-neighborhood figure
#   decomp_gatksig_clean(pairs); decomp_gatksig_delta(pairs)   # the two gatk-signature tables
# Works identically for hom-div (s01/s02) and het (s03) pileup CSVs -- the columns are the same; only
# the interpretation TEXT (written in each notebook) differs.
# ===================================================================================================
suppressMessages({library(ggplot2); library(dplyr); library(tidyr); library(readr)})

.FATES <- c("stayed not-called", "flip (→called)", "stayed called")
.BINS  <- c("0.0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0")
.CLR   <- function() corecol(pal = "wilke", numcol = 3)

# which failure mode each pileup metric diagnoses: the read-placement metrics say whether the
# alt reads ARRIVED (mapping); the assemblable-neighborhood load says whether gatk's reassembler
# could act (calling). Used only when a figure is called with tag = TRUE (default off, so
# s01/s02/s03 are unchanged).
.MCAT <- c("depth at p" = "mapping", "f_alt" = "mapping", "ref reads at p" = "mapping",
           "mismapped-in ref" = "mapping", "mismapped-out" = "mapping",
           "mean MAPQ at p" = "mapping", "neighbor unrec SNPs (±200bp)" = "calling",
           "neighbor unrec indels (±200bp)" = "calling")
.tagmet <- function(m) paste0(m, ifelse(.MCAT[m] == "calling", "  · CALLING ·", "  · mapping ·"))
.sep <- function() theme(strip.text.y = element_text(angle = 0), panel.spacing = unit(0.85, "lines"),
                         panel.border = element_rect(color = "grey75", fill = NA, linewidth = 0.4),
                         strip.background = element_rect(fill = "grey92", color = "grey75", linewidth = 0.4))

# --- load: per-SNP pileup CSVs joined to the per-SNP neighborhood load, then paired across iters -----
decomp_load <- function(pile_dir, neigh_per_snp) {
  nbps <- read_csv(neigh_per_snp, col_types = cols())
  pl <- bind_rows(lapply(list.files(pile_dir, "vrf-pileup-.*\\.csv$", full.names = TRUE),
                         read_csv, col_types = cols())) %>%
    left_join(select(nbps, caller, divergence, snp_pos, iteration, nb_unrec = n_unrecovered,
                     dplyr::any_of(c(nb_unrec_indel = "n_unrecovered_indel"))),
              by = c("caller", "divergence", "snp_pos", "iteration")) %>%
    mutate(caller = factor(caller, levels = c("gatk", "bcftools", "freebayes")),
           f_alt = ifelse(n_true_alt > 0, (ex_alt + sh_alt) / n_true_alt, NA_real_),
           called = recovered == 1)
  mk <- function(k) inner_join(filter(pl, iteration == k), filter(pl, iteration == k + 1),
                               by = c("caller", "divergence", "snp_pos"), suffix = c("_a", "_b"))
  bind_rows(mk(1), mk(2)) %>% filter(!is.na(f_alt_a), !is.na(f_alt_b)) %>%
    mutate(grp = case_when(!called_a &  called_b ~ "flip (→called)",
                            called_a &  called_b ~ "stayed called",
                           !called_a & !called_b ~ "stayed not-called",
                            TRUE ~ "lost"),
           grp = factor(grp, levels = .FATES),
           fbin_a = cut(f_alt_a, breaks = seq(0, 1, 0.2), include.lowest = TRUE, labels = .BINS),
           mmout_a = mm_clip_a + mm_noclip_a + dc_a + un_a,     # mismapped-OUT (mirror of mmin)
           mmout_b = mm_clip_b + mm_noclip_b + dc_b + un_b)
}

# --- per-SNP DISTRIBUTIONS by fate (boxplots), one divergence -----------------------------------------
decomp_distfig <- function(pairs, dd, tag = FALSE) {
  METS <- c("depth at p", "f_alt", "ref reads at p", "mismapped-in ref", "mismapped-out",
            "mean MAPQ at p", "neighbor unrec SNPs (±200bp)")
  d1 <- filter(pairs, grp %in% .FATES, divergence == dd) %>% mutate(id = row_number())
  bind_rows(
    d1 %>% transmute(id, caller, grp, metric = "depth at p",       before = depth_p_a,   after = depth_p_b),
    d1 %>% transmute(id, caller, grp, metric = "f_alt",            before = f_alt_a,      after = f_alt_b),
    d1 %>% transmute(id, caller, grp, metric = "ref reads at p",   before = ref_p_a,      after = ref_p_b),
    d1 %>% transmute(id, caller, grp, metric = "mismapped-in ref", before = mmin_ref_p_a, after = mmin_ref_p_b),
    d1 %>% transmute(id, caller, grp, metric = "mismapped-out",    before = mmout_a,      after = mmout_b),
    d1 %>% transmute(id, caller, grp, metric = "mean MAPQ at p",   before = mapq_p_a,     after = mapq_p_b),
    d1 %>% transmute(id, caller, grp, metric = "neighbor unrec SNPs (±200bp)", before = nb_unrec_a, after = nb_unrec_b)) %>%
    pivot_longer(c(before, after), names_to = "state", values_to = "value") %>%
    mutate(state = factor(state, levels = c("before", "after")),
           metric = factor(if (tag) .tagmet(metric) else metric,
                           levels = if (tag) .tagmet(METS) else METS)) %>%
    ggplot(aes(state, value, fill = caller)) +
    geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.2, linewidth = 0.28,
                 position = position_dodge(0.8), width = 0.7) +
    facet_grid(metric ~ grp, scales = "free_y") +
    scale_fill_manual(values = .CLR(), name = "caller") +
    labs(x = "earlier → later iteration", y = NULL,
         title = sprintf("Per-SNP pileup DISTRIBUTIONS by fate — D=%s", format(dd))) +
    bartheme() + .sep()
}

# --- mean trend across divergences for ONE fate ------------------------------------------------------
decomp_meantrend <- function(pairs, g, tag = FALSE) {
  labA <- c(depth = "depth at p", ref = "ref reads at p", mminref = "mismapped-in ref",
            mmout = "mismapped-out", mapq = "mean MAPQ at p", falt = "f_alt",
            nbunrec = "neighbor unrec SNPs (±200bp)")
  filter(pairs, grp == g) %>% group_by(caller, divergence) %>%
    summarize(depth_a = mean(depth_p_a), depth_b = mean(depth_p_b), ref_a = mean(ref_p_a), ref_b = mean(ref_p_b),
              mminref_a = mean(mmin_ref_p_a), mminref_b = mean(mmin_ref_p_b),
              mmout_a = mean(mmout_a), mmout_b = mean(mmout_b), mapq_a = mean(mapq_p_a), mapq_b = mean(mapq_p_b),
              falt_a = mean(f_alt_a), falt_b = mean(f_alt_b),
              nbunrec_a = mean(nb_unrec_a), nbunrec_b = mean(nb_unrec_b), .groups = "drop") %>%
    pivot_longer(matches("_(a|b)$"), names_pattern = "(.+)_([ab])$", names_to = c("metric", ".value")) %>%
    mutate(metric = factor(if (tag) .tagmet(labA[metric]) else labA[metric],
                           levels = if (tag) .tagmet(unname(labA)) else unname(labA))) %>%
    pivot_longer(c(a, b), names_to = "state", values_to = "value") %>%
    mutate(state = factor(ifelse(state == "a", "earlier", "later"), levels = c("earlier", "later"))) %>%
    ggplot(aes(state, value, group = caller, color = caller)) +
    geom_line(linewidth = 0.85) + geom_point(size = 2.2) +
    facet_grid(metric ~ divergence, scales = "free_y", labeller = labeller(divergence = label_both)) +
    scale_color_manual(values = .CLR(), name = "caller") +
    labs(x = NULL, y = "group mean", title = sprintf("MEAN pileup by divergence — %s", g)) +
    bartheme() + .sep()
}

# --- neighbour unrecovered SNP vs INDEL load, per-site by fate (same format as decomp_distfig /
#     decomp_meantrend, restricted to the two assemblable-neighbourhood rows). Needs nb_unrec_indel,
#     so decomp_load must have been given an s02-style neighborhood file with the indel column. -------
decomp_neighbox <- function(pairs, dd, tag = FALSE) {
  METS <- c("neighbor unrec SNPs (±200bp)", "neighbor unrec indels (±200bp)")
  d1 <- filter(pairs, grp %in% .FATES, divergence == dd) %>% mutate(id = row_number())
  bind_rows(
    d1 %>% transmute(id, caller, grp, metric = "neighbor unrec SNPs (±200bp)",   before = nb_unrec_a,       after = nb_unrec_b),
    d1 %>% transmute(id, caller, grp, metric = "neighbor unrec indels (±200bp)", before = nb_unrec_indel_a, after = nb_unrec_indel_b)) %>%
    pivot_longer(c(before, after), names_to = "state", values_to = "value") %>%
    mutate(state = factor(state, levels = c("before", "after")),
           metric = factor(if (tag) .tagmet(metric) else metric,
                           levels = if (tag) .tagmet(METS) else METS)) %>%
    ggplot(aes(state, value, fill = caller)) +
    geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.2, linewidth = 0.28,
                 position = position_dodge(0.8), width = 0.7) +
    facet_grid(metric ~ grp, scales = "free_y") +
    scale_fill_manual(values = .CLR(), name = "caller") +
    labs(x = "earlier → later iteration", y = NULL,
         title = sprintf("Neighbour unrecovered SNPs vs indels — per-site DISTRIBUTIONS by fate, D=%s (s02)", format(dd))) +
    bartheme() + .sep()
}

decomp_neighmean <- function(pairs, g, tag = FALSE) {
  labA <- c(nbsnp = "neighbor unrec SNPs (±200bp)", nbindel = "neighbor unrec indels (±200bp)")
  filter(pairs, grp == g) %>% group_by(caller, divergence) %>%
    summarize(nbsnp_a = mean(nb_unrec_a), nbsnp_b = mean(nb_unrec_b),
              nbindel_a = mean(nb_unrec_indel_a), nbindel_b = mean(nb_unrec_indel_b), .groups = "drop") %>%
    pivot_longer(matches("_(a|b)$"), names_pattern = "(.+)_([ab])$", names_to = c("metric", ".value")) %>%
    mutate(metric = factor(if (tag) .tagmet(labA[metric]) else labA[metric],
                           levels = if (tag) .tagmet(unname(labA)) else unname(labA))) %>%
    pivot_longer(c(a, b), names_to = "state", values_to = "value") %>%
    mutate(state = factor(ifelse(state == "a", "earlier", "later"), levels = c("earlier", "later"))) %>%
    ggplot(aes(state, value, group = caller, color = caller)) +
    geom_line(linewidth = 0.85) + geom_point(size = 2.2) +
    facet_grid(metric ~ divergence, scales = "free_y", labeller = labeller(divergence = label_both)) +
    scale_color_manual(values = .CLR(), name = "caller") +
    labs(x = NULL, y = "group mean",
         title = sprintf("MEAN neighbour unrecovered SNP vs indel load by divergence — %s (s02)", g)) +
    bartheme() + .sep()
}

# --- flips by starting f_alt bin (D=0.10) ------------------------------------------------------------
decomp_bybin <- function(pairs) {
  labB <- c(falt = "f_alt (before vs after)", mminref = "mismapped-in ref reads at p")
  filter(pairs, grp == "flip (→called)", divergence == 0.10) %>% group_by(caller, fbin_a) %>%
    summarize(falt_a = mean(f_alt_a), falt_b = mean(f_alt_b),
              mminref_a = mean(mmin_ref_p_a), mminref_b = mean(mmin_ref_p_b), .groups = "drop") %>%
    pivot_longer(c(falt_a, falt_b, mminref_a, mminref_b),
                 names_pattern = "(.+)_([ab])$", names_to = c("metric", ".value")) %>%
    mutate(metric = factor(labB[metric], levels = unname(labB))) %>%
    pivot_longer(c(a, b), names_to = "state", values_to = "value") %>%
    mutate(state = factor(ifelse(state == "a", "before (uncalled)", "after (called)"),
                          levels = c("before (uncalled)", "after (called)"))) %>%
    ggplot(aes(fbin_a, value, color = caller, linetype = state, group = interaction(caller, state))) +
    geom_line(linewidth = 0.85) + geom_point(size = 2.1) + facet_wrap(~ metric, scales = "free_y") +
    scale_color_manual(values = .CLR(), name = "caller") +
    scale_linetype_manual(values = c("before (uncalled)" = "22", "after (called)" = "solid"), name = NULL) +
    labs(x = "starting f_alt bin (iteration where still uncalled)", y = NULL,
         title = "Flip variants by starting f_alt bin (D=0.10, group means)") +
    bartheme() + theme(panel.spacing = unit(1.4, "lines"),
                       panel.border = element_rect(color = "grey75", fill = NA, linewidth = 0.4),
                       strip.background = element_rect(fill = "grey92", color = "grey75", linewidth = 0.4))
}

# --- STATIC calling-side covariates by single-shot fate (sibling of the mean-trend plot) --------------
# The decomp_distfig/decomp_meantrend plots above are an ITERATION analysis: a metric watched
# before -> after as iteration cleans the pileup. The metrics here are STATIC site properties
# (local sequence/variant context) that do NOT change across iterations, so they are shown on the
# single-shot three-way fate instead: recovered / calling-limited (alt reads present, uncalled) /
# mapping-limited (alt reads never placed). Same visual family (bartheme, caller colours, mapping vs
# CALLING tag) so it reads as a companion to the decomposition means. Input is the aggregated
# calling-metrics-by-fate.csv (divergence, caller, outcome, metric, n, mean, se).
#
# Empirically confirmed roles (used for the tag): homopolymer_frac separates gatk's CALLING-limited
# sites (the assembler-resolvability signal); indel/snp density separate MAPPING-limited sites for
# every caller; het_density is flat at these het rates.
.CMET <- c(homopolymer_frac = "homopolymer frac (≥5bp)", indel_density = "indel density (/bp)",
           snp_density = "div-SNP density (/bp)", het_density = "het-SNP density (/bp)")
.CMCAT <- c(homopolymer_frac = "calling", het_density = "calling",
            indel_density = "mapping", snp_density = "mapping")
.ctagmet <- function(m) paste0(.CMET[m], ifelse(.CMCAT[m] == "calling", "  · CALLING ·", "  · mapping ·"))

decomp_callingfig <- function(byfate_csv, tag = TRUE) {
  ord_fate <- c("mapping-limited", "calling-limited", "recovered")
  ord_met  <- if (tag) .ctagmet(names(.CMET)) else unname(.CMET)
  read_csv(byfate_csv, col_types = cols()) %>%
    filter(metric %in% names(.CMET)) %>%
    mutate(caller = factor(caller, levels = c("gatk", "bcftools", "freebayes")),
           outcome = factor(outcome, levels = ord_fate),
           metric = factor(if (tag) .ctagmet(metric) else .CMET[metric], levels = ord_met)) %>%
    ggplot(aes(outcome, mean, group = caller, color = caller)) +
    geom_line(linewidth = 0.85) + geom_point(size = 2.2) +
    geom_errorbar(aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se), width = 0.2, linewidth = 0.3) +
    facet_grid(metric ~ divergence, scales = "free_y", labeller = labeller(divergence = label_both)) +
    scale_color_manual(values = .CLR(), name = "caller") +
    labs(x = "single-shot fate", y = "group mean (95% CI)",
         title = "Static calling-side covariates by SNP fate (s02 indels)") +
    bartheme() + .sep() + theme(axis.text.x = element_text(angle = 20, hjust = 1))
}

# --- STATIC calling covariates by ITERATION transition fate (companion to decomp_callingfig) ----------
# Same static metrics, but grouped by the iteration-transition fate (stayed not-called / flip->called
# / stayed called) instead of the single-shot fate. The metrics don't change across iterations, so
# there is no before/after -- but the transition-fate GROUPING is valid and answers a different
# question: does static sequence/variant context PREDICT which sites iteration can rescue? Input is
# the aggregated calling-metrics-by-iterfate.csv (caller, divergence, fate, metric, n, mean, se).
decomp_callingfig_iter <- function(byiterfate_csv, tag = TRUE) {
  ord_fate <- c("stayed not-called", "flip (→called)", "stayed called")
  ord_met  <- if (tag) .ctagmet(names(.CMET)) else unname(.CMET)
  read_csv(byiterfate_csv, col_types = cols()) %>%
    filter(metric %in% names(.CMET)) %>%
    mutate(caller = factor(caller, levels = c("gatk", "bcftools", "freebayes")),
           fate = factor(fate, levels = ord_fate),
           metric = factor(if (tag) .ctagmet(metric) else .CMET[metric], levels = ord_met)) %>%
    ggplot(aes(fate, mean, group = caller, color = caller)) +
    geom_line(linewidth = 0.85) + geom_point(size = 2.2) +
    geom_errorbar(aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se), width = 0.2, linewidth = 0.3) +
    facet_grid(metric ~ divergence, scales = "free_y", labeller = labeller(divergence = label_both)) +
    scale_color_manual(values = .CLR(), name = "caller") +
    labs(x = "iteration fate (iter1 → iter2)", y = "group mean (95% CI)",
         title = "Static calling-side covariates by ITERATION fate (s02 indels)") +
    bartheme() + .sep() + theme(axis.text.x = element_text(angle = 20, hjust = 1))
}

# --- assemblable-neighborhood figure (from the aggregate summary CSV) ---------------------------------
decomp_neighfig <- function(neigh_summary, dd = 0.10) {
  fcol <- c("stayed called" = "#4d9221", "flip (→called)" = "#e08214", "stayed not-called" = "#b2182b")
  read_csv(neigh_summary, col_types = cols()) %>%
    filter(divergence == dd, fate %in% c("flip", "stayed_notcalled", "stayed_called")) %>%
    transmute(caller = factor(caller, levels = c("gatk", "bcftools", "freebayes")),
              fate = factor(fate, levels = c("stayed_called", "flip", "stayed_notcalled"), labels = names(fcol)),
              before = mean_unrecovered_before, after = mean_unrecovered_after) %>%
    pivot_longer(c(before, after), names_to = "state", values_to = "unrec") %>%
    mutate(state = factor(state, levels = c("before", "after"))) %>%
    ggplot(aes(state, unrec, color = fate, group = fate)) +
    geom_line(linewidth = 0.9) + geom_point(size = 2.8) + facet_wrap(~ caller) +
    scale_color_manual(values = fcol, name = "fate") +
    labs(x = "earlier → later iteration", y = "unrecovered neighbor\nhom-div SNPs (±200 bp)",
         title = sprintf("Assemblable neighborhood: local unrecovered-SNP load, D=%s", format(dd))) +
    bartheme() + theme(plot.title = element_text(size = 11))
}

# --- the two gatk-signature tables (return data frames; the notebook kables them) --------------------
decomp_gatksig_clean <- function(pairs, dd = 0.10, het = FALSE) {
  df <- filter(pairs, divergence == dd, grp == "flip (→called)")
  # "call-ready before the flip": for hom-div, ~0 ref reads (ref = contamination); for HET the ref
  # haplotype supplies ~half the reads by design, so ref_p is not a cleanliness signal -- use the
  # foreign mismapped-in ref instead and a lower alt-support bar.
  df <- if (het)
    dplyr::mutate(df, already_clean = f_alt_a >= 0.6 & mmin_ref_p_a <= 2 & mapq_p_a >= 50)
  else
    dplyr::mutate(df, already_clean = f_alt_a >= 0.8 & ref_p_a <= 2 & mapq_p_a >= 50 & mmin_ref_p_a <= 2)
  df %>% group_by(caller) %>%
    summarize(`n flips` = n(),
              `pileup already clean before flip` = sprintf("%.0f%%", 100 * mean(already_clean)),
              .groups = "drop")
}
decomp_gatksig_delta <- function(pairs, dd = 0.10) {
  filter(pairs, divergence == dd, grp == "flip (→called)") %>% group_by(caller) %>%
    summarize(`Δ depth` = round(mean(abs(depth_p_b - depth_p_a)), 1),
              `Δ ref` = round(mean(abs(ref_p_b - ref_p_a)), 1),
              `Δ mismap-in` = round(mean(abs(mmin_ref_p_b - mmin_ref_p_a)), 1),
              `Δ MAPQ` = round(mean(abs(mapq_p_b - mapq_p_a)), 1),
              `Δ f_alt` = round(mean(abs(f_alt_b - f_alt_a)), 2), .groups = "drop")
}
