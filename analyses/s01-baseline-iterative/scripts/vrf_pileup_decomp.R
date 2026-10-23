#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------------------------------
# s01 COMPANION DECOMPOSITION -- what changes when a SNP flips not-called -> called?
#
# f_alt (alt-support survival) is one-dimensional; it says nothing about the reference/contaminant
# reads at p, the depth, or the MAPQ, nor whether contaminant reads are MISMAPPED IN vs whether the
# SNP's own reads are MISMAPPED OUT. variant_read_fate.py --pileup adds the position-centric pileup at
# p (caller's-eye view) and a direct mismapped-IN count (mmin_*: reads aligned at p whose true span is
# elsewhere). Its mirror, mismapped-OUT, comes from the truth-centric ledger cells
# (mm_clip+mm_noclip+dc+un: reads that belong at p but map elsewhere / not at all).
#
# Crucially these groups are RIGHT-SKEWED, so this script reports and plots PER-SNP DISTRIBUTIONS
# (boxplots + medians), not just means -- a "mean 15 ref reads" flip SNP is really a median of ~5 with
# a long tail. It compares three fates (flip, stayed-called, stayed-not-called), by f_alt bin, and by
# divergence. Run from project root:
#   Rscript analyses/s01-baseline-iterative/scripts/vrf_pileup_decomp.R
# ---------------------------------------------------------------------------------------------------
suppressMessages({library(dplyr); library(tidyr); library(readr); library(ggplot2)})
source("analyses/lib/design.r")

PILE  <- "analyses/s01-baseline-iterative/results/vrf-pileup"
OUT   <- "analyses/s01-baseline-iterative/results/vrf-pileup-viz"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
lv    <- c("gatk", "bcftools", "freebayes")
BINS  <- c("0.0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0")
FATES <- c("stayed not-called", "flip (->called)", "stayed called")

files <- list.files(PILE, "vrf-pileup-.*\\.csv$", full.names = TRUE)
if (length(files) == 0) stop("no vrf-pileup CSVs in ", PILE)
nbps <- read_csv(file.path(OUT, "..", "neighborhood", "neighborhood-per-snp.csv"),
                 show_col_types = FALSE)              # per-SNP unrecovered-neighbor load
pl <- bind_rows(lapply(files, read_csv, show_col_types = FALSE)) %>%
  left_join(select(nbps, caller, divergence, snp_pos, iteration, nb_unrec = n_unrecovered),
            by = c("caller", "divergence", "snp_pos", "iteration")) %>%
  mutate(caller = factor(caller, levels = lv),
         f_alt  = ifelse(n_true_alt > 0, (ex_alt + sh_alt) / n_true_alt, NA_real_),
         called = recovered == 1)

mk <- function(k) inner_join(filter(pl, iteration == k), filter(pl, iteration == k + 1),
                             by = c("caller", "divergence", "snp_pos"), suffix = c("_a", "_b"))
pairs <- bind_rows(mk(1), mk(2)) %>% filter(!is.na(f_alt_a), !is.na(f_alt_b)) %>%
  mutate(grp = case_when(!called_a &  called_b ~ "flip (->called)",
                          called_a &  called_b ~ "stayed called",
                         !called_a & !called_b ~ "stayed not-called",
                          TRUE ~ "lost"),
         grp = factor(grp, levels = FATES),
         fbin_a = cut(f_alt_a, breaks = seq(0, 1, 0.2), include.lowest = TRUE, labels = BINS),
         mmout_a = mm_clip_a + mm_noclip_a + dc_a + un_a,     # mismapped-OUT (mirror of mmin)
         mmout_b = mm_clip_b + mm_noclip_b + dc_b + un_b)
flip <- filter(pairs, grp == "flip (->called)")

sep <- theme(strip.text.y = element_text(angle = 0), panel.spacing = unit(0.85, "lines"),
             panel.border = element_rect(color = "grey75", fill = NA, linewidth = 0.4),
             strip.background = element_rect(fill = "grey92", color = "grey75", linewidth = 0.4))

# --- Figure 1: PER-SNP DISTRIBUTIONS by fate (boxplots), D=0.10 ------------------------------------
METS <- c("depth at p", "f_alt", "ref reads at p", "mismapped-in ref", "mismapped-out",
          "mean MAPQ at p", "neighbor unrec SNPs (+/-200bp)")
d10 <- filter(pairs, grp %in% FATES, divergence == 0.10) %>% mutate(id = row_number())
dist <- bind_rows(
  d10 %>% transmute(id, caller, grp, metric = "depth at p",       before = depth_p_a,   after = depth_p_b),
  d10 %>% transmute(id, caller, grp, metric = "f_alt",            before = f_alt_a,      after = f_alt_b),
  d10 %>% transmute(id, caller, grp, metric = "ref reads at p",   before = ref_p_a,      after = ref_p_b),
  d10 %>% transmute(id, caller, grp, metric = "mismapped-in ref", before = mmin_ref_p_a, after = mmin_ref_p_b),
  d10 %>% transmute(id, caller, grp, metric = "mismapped-out",    before = mmout_a,      after = mmout_b),
  d10 %>% transmute(id, caller, grp, metric = "mean MAPQ at p",   before = mapq_p_a,     after = mapq_p_b),
  d10 %>% transmute(id, caller, grp, metric = "neighbor unrec SNPs (+/-200bp)", before = nb_unrec_a, after = nb_unrec_b)) %>%
  pivot_longer(c(before, after), names_to = "state", values_to = "value") %>%
  mutate(state = factor(state, levels = c("before", "after")), metric = factor(metric, levels = METS))
pDist <- ggplot(dist, aes(state, value, fill = caller)) +
  geom_boxplot(outlier.size = 0.2, outlier.alpha = 0.2, linewidth = 0.28,
               position = position_dodge(0.8), width = 0.7) +
  facet_grid(metric ~ grp, scales = "free_y") +
  scale_fill_manual(values = corecol(pal = "wilke", numcol = 3), name = "caller") +
  labs(x = "earlier -> later iteration", y = NULL,
       title = "Per-SNP pileup DISTRIBUTIONS by fate (s01, D=0.10)") + bartheme() + sep
ggsave(file.path(OUT, "pileup_decomp_dist.png"), pDist, width = 11, height = 9, dpi = 150)

# --- Figure 2: flip MEAN trend across divergences -------------------------------------------------
labA <- c(depth = "depth at p", ref = "ref reads at p", mminref = "mismapped-in ref",
          mmout = "mismapped-out", mapq = "mean MAPQ at p", falt = "f_alt")
pDiv <- flip %>% group_by(caller, divergence) %>%
  summarize(depth_a = mean(depth_p_a), depth_b = mean(depth_p_b), ref_a = mean(ref_p_a), ref_b = mean(ref_p_b),
            mminref_a = mean(mmin_ref_p_a), mminref_b = mean(mmin_ref_p_b),
            mmout_a = mean(mmout_a), mmout_b = mean(mmout_b), mapq_a = mean(mapq_p_a), mapq_b = mean(mapq_p_b),
            falt_a = mean(f_alt_a), falt_b = mean(f_alt_b), .groups = "drop") %>%
  pivot_longer(matches("_(a|b)$"), names_pattern = "(.+)_([ab])$", names_to = c("metric", ".value")) %>%
  mutate(metric = factor(labA[metric], levels = unname(labA))) %>%
  pivot_longer(c(a, b), names_to = "state", values_to = "value") %>%
  mutate(state = factor(ifelse(state == "a", "before", "after"), levels = c("before", "after"))) %>%
  ggplot(aes(state, value, group = caller, color = caller)) + geom_line(linewidth = 0.85) + geom_point(size = 2.2) +
  facet_grid(metric ~ divergence, scales = "free_y", labeller = labeller(divergence = label_both)) +
  scale_color_manual(values = corecol(pal = "wilke", numcol = 3), name = "caller") +
  labs(x = NULL, y = "group mean", title = "Flip SNPs: MEAN pileup before vs after, by divergence") + bartheme() + sep
ggsave(file.path(OUT, "pileup_decomp_flip.png"), pDiv, width = 9.5, height = 9, dpi = 150)

# --- Figure 3: flips by starting f_alt bin (means) ------------------------------------------------
labB <- c(falt = "f_alt (before vs after)", mminref = "mismapped-in ref reads at p")
pBin <- flip %>% filter(divergence == 0.10) %>% group_by(caller, fbin_a) %>%
  summarize(falt_a = mean(f_alt_a), falt_b = mean(f_alt_b),
            mminref_a = mean(mmin_ref_p_a), mminref_b = mean(mmin_ref_p_b), .groups = "drop") %>%
  pivot_longer(c(falt_a, falt_b, mminref_a, mminref_b), names_pattern = "(.+)_([ab])$",
               names_to = c("metric", ".value")) %>%
  mutate(metric = factor(labB[metric], levels = unname(labB))) %>%
  pivot_longer(c(a, b), names_to = "state", values_to = "value") %>%
  mutate(state = factor(ifelse(state == "a", "before", "after"), levels = c("before", "after"))) %>%
  ggplot(aes(fbin_a, value, color = caller, linetype = state, group = interaction(caller, state))) +
  geom_line(linewidth = 0.85) + geom_point(size = 2.1) + facet_wrap(~ metric, scales = "free_y") +
  scale_color_manual(values = corecol(pal = "wilke", numcol = 3), name = "caller") +
  scale_linetype_manual(values = c(before = "22", after = "solid"), name = NULL) +
  labs(x = "starting f_alt bin", y = NULL, title = "Flip SNPs by starting f_alt bin (D=0.10, means)") + bartheme() +
  theme(panel.spacing = unit(1.4, "lines"), panel.border = element_rect(color = "grey75", fill = NA, linewidth = 0.4),
        strip.background = element_rect(fill = "grey92", color = "grey75", linewidth = 0.4))
ggsave(file.path(OUT, "pileup_decomp_bybin.png"), pBin, width = 10, height = 4.2, dpi = 150)

# --- console report: three-fate MEDIAN [IQR] at D=0.10 (bcftools), showing the spread -------------
q3 <- function(x) sprintf("%.0f [%.0f-%.0f]", median(x), quantile(x, .25), quantile(x, .75))
cat("\n=== D=0.10 bcftools: three fates, median [IQR] (means overstate -- distributions are skewed) ===\n")
d10 %>% filter(caller == "bcftools") %>% group_by(grp) %>%
  summarize(n = n(), ref_before = q3(ref_p_a), ref_after = q3(ref_p_b),
            mmout_before = q3(mmout_a), mmout_after = q3(mmout_b),
            falt_before = sprintf("%.2f", median(f_alt_a)), .groups = "drop") %>%
  as.data.frame() %>% print(row.names = FALSE)
cat("\nwrote", OUT, "/{pileup_decomp_dist.png, pileup_decomp_flip.png, pileup_decomp_bybin.png}\n")
