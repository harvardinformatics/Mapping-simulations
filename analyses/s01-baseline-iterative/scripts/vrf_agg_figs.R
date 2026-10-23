# Whole-chromosome figures from the aggregate variant-read-fate tables (run_vrf_aggregate.sh output).
# Reads results/vrf-agg/vrf-agg-<caller>-<div>-hcomp.csv (at-p composition) and -gtraj.csv (G flow).
# Produces, into results/vrf-viz-full/:
#   full_altatp_notcalled.png  -- headline: why SNPs stay uncalled (assembler vs alignment), all SNPs
#   full_atp_composition.png   -- at-p read-fate composition at D=0.10, caller x call-status
suppressMessages({library(ggplot2);library(dplyr);library(tidyr);library(readr)})
R   <- "analyses/s01-baseline-iterative/results/vrf-agg"
OUT <- "analyses/s01-baseline-iterative/results/vrf-viz-full"
lv  <- c("gatk","bcftools","freebayes")
h <- bind_rows(lapply(list.files(R, "hcomp.csv$", full.names=TRUE), read_csv, col_types=cols())) %>%
  mutate(caller=factor(caller, levels=lv))

## Headline: alt-maps-at-p fraction for NOT-CALLED SNPs vs divergence, by caller (iter1 solid, iter3 dashed)
pA <- h %>% filter(status=="not_called", iteration %in% c(1,3)) %>%
  mutate(alt_frac=alt_at_p/n_reads) %>%
  ggplot(aes(divergence, alt_frac, color=caller, linetype=factor(iteration), group=interaction(caller,iteration))) +
  geom_line(linewidth=0.9) + geom_point(size=2.4) +
  scale_y_continuous(labels=scales::percent, limits=c(0,1)) +
  scale_color_manual(values=corecol(pal="wilke", numcol=3), name="caller") +
  scale_linetype_manual(values=c(`1`="solid",`3`="22"), name="iteration") +
  labs(x="divergence (D)", y="alt maps at p  (fraction of covering reads)",
       title="Why do SNPs stay uncalled? (whole chromosome, every hom-div SNP)",
       subtitle="not-called SNPs: gatk ~80% (base IS presented -> assembler-limited); bcftools/freebayes ~20-30% (base absent -> alignment-limited)") +
  bartheme() + theme(plot.subtitle=element_text(size=9, colour="#555555"))
ggsave(file.path(OUT,"full_altatp_notcalled.png"), pA, width=9, height=5.2, dpi=125)

## Composition at D=0.10: at-p read fate, caller (cols) x call status (rows), across iterations
cats <- c("alt_at_p","ref_at_p","err_at_p","clipped","deleted","mapped_elsewhere","unmapped")
nice <- c(alt_at_p="alt maps at p", ref_at_p="ref maps at p", err_at_p="other base (error)",
          clipped="clipped over p", deleted="deleted at p", mapped_elsewhere="mapped elsewhere", unmapped="unmapped")
fcol <- setNames(c("#1b7837","#2166ac","#92c5de","#e08214","#8073ac","#b2182b","#4d4d4d"), unname(nice[cats]))
pB <- h %>% filter(divergence==0.10) %>%
  pivot_longer(all_of(cats), names_to="atp", values_to="reads") %>%
  mutate(frac=reads/n_reads, atp=factor(nice[atp], levels=unname(nice[cats])),
         status=recode(status, called="called", not_called="not called")) %>%
  ggplot(aes(factor(iteration), frac, fill=atp)) + geom_col(width=0.75, color="grey30", linewidth=0.12) +
  facet_grid(status~caller) + scale_fill_manual(values=fcol, name="what is at position p") +
  scale_y_continuous(labels=scales::percent) +
  labs(x="iteration", y="fraction of covering reads",
       title="What maps at position p — D=0.10, whole chromosome, by caller and call outcome") +
  bartheme()
ggsave(file.path(OUT,"full_atp_composition.png"), pB, width=10, height=5.5, dpi=125)
cat("done\n")
