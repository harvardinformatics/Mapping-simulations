############################################################
# For Mapping sims Evolution talk
# June 2023
# Check simulated variants
# Gregg Thomas
############################################################

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

library(tidyverse)
library(here)
source(here("evolution/scripts/lib/design.r"))

############################################################

read_data = T

save_fig = T

total_len = 90720763

############################################################

vcf_comp_file = here("data", "mm39-30X-0.005h-3i-vcf-comparison.tsv.gz")
vcf_comp = read_tsv(vcf_comp_file, col_names=c("coverage", "divergence", "heterozygosity", "iteration", "chr", "pos", "ref", "golden.div", "prev.iter", "iter.gt", "iter.allele.1", "iter.allele.2", "minimap.gt", "minimap.allele.1", "minimap.allele.2"))
# Load the SNP data (~15gb, takes time)

sim_vars_div_iter1 = vcf_comp %>% filter(iteration == 1 & ref != golden.div) %>%
  group_by(divergence) %>%
  summarize(n=n()) %>%
  mutate(perc.var.sites=n/total_len)
# Subset to only the first iteration data

############################################################

div_sim_p = ggplot(sim_vars_div_iter1, aes(x=divergence, y=perc.var.sites)) +
  geom_point(size=4, color="#666666") +
  geom_abline(slope=1, intercept=0, size=1, linetype="dashed", color=corecol(numcol=1, pal="wilke")) +
  scale_y_continuous(limits=c(0,0.10), breaks=seq(0,0.10,0.02)) +
  xlab("Simulated divergence") +
  ylab("% of sites simulated with SNPs") +
  bartheme()
print(div_sim_p)

############################################################