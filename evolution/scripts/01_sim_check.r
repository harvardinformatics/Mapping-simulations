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

if(read_data){
  vcf_comp_file = here("data", "mm39-30X-0.005h-3i-vcf-comparison-iter1.tsv.gz")
  sim_vars_div_iter1 = read_tsv(vcf_comp_file, 
                      col_names=c("coverage", "divergence", "heterozygosity", "iteration", "chr", "pos", "ref", "golden.div", "prev.iter", 
                                  "iter.gt", "iter.allele.1", "iter.allele.2", "minimap.gt", "minimap.allele.1", "minimap.allele.2")) %>%
    filter(ref != golden.div) %>%
    group_by(divergence) %>%
    summarize(n=n()) %>%
    mutate(perc.var.sites=n/total_len)
  # Load the SNP data from first iteration
}

############################################################

div_sim_p = ggplot(sim_vars_div_iter1, aes(x=divergence, y=perc.var.sites)) +
  geom_point(aes(color="Simulations"), size=4) +
  geom_line(aes(color="Simulations"), show.legend=F) +
  geom_abline(aes(slope=1, intercept=0, color="1:1 line"), size=1, linetype="dashed") +
  scale_x_continuous(limits=c(0,0.105), breaks=seq(0,0.10,0.02), expand=c(0,0)) +
  scale_y_continuous(limits=c(0,0.105), breaks=seq(0,0.10,0.02), expand=c(0,0)) +
  scale_color_manual(values=c("1:1 line"="#666666", "Simulations"="#920000")) +
  xlab("Specified simulated divergence") +
  ylab("Actual % of sites with SNPs") +
  bartheme() +
  theme(legend.spacing.y = unit(0.5, 'cm'),
        axis.text=element_text(size=12), 
        axis.title=element_text(size=14)) +
  guides(color = guide_legend(byrow = TRUE, 
                override.aes = list(linetype = c("dashed", "blank"),
                                    shape = c(NA, 19))))
print(div_sim_p)

############################################################

if(save_fig){
  outfilename = here("evolution", "figs", "01-sim-check.png")
  ggsave(filename=outfilename, div_sim_p, width=6, height=5, units="in")
}


