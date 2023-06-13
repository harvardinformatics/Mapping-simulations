############################################################
# For Mapping sims Evolution talk
# June 2023
# SNP accuracy after 1 iteration
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
  vcf_comp_iter1 = read_tsv(vcf_comp_file, 
                            col_names=c("coverage", "divergence", "heterozygosity", "iteration", "chr", "pos", "ref", "golden.div", "prev.iter", 
                                            "iter.gt", "iter.allele.1", "iter.allele.2", "minimap.gt", "minimap.allele.1", "minimap.allele.2")) %>%
    filter(divergence != 0.00)
  # Load the SNP data from first iteration
}

############################################################

tp = vcf_comp_iter1 %>% filter(golden.div != ref & iter.gt == "hom.alt" & iter.allele.1 == golden.div) %>% 
  group_by(divergence, iteration) %>% 
  summarize(n=n(), class="tp")
fp = vcf_comp_iter1 %>% filter(golden.div == ref & iter.gt == "hom.alt") %>% 
  group_by(divergence, iteration) %>% 
  summarize(n=n(), class="fp")
fn = vcf_comp_iter1 %>% filter(golden.div != ref & iter.gt == "hom.prev" & iter.allele.1 == ref) %>% 
  group_by(divergence, iteration) %>% 
  summarize(n=n(), class="fn")
# Get sites of different classes from the iteration VCF (GATK)

snp_classes = rbind(tp, fp, fn)
# Combine the various snp classes

snp_classes = snp_classes %>% group_by(divergence, iteration) %>% 
  mutate(total.sites = sum(n)) %>% 
  mutate(prop = n / total.sites)
# Calculate the proportion of each type of site

snp_classes$class = factor(snp_classes$class, levels=c("fn", "tp", "fp"))
cols = c("fp"=corecol(pal="wilke", numcol=1, offset=1), "tp"="#bda988", "fn"=corecol(pal="wilke", numcol=1))
# Factorize the classes and get some colors

snps_p = ggplot(snp_classes, aes(x=divergence, y=prop, fill=class)) +
  geom_bar(stat="identity", position = position_fill(reverse = TRUE)) +
  #geom_text(size=2.3, position=position_stack(vjust=0.5), color="#f2f2f2") +
  xlab("Simulated divergence") +
  ylab("Proportion of variants") +
  scale_x_continuous(breaks=seq(0, 0.1, by=0.02)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_fill_manual(labels=c("False negative", "True positive", "False positive"), values=cols) +
  bartheme() +
  theme(legend.position="bottom") +
  coord_flip()
print(snps_p)

############################################################

if(save_fig){
  outfilename = here("evolution", "figs", "03-snp-accuracy-iter1.png")
  ggsave(filename=outfilename, snps_p, width=8, height=5, units="in")
}

