############################################################
# For Mapping sims Evolution talk
# June 2023
# SNP accuracy after 3 iterations
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
  vcf_comp_file = here("data", "mm39-30X-0.005h-3i-vcf-comparison.tsv.gz")
  vcf_comp = read_tsv(vcf_comp_file, 
                            col_names=c("coverage", "divergence", "heterozygosity", "iteration", "chr", "pos", "ref", "golden.div", "prev.iter", 
                                        "iter.gt", "iter.allele.1", "iter.allele.2", "minimap.gt", "minimap.allele.1", "minimap.allele.2")) %>%
    filter(divergence != 0.00)
  # Load the SNP data from first iteration
}

############################################################

tp = vcf_comp %>% filter(golden.div != ref & minimap.gt == "hom.alt" & minimap.allele.1 == golden.div) %>% 
  group_by(divergence, iteration) %>% 
  summarize(n=n(), class="tp")
fp = vcf_comp %>% filter(golden.div == ref & minimap.gt == "hom.alt") %>% 
  group_by(divergence, iteration) %>% 
  summarize(n=n(), class="fp")
fn = vcf_comp %>% filter(golden.div != ref & minimap.gt == "hom.ref" & minimap.allele.1 == ref) %>% 
  group_by(divergence, iteration) %>% 
  summarize(n=n(), class="fn")

snp_classes = rbind(tp, fp, fn)
snp_classes = snp_classes %>% group_by(divergence, iteration) %>% mutate(total.sites = sum(n)) %>% mutate(prop = n / total.sites)
snp_classes$class = factor(snp_classes$class, levels=c("fn", "tp", "fp"))

cols = c("fp"=corecol(pal="wilke", numcol=1, offset=1), "tp"="#bda988", "fn"=corecol(pal="wilke", numcol=1))

snps_p = ggplot(snp_classes, aes(x=iteration, y=prop, fill=class)) +
  geom_bar(stat="identity", position = position_fill(reverse = TRUE)) +
  #geom_text(size=2.3, position=position_stack(vjust=0.5), color="#f2f2f2") +
  xlab("Mapping iteration") +
  ylab("Proportion of variants") +
  ggtitle("By divergence") +
  scale_y_continuous(expand=c(0,0), limits=c(0,1), breaks=seq(0, 1, by=0.2)) +
  scale_fill_manual(labels=c("False negative", "True positive", "False positive"), values=cols) +
  facet_wrap(~divergence) +
  bartheme() +
  theme(legend.position="bottom") +
  coord_flip()
print(snps_p)

############################################################