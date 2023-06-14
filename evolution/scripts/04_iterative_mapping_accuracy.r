############################################################
# For Mapping sims Evolution talk
# June 2023
# Mapping accuracy after 3 iterations
# Gregg Thomas
############################################################

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

library(tidyverse)
library(here)
source(here("evolution/scripts/lib/design.r"))

############################################################

read_data = T

save_fig = F

total_len = 90720763

############################################################

if(read_data){
  mapping_file = here("data", "mm39-30X-0.005h-3i-bam-summary.csv")
  bam_stats = read.csv(mapping_file, header=T) %>% filter(divergence != 0.00) %>% mutate(total.reads = exact.map + close.map + mismapped + diff.chr + unmapped + missing.in.query) %>%
    mutate(exact.map.prop = exact.map / total.reads, close.map.prop = close.map / total.reads, mismapped.prop = mismapped / total.reads, 
           diff.chr.prop = diff.chr / total.reads, unmapped.prop = unmapped / total.reads, missing.prop = missing.in.query / total.reads)
  # Read the bam compare file and calculate the percent of each type of read from the first iteration of mappping
}

############################################################

read_classes = c("exact.map.prop", "close.map.prop", "mismapped.prop", "diff.chr.prop", "unmapped.prop", "missing.prop")
read_class_names = c("Exact map", "Close map", "Mismapped", "Different chrome", "Unmapped", "Missing")
# The various classes of reads

long_bam_data = bam_stats %>% 
  select(iteration, divergence, exact.map.prop, close.map.prop, mismapped.prop, diff.chr.prop, unmapped.prop, missing.prop) %>% 
  pivot_longer(cols=c("exact.map.prop", "close.map.prop", "mismapped.prop", "diff.chr.prop", "unmapped.prop", "missing.prop"))
# Convert to long format

long_bam_data$name = factor(long_bam_data$name, levels=read_classes)
# Factorize the read classes

cols = corecol(pal="wilke", numcol=6, offset=1)
# Get some colors

reads_p = ggplot(long_bam_data, aes(x=iteration, y=value, fill=name)) +
  geom_bar(stat="identity", position=position_stack(reverse=TRUE)) +
  #geom_text(size=2.3, position=position_stack(vjust=0.5), color="#f2f2f2") +
  xlab("Mapping iteration") +
  ylab("Proportion of reads") +
  scale_y_continuous(expand=c(0,0), limits=c(0,1), breaks=seq(0, 1, by=0.2)) +
  scale_fill_manual(labels=read_class_names, values=cols) +
  facet_wrap(~divergence) +
  bartheme() +
  theme(legend.position="bottom",
        panel.spacing=unit(1.5, "lines"),
        axis.text=element_text(size=12), 
        axis.title=element_text(size=14),
        plot.margin = unit(c(1,0.5,0.25,0.25), "cm"),
        #strip.background=element_blank(),
        strip.text=element_text(size=14)) +
  coord_flip()
print(reads_p)

############################################################

if(save_fig){
  outfilename = here("evolution", "figs", "04-iterative-mapping-accuracy.png")
  ggsave(filename=outfilename, reads_p, width=8, height=5, units="in")
}
