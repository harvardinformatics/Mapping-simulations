############################################################
# For Mapping sims Evolution talk
# June 2023
# Divergence est. across genome
# Gregg Thomas
############################################################

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

library(tidyverse)
library(here)
source(here("evolution/scripts/lib/design.r"))
library(zoo)
# For the rolling mean function

############################################################

read_data = T

save_fig = T

############################################################

window_file = here("data", "window-snps-30X-summary.tab.gz")
windows = read_tsv(window_file, col_names=c("chr", "start", "end", "snps", "coverage", "divergence", "heterozygosity", "iteration"))
windows = windows %>% filter(divergence != 0.00 & chr %in% c("18")) %>% 
  mutate(length=end-start) %>% 
  mutate(est.div=snps/length) %>% 
  mutate(iteration = replace(iteration, iteration == "golden", "Golden"))
# Read the data and calculate divergence

############################################################

div_p = ggplot(windows, aes(x=start, y=rollmean(est.div, 75, na.pad = TRUE, align = "right"), color=iteration)) +
  geom_line(linewidth=0.75, alpha=0.75) +
  facet_wrap(~divergence) +
  xlab("Position along chromosome") +
  ylab("Avg. divergence from reference genome") +
  scale_color_manual(name="Iteration", values=corecol(pal="wilke", numcol=4)) +
  bartheme() +
  theme(legend.title=element_text(size=12),
        axis.text.x = element_text(angle=25, hjust=0.5),
        panel.spacing=unit(1.5, "lines"),
        axis.text=element_text(size=10), 
        axis.title=element_text(size=14),
        plot.margin = unit(c(1,0.5,0.25,0.25), "cm"),
        strip.text=element_text(size=14)) +
  guides(color = guide_legend(override.aes = list(linewidth = 2, alpha=1)))

print(div_p)

############################################################

if(save_fig){
  outfilename = here("evolution", "figs", "07-iterative-divergence-genome.png")
  ggsave(filename=outfilename, div_p, width=8, height=5, units="in")
}

