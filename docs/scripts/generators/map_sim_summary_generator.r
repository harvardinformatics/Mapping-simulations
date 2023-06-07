############################################################
# For mapsims web, 10.21
# This generates the file "map_sim_summary.html"
############################################################

cat("Rendering map_sim_summary.rmd/html\n")

setwd("../generators/")
# Set the working directory.

Sys.setenv(RSTUDIO_PANDOC="C:/Program Files/RStudio/bin/pandoc/")
# Set the pandoc path

library(rmarkdown)
# Load markdown

output_dir = "../.."
# Set the output directory to be the docs directory, two up

render("../markdown/map_sim_summary.rmd", output_dir = output_dir, params = list(output_dir = output_dir), quiet = TRUE)
# Knit the html page