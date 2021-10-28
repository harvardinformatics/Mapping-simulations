############################################################
# For mapsims web, 10.21
# This generates the file "variants.html"
############################################################

cat("Rendering variants.rmd/html\n")

setwd("../generators/")
# Set the working directory.

Sys.setenv(RSTUDIO_PANDOC="C:/Program Files/RStudio/bin/pandoc/")
# Set the pandoc path

library(rmarkdown)
# Load markdown

output_dir = "../.."
# Set the output directory to be the docs directory, two up

render("../markdown/variants.rmd", output_dir = output_dir, params = list(output_dir = output_dir), quiet = TRUE)
# Knit the html page