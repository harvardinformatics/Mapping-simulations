#!/bin/bash
# Convert a mason-format truth VCF (FORMAT='.', with the genotype value sitting in the
# sample column but no GT tag declared) into a valid GT-tagged VCF, so bcftools
# consensus / `bcftools view -g` can actually use the genotypes. mason's -iv reader
# tolerates the '.' FORMAT, but bcftools needs a declared GT.
#
# Usage: mason2gt.sh IN.vcf.gz OUT.vcf.gz [hom]
#   hom -> keep only homozygous-alt (1|1 or 1/1) records
set -euo pipefail
in=$1; out=$2; mode=${3:-all}
zcat "$in" | awk -v mode="$mode" 'BEGIN{OFS="\t"}
  /^##/ { print; next }
  /^#CHROM/ { print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"; print; next }
  {
    gt=$10
    if (mode=="hom" && gt!="1|1" && gt!="1/1") next
    $9="GT"; print
  }' | bgzip > "$out"
tabix -f -p vcf "$out"
