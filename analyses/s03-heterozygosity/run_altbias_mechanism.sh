#!/bin/bash
# s03 mechanism dive: alt-read deficit vs divergence, from the s01 (SNP-only, identity-coord)
# BAMs. iter1 is single-shot to the reference (caller-independent -> one caller suffices);
# iter3 is measured per caller to show the correction (and its caller-dependence). Emits a
# tidy CSV consumed by the "mechanism" subsection of s03-results.rmd. Existing BAMs only.
set -uo pipefail
cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations
LIB=analyses/lib/altbias_mechanism.py
DATA=data/s01-baseline-iterative
OUT=analyses/s03-heterozygosity/results/s03-altbias-mechanism.csv
STRIDE=25
DIVS="0.00 0.02 0.04 0.06 0.08 0.10"

echo "divergence,iteration,caller,n_sites,n_ref,n_alt,alt_ref_ratio,mapq_ref,mapq_alt,pct_lowmapq_ref,pct_lowmapq_alt,pct_clip_ref,pct_clip_alt" > "$OUT"
run() { # div iter caller
  local D=$1 IT=$2 C=$3
  local B=$DATA/iter/$C/30X/19/$D/0.005/42/iter$IT/mm39_mapped.bam
  local T=$DATA/sim/30X/19/$D/0.005/42/mm39-truth.vcf.gz
  [ -f "$B" ] || { echo "  MISSING $B" >&2; return; }
  python "$LIB" "$B" "$T" 19 --stride $STRIDE --label "$D,$IT,$C" --csv >> "$OUT"
  echo "  done D=$D iter$IT $C" >&2
}
for D in $DIVS; do run "$D" 1 bcftools; done           # single-shot: caller-independent
for D in $DIVS; do for C in bcftools gatk freebayes; do run "$D" 3 "$C"; done; done
echo "WROTE $OUT ($(($(wc -l < "$OUT")-1)) rows)" >&2
