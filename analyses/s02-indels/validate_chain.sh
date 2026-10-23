#!/bin/bash
#SBATCH --job-name=s02_valchain
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=1-00:00:00
#SBATCH --output=../../snakemake/slurm-logs/s02_valchain-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gthomas@g.harvard.edu

# GO/NO-GO validation of ONE full recursive chain (iter1->2->3) before the s02 sweep.
# Tests the recursive rebuild-from-reference liftover end to end for a single condition
# (bcftools, D=0.06, seed 42 -- fast, and the liftover is caller-independent). Confirms the
# accumulation + chains don't drift and that SNP+indel recovery scores sanely each round.
# Launch from this directory:  sbatch validate_chain.sh

set -euo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims

DIR=/n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations/analyses/s02-indels
D=/n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations/data/s02-indels
cd $DIR

COND=bcftools/30X/19/0.06/0.005/42
T=""
for k in 1 2 3; do
  T="$T $D/iter/$COND/iter$k/mm39-recovery-vcf-summary.csv"
  T="$T $D/iter/$COND/iter$k/mm39-het-vcf-summary.csv"
done

snakemake -s workflow.smk --configfile config/config.yaml --unlock || true
snakemake -p -s workflow.smk --configfile config/config.yaml \
    --profile ../../snakemake/profiles/slurm_profile/ \
    --rerun-incomplete --keep-going $T

echo
echo "############## S02 SINGLE-CHAIN VALIDATION ($COND) ##############"
echo "### divergence RECOVERY (accumulated recovered vs truth), by iteration:"
echo "### watch div_tp (SNP) and indel_div_tp rise across iters; consensus length should"
echo "### track total recovered indel bases; nothing should crash or lose the truth."
for k in 1 2 3; do
  f=$D/iter/$COND/iter$k/mm39-recovery-vcf-summary.csv
  [ -f "$f" ] && { echo "--- iter$k recovery ---"; awk -F, 'NR==1{for(i=1;i<=NF;i++)h[$i]=i} NR==2{
     printf "  SNP:   n_truth_div=%s div_tp=%s (%.1f%%)  ref_or_absent=%s\n", $h["n_truth_div"],$h["div_tp"],100*$h["div_tp"]/$h["n_truth_div"],$h["div_ref_or_absent"];
     printf "  INDEL: n_truth=%s tp=%s (%.1f%%)  ref_or_absent=%s\n", $h["n_truth_indel_div"],$h["indel_div_tp"],($h["n_truth_indel_div"]>0?100*$h["indel_div_tp"]/$h["n_truth_indel_div"]:0),$h["indel_div_ref_or_absent"] }' "$f"; }
  fh=$D/iter/$COND/iter$k/mm39-het-vcf-summary.csv
  [ -f "$fh" ] && awk -F, 'NR==1{for(i=1;i<=NF;i++)h[$i]=i} NR==2{printf "  HET SNP: n_truth_het=%s het_tp=%s (%.1f%%)\n", $h["n_truth_het"],$h["het_tp"],($h["n_truth_het"]>0?100*$h["het_tp"]/$h["n_truth_het"]:0)}' "$fh"
done
echo "### consensus lengths (ref=$(cut -f2 $D/ref-regions/mm39-19.fa.fai 2>/dev/null|head -1)):"
for k in 1 2; do
  c=$D/iter/$COND/iter$k/mm39-consensus.fa.fai
  [ -f "$c" ] && echo "  iter$k cons length: $(cut -f2 $c|head -1)"
done
echo "############## END ####"
