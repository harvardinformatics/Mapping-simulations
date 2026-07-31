#!/bin/bash
#SBATCH --job-name=caller_cmp
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --mem=16g
#SBATCH --time=6:00:00
#SBATCH --output=snakemake/slurm-logs/caller_cmp-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gthomas@g.harvard.edu

# One-off: does a PILEUP-based caller (bcftools) recover the divergence sites that
# GATK HaplotypeCaller (assembly-based) lost, at div=0.10 (SNP-only condition)?
# Scored by outcome class from the site-linkage analysis. --min-MQ 20 matches GATK.
# (Superseded by the full caller benchmark in caller_benchmark.smk /
# run_caller_benchmark.sh; kept for the single-condition recovery check.)

set -euo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims

PROJ=/n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations
cd "$PROJ"
BASE=data/mm39-mason-chr19
REF=$BASE/ref-regions/mm39-19.fa
MAPPED=$BASE/mapped-reads/30X/19/0.10/0.005/mm39_mapped.bam
TRUTH=$BASE/truth/19/0.10/0.005/mm39-19-0.10d-0.005h-truth.vcf.gz
OUT=$BASE/caller-comparison
mkdir -p "$OUT"

# pileup-based calls
bcftools mpileup -f "$REF" --min-MQ 20 -a FORMAT/AD,FORMAT/DP --threads 4 "$MAPPED" \
  | bcftools call -m -v --threads 4 -Oz -o "$OUT/bcftools.vcf.gz"
bcftools norm -f "$REF" -m- "$OUT/bcftools.vcf.gz" -Oz -o "$OUT/bcftools.norm.vcf.gz"
tabix -f -p vcf "$OUT/bcftools.norm.vcf.gz"

# recovery by site class (uses the corrected comparison script)
python scripts/compare_callers.py \
    "$OUT/bcftools.norm.vcf.gz" "$TRUTH" \
    summary-data/mm39-mason-chr19-site-linkage.csv 0.10 \
    "$OUT/caller-comparison-0.10.txt"
echo "DONE -> $OUT/caller-comparison-0.10.txt"
