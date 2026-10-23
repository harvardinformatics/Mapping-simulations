#!/bin/bash
#SBATCH --job-name=s07_sv_val
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=2-00:00:00
#SBATCH --output=../../snakemake/slurm-logs/s07_sv_val-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gthomas@g.harvard.edu

# Validation gate for s07 (SVs as confounder): one recursive chain (bcftools iter1->3)
# end-to-end, confirming the SV-bearing truth survives normalize_truth + the s02 liftover
# (SVs are never recovered, so they must not break the small-variant chain). bcftools is
# fast and exercises the whole path. Launch from the s07 directory: sbatch run_validate_sv.sh

set -euo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims

cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations/analyses/s07-svs

DATA=/n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations/data/s07-svs
D=$DATA/iter/bcftools/30X/19/0.02/0.005/42/iter3
TARGETS="$D/mm39-recovery-vcf-summary.csv $D/mm39-het-vcf-summary.csv"

snakemake -s workflow.smk --configfile config/config-sv.yaml --unlock || true

snakemake -p -s workflow.smk --configfile config/config-sv.yaml \
    --profile ../../snakemake/profiles/slurm_profile/ \
    --rerun-incomplete --keep-going \
    $TARGETS
