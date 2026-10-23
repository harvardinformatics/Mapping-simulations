#!/bin/bash
#SBATCH --job-name=s03_ab
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=1-00:00:00
#SBATCH --output=../../snakemake/slurm-logs/s03_ab-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gthomas@g.harvard.edu

# Orchestrator for the s03 allele-balance analysis (Phase 1). Runs allele_balance.py over
# the EXISTING s01 and s02 BAMs at true het sites (no new simulation). Runs as its own job
# so it survives a disconnect. Launch from this directory:  sbatch run_allele_balance.sh

set -euo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims

cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations/analyses/s03-heterozygosity

for cfg in config/config-s01.yaml config/config-s02.yaml; do
    echo "=== $cfg ==="
    snakemake -s het_allele_balance.smk --configfile "$cfg" --unlock || true
    snakemake -p -s het_allele_balance.smk --configfile "$cfg" \
        --profile ../../snakemake/profiles/slurm_profile/ \
        --rerun-incomplete --keep-going
done
