#!/bin/bash
#SBATCH --job-name=s05_realistic
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=3-00:00:00
#SBATCH --output=../../snakemake/slurm-logs/s05_realistic-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gthomas@g.harvard.edu

# Orchestrator for the s05 REALISTIC scenario (config-realistic.yaml): mean D=0.03,
# sigma=1, local-D capped at 0.20, heterozygosity on the same field as divergence,
# plus a mean-0.03 homogeneous control. Runs as its own SLURM job so it survives a
# disconnected session. Launch from the s05 directory:  sbatch run_realistic.sh

set -euo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims

cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations/analyses/s05-heterogeneous-divergence

snakemake -s workflow.smk --configfile config/config-realistic.yaml --unlock || true

snakemake -p -s workflow.smk --configfile config/config-realistic.yaml \
    --profile ../../snakemake/profiles/slurm_profile/ \
    --rerun-incomplete \
    --keep-going
