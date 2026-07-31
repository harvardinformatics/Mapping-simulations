#!/bin/bash
#SBATCH --job-name=caller_bench
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=2-00:00:00
#SBATCH --output=snakemake/slurm-logs/caller_bench-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gthomas@g.harvard.edu

# Orchestrator for the variant-caller benchmark (gatk vs bcftools vs freebayes).
# Runs snakemake, which submits the per-caller/per-condition jobs via the SLURM
# profile. Must run from the snakemake/ dir so the relative -s / --configfile /
# --profile paths resolve.

set -euo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims

cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations/snakemake

snakemake -p -s caller_benchmark.smk \
    --configfile ../simulation-configs/mm39-mason.yaml \
    --profile profiles/slurm_profile/ \
    --rerun-incomplete \
    --keep-going
