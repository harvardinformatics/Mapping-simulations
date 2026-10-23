#!/bin/bash
#SBATCH --job-name=s05_sweep
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=3-00:00:00
#SBATCH --output=../../snakemake/slurm-logs/s05_sweep-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gthomas@g.harvard.edu

# Orchestrator for the s05 heterogeneous-divergence sweep. Runs snakemake, which submits
# the per-rule jobs via the SLURM profile. Runs as its own job so the sweep survives a
# disconnected session. Launch from the s05 directory:  sbatch run_sweep.sh
#
# Uses the default target (rule all -> combined summary + per-site local-divergence
# tables). --keep-going lets one failed configuration not abort the rest.
#
# NOTE: to run only ONE field configuration first (validation gate step 2), pass its
# target explicitly, e.g.:
#   snakemake -s workflow.smk --configfile config/config.yaml \
#       --profile ../../snakemake/profiles/slurm_profile/ \
#       <path to that config's iter.../iter3 -vcf-summary.csv>

set -euo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims

cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations/analyses/s05-heterogeneous-divergence

# clear any stale lock from a previously-interrupted run, then run
snakemake -s workflow.smk --configfile config/config.yaml --unlock || true

snakemake -p -s workflow.smk --configfile config/config.yaml \
    --profile ../../snakemake/profiles/slurm_profile/ \
    --rerun-incomplete \
    --keep-going
