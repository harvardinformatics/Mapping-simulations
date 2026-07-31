#!/bin/bash
#SBATCH --job-name=s01_sweep
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=3-00:00:00
#SBATCH --output=../../snakemake/slurm-logs/s01_sweep-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gthomas@g.harvard.edu

# Orchestrator for the s01 baseline-vs-iterative sweep. Runs snakemake, which submits
# the per-rule jobs via the SLURM profile. Runs as its own job so the sweep survives a
# disconnected session. Launch from the s01 directory:  sbatch run_sweep.sh
#
# Uses the plain default target (rule all -> the combined summary); NO --forcerun (that
# forces every instance of a rule and would balloon the DAG). --keep-going lets one
# failed condition not abort the rest.

set -euo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims

cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations/analyses/s01-baseline-iterative

# clear any stale lock from a previously-interrupted run, then run
snakemake -s workflow.smk --configfile config/config.yaml --unlock || true

snakemake -p -s workflow.smk --configfile config/config.yaml \
    --profile ../../snakemake/profiles/slurm_profile/ \
    --rerun-incomplete \
    --keep-going
