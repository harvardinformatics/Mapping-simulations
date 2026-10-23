#!/bin/bash
#SBATCH --job-name=s02_sweep
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=3-00:00:00
#SBATCH --output=../../snakemake/slurm-logs/s02_sweep-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gthomas@g.harvard.edu

# Orchestrator for the full s02 (indels + liftover) sweep. Runs as its own job so it
# survives a disconnected session; snakemake submits per-rule jobs via the SLURM profile.
# NO --forcerun (forces every instance of a rule and balloons the DAG).
# Launch from this directory:  sbatch run_sweep.sh

set -euo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims

cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations/analyses/s02-indels

snakemake -s workflow.smk --configfile config/config.yaml --unlock || true

snakemake -p -s workflow.smk --configfile config/config.yaml \
    --profile ../../snakemake/profiles/slurm_profile/ \
    --rerun-incomplete \
    --keep-going
# NOTE: default (code) rerun-triggers. After editing only a rule's RESOURCES (e.g. Picard
# heap), this re-runs that rule's outputs but REUSES upstream map/call. Do NOT use
# --rerun-triggers mtime here: retry rounds updated consensus mtimes, so mtime spuriously
# re-runs the expensive gatk map/call.
