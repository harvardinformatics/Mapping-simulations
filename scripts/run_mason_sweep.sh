#!/bin/bash
#SBATCH --job-name=mason_sweep
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=3-00:00:00
#SBATCH --output=snakemake/slurm-logs/mason_sweep-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gthomas@g.harvard.edu

# Orchestrator for the v2 (mason2) reference-bias sweep. This long, low-resource
# job runs snakemake, which submits the per-rule cluster jobs via the SLURM
# profile. Running the orchestrator as its own job lets the sweep survive a
# disconnected login session.
#
# --keep-going : a single failed condition does not abort the others, so we can
#                inspect partial results. (The combine step still needs all
#                conditions; re-run after any fixes to finish it.)

set -euo pipefail

source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims

cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations/snakemake

snakemake -p -s simulate_and_map.smk \
    --configfile ../simulation-configs/mm39-mason.yaml \
    --profile profiles/slurm_profile/ \
    --rerun-incomplete \
    --keep-going
