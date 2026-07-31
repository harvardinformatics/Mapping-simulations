#!/bin/bash
#SBATCH --job-name=s01_hetrec
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=12:00:00
#SBATCH --output=../../snakemake/slurm-logs/s01_hetrec-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gthomas@g.harvard.edu

# Orchestrator for the s01 heterozygous-site-recovery add-on. Scores each EXISTING
# per-iteration caller VCF (calls.norm) against truth with compare_vcfs and reads the
# het_* columns. Cheap (no simulation); only reads outputs of run_sweep.sh.
#
# Launch from this directory:  sbatch run_het_recovery.sh
# NO --forcerun (that forces every instance of a rule and balloons the DAG).

set -euo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims

cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations/analyses/s01-baseline-iterative

snakemake -s het_recovery.smk --configfile config/config.yaml --unlock || true

snakemake -p -s het_recovery.smk --configfile config/config.yaml \
    --profile ../../snakemake/profiles/slurm_profile/ \
    --rerun-incomplete \
    --keep-going
