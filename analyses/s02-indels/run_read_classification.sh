#!/bin/bash
#SBATCH --job-name=s02_readclass
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=1-00:00:00
#SBATCH --output=../../snakemake/slurm-logs/s02_readclass-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gthomas@g.harvard.edu

# Orchestrator for the s02 read-classification-by-iteration add-on. Classifies each
# EXISTING per-iteration mapped BAM against the mason -oa truth, point-lifting iter>1 BAMs
# (consensus coords) back to reference via the cons->ref chains. Consumes outputs of
# run_sweep.sh; run AFTER it. Launch from this directory:  sbatch run_read_classification.sh
# Extra arguments go to snakemake, e.g. `sbatch run_read_classification.sh --forcerun classify_reads`
# (2026-10-27: reclassify all 72 after the compare_bams.py true-span fix; one table was made
# outside snakemake, so its code-change trigger cannot see the fix).

set -euo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims

cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations/analyses/s02-indels

snakemake -s read_classification.smk --configfile config/config.yaml --unlock || true

snakemake -p -s read_classification.smk --configfile config/config.yaml \
    --profile ../../snakemake/profiles/slurm_profile/ \
    --rerun-incomplete \
    --keep-going "$@"
