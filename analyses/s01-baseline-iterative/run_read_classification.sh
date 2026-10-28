#!/bin/bash
#SBATCH --job-name=s01_readclass
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=1-00:00:00
#SBATCH --output=../../snakemake/slurm-logs/s01_readclass-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gthomas@g.harvard.edu

# Orchestrator for the s01 read-classification-by-iteration add-on. Regenerates the
# mason -oa truth BAM (deterministic, same seed => identical reads) for each div x seed,
# then classifies every EXISTING per-iteration mapped BAM against it with compare_bams.
# Reads the mapped BAMs produced by run_sweep.sh; does not rebuild them.
#
# Runs as its own job so it survives a disconnected session. Launch from this directory:
#   sbatch run_read_classification.sh
#
# NO --forcerun (that forces every instance of a rule and balloons the DAG).
# Extra arguments go to snakemake, e.g. `sbatch run_read_classification.sh --forcerun classify_reads`
# (2026-10-28: reclassify all 216 after the anchored placement start in truth_haplotypes.py, which
# snakemake cannot see -- it is imported by compare_bams.py, not a rule input).

set -euo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims

cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations/analyses/s01-baseline-iterative

# clear any stale lock from a previously-interrupted run, then run
snakemake -s read_classification.smk --configfile config/config.yaml --unlock || true

snakemake -p -s read_classification.smk --configfile config/config.yaml \
    --profile ../../snakemake/profiles/slurm_profile/ \
    --rerun-incomplete \
    --keep-going "$@"
