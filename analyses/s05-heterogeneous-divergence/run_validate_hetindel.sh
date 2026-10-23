#!/bin/bash
#SBATCH --job-name=s05_hetind_val
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=2-00:00:00
#SBATCH --output=../../snakemake/slurm-logs/s05_hetind_val-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gthomas@g.harvard.edu

# Validation gate for s05 heterogeneous-divergence WITH indels: run ONE recursive chain
# (bcftools, iter1->3) end-to-end before the full 3-caller run, confirming the field truth
# VCF flows through the s02 liftover machinery (Picard lift, accumulate, chainSwap consensus)
# with indels co-clustered. bcftools is fast and exercises the whole liftover path.
# Launch from the s05 directory:  sbatch run_validate_hetindel.sh

set -euo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims

cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations/analyses/s05-heterogeneous-divergence

DATA=/n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations/data/s05-heterogeneous-divergence/indels
D=$DATA/iter/bcftools/30X/19/0.06/0.005/42/iter3
# iter3 recovery + het summaries pull the whole bcftools chain (iter1->2->3)
TARGETS="$D/mm39-recovery-vcf-summary.csv $D/mm39-het-vcf-summary.csv"

snakemake -s workflow-indels.smk --configfile config/config-hetindel.yaml --unlock || true

snakemake -p -s workflow-indels.smk --configfile config/config-hetindel.yaml \
    --profile ../../snakemake/profiles/slurm_profile/ \
    --rerun-incomplete --keep-going \
    $TARGETS
