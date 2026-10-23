#!/bin/bash
#SBATCH --job-name=s05_validate
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=3-00:00:00
#SBATCH --output=../../snakemake/slurm-logs/s05_validate-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gthomas@g.harvard.edu

# Validation gate step 2: run ONE full field configuration (L20000-S1.0) end-to-end --
# all three callers, all iterations -- before committing the whole sweep. This confirms
# the field truth VCF flows through mapping/calling/consensus/scoring and that recovery
# varies with LOCAL divergence. Orchestrator runs as its own SLURM job (submits per-rule
# jobs via the profile) so it survives a disconnected session.
#
# Launch from the s05 directory:  sbatch run_validate_chain.sh

set -euo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims

cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations/analyses/s05-heterogeneous-divergence

DATA=/n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations/data/s05-heterogeneous-divergence
FIELD=L20000-S1.0
# the three iter3 summaries pull the entire chain (iter1->2->3) for each caller
TARGETS=""
for c in gatk bcftools freebayes; do
    TARGETS="$TARGETS $DATA/iter/$c/30X/19/0.06/0.005/$FIELD/42/iter3/mm39-vcf-summary.csv"
done

snakemake -s workflow.smk --configfile config/config.yaml --unlock || true

snakemake -p -s workflow.smk --configfile config/config.yaml \
    --profile ../../snakemake/profiles/slurm_profile/ \
    --rerun-incomplete --keep-going \
    $TARGETS
