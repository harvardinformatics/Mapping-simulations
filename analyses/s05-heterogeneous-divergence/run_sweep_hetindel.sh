#!/bin/bash
#SBATCH --job-name=s05_hetindel
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=3-00:00:00
#SBATCH --output=../../snakemake/slurm-logs/s05_hetindel-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gthomas@g.harvard.edu

# Full s05 heterogeneous-divergence-with-indels run: all 3 callers x 3 iterations for the
# isolation config (mean 0.06, sigma=1, indels on field, het uniform, uncapped). The
# bcftools validation chain already passed (liftover + field verified), so this reuses it
# and fills in gatk/freebayes. Homogeneous baseline for comparison = s02 D=0.06.
# Launch from the s05 directory:  sbatch run_sweep_hetindel.sh

set -euo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims

cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations/analyses/s05-heterogeneous-divergence

snakemake -s workflow-indels.smk --configfile config/config-hetindel.yaml --unlock || true

snakemake -p -s workflow-indels.smk --configfile config/config-hetindel.yaml \
    --profile ../../snakemake/profiles/slurm_profile/ \
    --rerun-incomplete --keep-going
