#!/bin/bash
#SBATCH --job-name=s07_svs
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=3-00:00:00
#SBATCH --output=../../snakemake/slurm-logs/s07_svs-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gthomas@g.harvard.edu

# Full s07 SV-confounder run: all 3 callers x 3 iterations (SVs + SNP/indel at D=0.02).
# The bcftools validation chain already passed (SV truth flows through the liftover), so
# this reuses it and fills in gatk/freebayes.
# Launch from the s07 directory:  sbatch run_sweep_sv.sh

set -euo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims

cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations/analyses/s07-svs

snakemake -s workflow.smk --configfile config/config-sv.yaml --unlock || true

snakemake -p -s workflow.smk --configfile config/config-sv.yaml \
    --profile ../../snakemake/profiles/slurm_profile/ \
    --rerun-incomplete --keep-going
