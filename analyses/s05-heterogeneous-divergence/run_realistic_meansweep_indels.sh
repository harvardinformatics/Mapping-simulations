#!/bin/bash
#SBATCH --job-name=s05_real_mean_indel
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=3-00:00:00
#SBATCH --output=../../snakemake/slurm-logs/s05_real_mean_indel-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gthomas@g.harvard.edu

# s05 REALISTIC scenario across MEAN divergence, WITH INDELS (the indel counterpart of
# run_realistic_meansweep.sh). Realistic field (sigma=1, ell=20kb, cap 0.20, het-on-field)
# swept over mean divergence 0.01-0.06 with indels on the field. This lightweight coordinator
# runs snakemake, which dispatches the per-condition jobs to SLURM via the cluster profile
# (GATK HaplotypeCaller is the slow arm, ~many hours per high-divergence call -> multi-day).
# Produces results/s05-realistic-meansweep-indels-recovery-summary.csv (combine rule).
# Launch from the s05 directory:  sbatch run_realistic_meansweep_indels.sh

set -euo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims

cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations/analyses/s05-heterogeneous-divergence

CFG=config/config-realistic-meansweep-indels.yaml
PROFILE=../../snakemake/profiles/slurm_profile/

snakemake -s workflow-indels.smk --configfile "$CFG" --unlock || true

snakemake -p -s workflow-indels.smk --configfile "$CFG" \
    --profile "$PROFILE" \
    --rerun-incomplete --keep-going

echo "### DONE: results/s05-realistic-meansweep-indels-recovery-summary.csv ###"
