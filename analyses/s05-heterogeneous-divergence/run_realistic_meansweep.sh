#!/bin/bash
#SBATCH --job-name=s05_real_mean
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=3-00:00:00
#SBATCH --output=../../snakemake/slurm-logs/s05_real_mean-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gthomas@g.harvard.edu

# s05 realistic scenario across MEAN divergence (SNP-only). Runs the realistic field config
# (sigma=1, ell=20kb, cap 0.20, het-on-field) at several genome-wide mean divergences, each
# written to its own data branch (DIV is in the path). Sequential per mean; mean 0.03 is
# already done (config-realistic) and reused. After all means, combine every per-condition
# summary across means into one table for the "reference bias vs mean divergence, realistic"
# curve (homogeneous baseline for comparison = s01).
# Launch from the s05 directory:  sbatch run_realistic_meansweep.sh

set -euo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims

cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations/analyses/s05-heterogeneous-divergence

CFG=config/config-realistic-meansweep.yaml
PROFILE=../../snakemake/profiles/slurm_profile/
DATA=/n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations/data/s05-heterogeneous-divergence
RES=/n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations/analyses/s05-heterogeneous-divergence/results
LIB=/n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations/analyses/lib

# each mean written to its own DIV path branch; 0.03 already exists and is skipped
for D in 0.01 0.02 0.03 0.04 0.06; do
    echo "### realistic mean sweep: div_target=$D ###"
    snakemake -s workflow.smk --configfile $CFG --config div_target=$D --unlock || true
    snakemake -p -s workflow.smk --configfile $CFG --config div_target=$D \
        --profile $PROFILE --rerun-incomplete --keep-going
done

# combine every realistic-config per-iteration summary across all means into one table
echo "### combining across means ###"
mapfile -t SUMS < <(find "$DATA/iter" -path "*L20000-S1.0-C0.20*/iter*/mm39-vcf-summary.csv" | sort)
python "$LIB/combine_summaries.py" --inputs "${SUMS[@]}" \
    --out "$RES/s05-realistic-meansweep-iterative-summary.csv" --kind vcf
echo "### DONE: $RES/s05-realistic-meansweep-iterative-summary.csv (${#SUMS[@]} rows) ###"
