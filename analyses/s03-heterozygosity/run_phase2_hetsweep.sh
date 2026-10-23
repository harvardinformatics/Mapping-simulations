#!/bin/bash
#SBATCH --job-name=s03_p2
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=3-00:00:00
#SBATCH --output=../../snakemake/slurm-logs/s03_p2-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gthomas@g.harvard.edu

# s03 Phase 2 orchestrator — het-density confirmation sweep.
#
# For each elevated het rate H, (1) run the s01 workflow to simulate SNP-only reads and
# map/call/iterate them into new BAMs under data/s03-hetsweep, then (2) run allele_balance
# over those BAMs. H=0.005 is NOT re-run here (Phase 1's s01 already has it at D=0.04/0.10);
# the notebook combines that with these H=0.01/0.02 summaries. Runs as its own SLURM job so
# it survives a disconnect. Launch:  sbatch run_phase2_hetsweep.sh
#
# NOTE: default (code) rerun-triggers — do NOT add --rerun-triggers mtime (consensus mtimes
# get bumped by SLURM retries and would spuriously redo map/call; see s02 run_sweep.sh).

set -euo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims

cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations/analyses/s03-heterozygosity
PROFILE=../../snakemake/profiles/slurm_profile/
S01=../s01-baseline-iterative/workflow.smk

for H in 0.01 0.02; do
    TAG="hetsweep-$H"
    echo "############## H=$H : simulate + map + call + iterate ##############"
    snakemake -s "$S01" --configfile config/config-s01-hetsweep.yaml \
        --config het="$H" run_tag="$TAG" --unlock || true
    snakemake -p -s "$S01" --configfile config/config-s01-hetsweep.yaml \
        --config het="$H" run_tag="$TAG" \
        --profile "$PROFILE" --rerun-incomplete --keep-going

    echo "############## H=$H : allele balance (RAF at het sites) ##############"
    snakemake -s het_allele_balance.smk --configfile config/config-hetsweep.yaml \
        --config het="$H" run_tag="$TAG" --unlock || true
    snakemake -p -s het_allele_balance.smk --configfile config/config-hetsweep.yaml \
        --config het="$H" run_tag="$TAG" \
        --profile "$PROFILE" --rerun-incomplete --keep-going
done

echo "############## Phase 2 complete ##############"
