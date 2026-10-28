#!/bin/bash
#SBATCH --job-name=s01_cluster
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=24g
#SBATCH --time=2:00:00
#SBATCH --output=../../snakemake/slurm-logs/s01_cluster-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gthomas@g.harvard.edu

# s01 recall of divergent SNPs by local SNP clustering: for every divergent truth SNP, the number of
# OTHER divergent SNPs within +/-W bp, and whether each caller recovered it at each iteration
# (analyses/lib/cluster_recall.py). Tests whether GATK's loss concentrates in SNP clusters (its
# k-mer re-assembly) while bcftools/freebayes are insensitive. Reads only VCFs already on disk.
# Every condition is checked against s01-base-iterative-summary.csv (sites = n_truth_div,
# recovered = div_tp) and the run stops on any mismatch.
# Usage (from this directory):  sbatch run_cluster_recall.sh    (or: bash run_cluster_recall.sh)
set -eo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims
cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations

# --- settings ----------------------------------------------------------------------------------
DIVS=0.001,0.005,0.02,0.04,0.06,0.08,0.10   # every s01 divergence with divergent SNPs (D = 0 has none)
SEEDS=42,43,44
WINDOWS=10,25                               # half-windows (bp); 25 = GATK's largest default k-mer
BASE=data/s01-baseline-iterative
OUT=analyses/s01-baseline-iterative/results/cluster-recall/s01-cluster-recall-w${WINDOWS//,/-w}.csv

python analyses/lib/cluster_recall.py \
    --truth     "$BASE/sim/30X/19/{D}/0.005/{seed}/mm39-truth.norm.vcf.gz" \
    --recovered "$BASE/iter/{caller}/30X/19/{D}/0.005/{seed}/iter{k}/mm39_recovered.vcf.gz" \
    --divs $DIVS --seeds $SEEDS --iterations 1,2,3 --windows $WINDOWS --threads 8 \
    --check-summary analyses/s01-baseline-iterative/results/s01-base-iterative-summary.csv \
    --out $OUT
