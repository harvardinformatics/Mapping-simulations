#!/bin/bash
#SBATCH --job-name=s02_ptlift
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=16g
#SBATCH --time=1:00:00
#SBATCH --output=../../snakemake/slurm-logs/s02_ptlift-%j.out
set -euo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims
BASE=/n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations
OUT=$BASE/data/s02-indels/validation; cd $OUT
echo "### point-lift read classification: map.cons.bam (cons coords) vs truth.oa.bam (ref), lifting via cons2ref chain ###"
python $BASE/analyses/lib/compare_bams.py \
    --truth-bam truth.oa.bam --mapped-bam map.cons.bam --chain cons2ref.chain \
    --mason-vcf truth.vcf.gz --reference $BASE/data/s01-baseline-iterative/ref-regions/mm39-19.fa \
    --region 19 --coverage 30 --divergence 0.06 --heterozygosity 0.005 --seed 42 \
    --caller bcftools --iteration 2 --out-prefix bamscore.pointlift 2> pointlift.log
echo "--- lift diagnostics ---"; grep -E 'loaded|not_in_truth|lift_fail' pointlift.log
echo "############## POINT-LIFT SUMMARY ##############"
column -s, -t bamscore.pointlift-bam-summary.csv
echo "############## END ####"
