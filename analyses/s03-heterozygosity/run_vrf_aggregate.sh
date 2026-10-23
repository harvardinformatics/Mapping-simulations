#!/bin/bash
#SBATCH --job-name=s03_vrf_agg
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=40g
#SBATCH --time=10:00:00
#SBATCH --output=../../snakemake/slurm-logs/s03_vrf_agg-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gthomas@g.harvard.edu

# Full-chromosome per-variant read-fate ledger for s03 (HETEROZYGOUS sites), aggregate mode.
# Runs on the s01 (SNP-only) data so coordinates are identity across iterations (no chain needed).
# At a het site mason puts the alt on ONE phased haplotype, so only ~half the covering reads truly
# carry the alt (n_true_alt, via the oH tag); f_alt = how much of that alt support aligns at p (the
# het-dropout mechanism). "Recovered" = the caller called the site correctly AS A HET (GT 0/1) in
# its per-iteration VCF (mm39_calls.norm.vcf.gz) -- so we pass that, not the hom recovered set.
# One condition per job:  sbatch run_vrf_aggregate.sh <CALLER> <DIV>
set -eo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims
cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations

C=$1; D=$2; SEED=42
TB=data/varfate/$D/truth.oa.bam                              # cached mason -oa truth BAM (has oH)
[ -s "$TB" ] || { echo "ERROR: truth BAM missing: $TB"; exit 1; }
BASE=data/s01-baseline-iterative
NTRUTH=$BASE/sim/30X/19/$D/0.005/$SEED/mm39-truth.norm.vcf.gz
OUT=analyses/s03-heterozygosity/results/vrf-agg
mkdir -p "$OUT"

# per-iteration caller VCF is the "called as het" set; identity coordinates => no chain
IT=""
for k in 1 2 3; do
    d=$BASE/iter/$C/30X/19/$D/0.005/$SEED/iter$k
    IT="$IT --iter $k:$d/mm39_mapped.bam:$d/mm39_calls.norm.vcf.gz"
done

echo "### s03 aggregate HET read-fate: caller=$C div=$D (full chromosome, s01 data) ###"
python analyses/lib/variant_read_fate.py --het --truth-bam "$TB" --truth-vcf "$NTRUTH" \
    --region 19 --divergence "$D" --caller "$C" --sample 0 --aggregate $IT \
    --out "$OUT/vrf-agg-$C-$D"
echo "### DONE caller=$C div=$D ###"
