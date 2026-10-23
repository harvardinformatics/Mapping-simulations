#!/bin/bash
#SBATCH --job-name=vrf_agg
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=40g
#SBATCH --time=10:00:00
#SBATCH --output=../../snakemake/slurm-logs/vrf_agg-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gthomas@g.harvard.edu

# Full-CHROMOSOME variant-read-fate ledger in AGGREGATE mode (figures G and H only).
# Scores every true hom-div SNP on chr19 (--sample 0) but writes only the two small summary
# tables (-gtraj.csv for G, -hcomp.csv for H) -- never the multi-GB per-SNP ledger.
# One condition (caller x divergence) per SLURM job:
#     sbatch run_vrf_aggregate.sh <CALLER> <DIV>
# The 12 conditions are submitted by submit_vrf_aggregate.sh. Truth BAMs (mason -oa, carrying the
# oH haplotype tag) must already exist under data/varfate/<DIV>/ -- this job does NOT regenerate
# them (so parallel jobs never race on the same file).
set -eo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims
cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations

C=$1; D=$2; SEED=42
TB=data/varfate/$D/truth.oa.bam
[ -s "$TB" ] || { echo "ERROR: truth BAM missing: $TB (generate it before submitting)"; exit 1; }

BASE=data/s01-baseline-iterative
NTRUTH=$BASE/sim/30X/19/$D/0.005/$SEED/mm39-truth.norm.vcf.gz
OUT=analyses/s01-baseline-iterative/results/vrf-agg
mkdir -p "$OUT"

IT=""
for k in 1 2 3; do
    d=$BASE/iter/$C/30X/19/$D/0.005/$SEED/iter$k
    IT="$IT --iter $k:$d/mm39_mapped.bam:$d/mm39_recovered.vcf.gz"
done

echo "### aggregate ledger: caller=$C div=$D (full chromosome, all SNPs) ###"
python analyses/lib/variant_read_fate.py --truth-bam "$TB" --truth-vcf "$NTRUTH" \
    --region 19 --divergence "$D" --caller "$C" --sample 0 --aggregate $IT \
    --out "$OUT/vrf-agg-$C-$D"
echo "### DONE caller=$C div=$D -> $OUT/vrf-agg-$C-$D-{gtraj,hcomp}.csv ###"
