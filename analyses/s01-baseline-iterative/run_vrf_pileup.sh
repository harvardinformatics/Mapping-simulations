#!/bin/bash
#SBATCH --job-name=s01_vrf_pile
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=32g
#SBATCH --time=8:00:00
#SBATCH --output=../../snakemake/slurm-logs/s01_vrf_pile-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gthomas@g.harvard.edu

# s01 companion decomposition: per-SNP read fate WITH the position-centric pileup at p (--pileup).
# Purpose -- for SNPs that flip not-called -> called while staying in the top f_alt bin, what actually
# changes? The truth-centric ledger (f_alt, ex_/sh_ cells) is blind to reads mismapped INTO p, to
# depth, and to MAPQ; --pileup adds the caller's-eye pileup (depth_p, ref_p, alt_p, oth_p, mapq_p,
# mapq20_p) so the decomposition can attribute the flip to ref/contaminant reads, depth, or MAPQ.
# Sampled (per-SNP mode, ~30k SNPs) because the pileup pass is per-position. Identity coords (s01).
# One condition per job:  sbatch run_vrf_pileup.sh <CALLER> <DIV>
set -eo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims
cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations

C=$1; D=$2; SEED=42; SAMPLE=30000
TB=data/varfate/$D/truth.oa.bam                              # cached mason -oa truth BAM (has oH)
[ -s "$TB" ] || { echo "ERROR: truth BAM missing: $TB"; exit 1; }
BASE=data/s01-baseline-iterative
NTRUTH=$BASE/sim/30X/19/$D/0.005/$SEED/mm39-truth.norm.vcf.gz
MTRUTH=$BASE/sim/30X/19/$D/0.005/$SEED/mm39-truth.vcf.gz      # mason-format (true read spans)
REF=$BASE/ref-regions/mm39-19.fa
OUT=analyses/s01-baseline-iterative/results/vrf-pileup
mkdir -p "$OUT"

# same iter specs as the aggregate run; recovered set = consensus-diff SNPs (mm39_recovered.vcf.gz)
IT=""
for k in 1 2 3; do
    d=$BASE/iter/$C/30X/19/$D/0.005/$SEED/iter$k
    IT="$IT --iter $k:$d/mm39_mapped.bam:$d/mm39_recovered.vcf.gz"
done

echo "### s01 per-SNP read-fate + pileup: caller=$C div=$D (sample $SAMPLE) ###"
python analyses/lib/variant_read_fate.py --truth-bam "$TB" --truth-vcf "$NTRUTH" \
    --mason-vcf "$MTRUTH" --reference "$REF" \
    --region 19 --divergence "$D" --caller "$C" --sample $SAMPLE --pileup $IT \
    --out "$OUT/vrf-pileup-$C-$D.csv"
echo "### DONE caller=$C div=$D ###"
