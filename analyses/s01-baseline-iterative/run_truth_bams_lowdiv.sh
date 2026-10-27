#!/bin/bash
#SBATCH --job-name=s01_truthbam_low
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=24g
#SBATCH --time=8:00:00
#SBATCH --output=../../snakemake/slurm-logs/s01_truthbam_low-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gthomas@g.harvard.edu

# Regenerate the mason -oa truth BAMs (with oH tags) for the near-0 / low divergences so the aggregate
# read-fate (G plots) can be extended to D=0.02 and 0.005. Same invocation as read_classification's
# truth_bam rule (deterministic: same seed -> identical reads), written to data/varfate/<D>/ where the
# aggregate driver looks. seed 42 (the read-fate seed).
# Usage:  cd analyses/s01-baseline-iterative && sbatch run_truth_bams_lowdiv.sh [D ...]
#         (default D: 0.02 0.005, the original run; 2026-10-26 also run with 0.00 0.001 for the
#         corrected truth-site scan, run_effective_truth.sh)
set -eo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims
cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations

REF=data/s01-baseline-iterative/ref-regions/mm39-19.fa
READLEN=150; FRAGMEAN=300; COV=30; SEED=42
LEN=$(cut -f2 "$REF".fai | head -1); NFRAG=$(( COV * LEN / (2 * READLEN) ))

DIVS="${*:-0.02 0.005}"                 # divergences from the command line, else the original pair
for D in $DIVS; do
    VCF=data/s01-baseline-iterative/sim/30X/19/$D/0.005/$SEED/mm39-truth.vcf.gz
    [ -s "$VCF" ] || { echo "ERROR: missing truth VCF $VCF"; exit 1; }
    OUT=data/varfate/$D
    mkdir -p "$OUT"
    if [ -s "$OUT/truth.oa.bam" ]; then echo "already have $OUT/truth.oa.bam"; continue; fi
    echo "### mason -oa truth BAM: D=$D ###"
    TMP=$(mktemp -d)
    mason_simulator -ir "$REF" -iv "$VCF" -n $NFRAG \
        --illumina-read-length $READLEN --fragment-mean-size $FRAGMEAN \
        --num-threads 8 --seed $SEED \
        -o "$TMP"/r1.fq.gz -or "$TMP"/r2.fq.gz -oa "$TMP"/truth.sam --version-check 0
    samtools view -b "$TMP"/truth.sam -o "$OUT/truth.oa.bam"
    rm -rf "$TMP"
    echo "wrote $OUT/truth.oa.bam"
done
