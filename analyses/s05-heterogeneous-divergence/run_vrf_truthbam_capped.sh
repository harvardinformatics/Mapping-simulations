#!/bin/bash
#SBATCH --job-name=s05_vrf_tbam_cap
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=24g
#SBATCH --time=6:00:00
#SBATCH --output=../../snakemake/slurm-logs/s05_vrf_tbam_cap-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gthomas@g.harvard.edu

# Regenerate the mason -oa TRUTH BAM for the REALISTIC CAPPED config (L20000-S1.0-C0.20,
# mean D=0.03, local divergence capped at 0.20) so the mapping-vs-calling ledger can be run
# on the actual capped data (not the uncapped characterization config). Deterministic: same
# reference + truth VCF + seed reproduces identical reads.
set -eo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims
cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations

FIELD=L20000-S1.0-C0.20; DIV=0.03; HET=0.005; SEED=42
READLEN=150; FRAGMEAN=300; COV=30
D=data/s05-heterogeneous-divergence
REF=$D/ref-regions/mm39-19.fa
VCF=$D/sim/${COV}X/19/$DIV/$HET/$FIELD/$SEED/mm39-truth.vcf.gz
OUT=$D/varfate/$FIELD
mkdir -p "$OUT"

[ -s "$VCF" ] || { echo "ERROR: truth VCF missing: $VCF"; exit 1; }
if [ -s "$OUT/truth.oa.bam" ]; then echo "already have $OUT/truth.oa.bam"; exit 0; fi

LEN=$(cut -f2 "$REF".fai | head -1)
NFRAG=$(( COV * LEN / (2 * READLEN) ))
echo "### mason -oa truth BAM: field=$FIELD D=$DIV  NFRAG=$NFRAG ###"

TMP=$(mktemp -d)
mason_simulator -ir "$REF" -iv "$VCF" -n $NFRAG \
    --illumina-read-length $READLEN --fragment-mean-size $FRAGMEAN \
    --num-threads 8 --seed $SEED \
    -o "$TMP"/r1.fq.gz -or "$TMP"/r2.fq.gz -oa "$TMP"/truth.sam --version-check 0
samtools view -b "$TMP"/truth.sam -o "$OUT/truth.oa.bam"
rm -rf "$TMP"
echo "wrote $OUT/truth.oa.bam"
