#!/bin/bash
#SBATCH --job-name=s05_vrf_tbam
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=24g
#SBATCH --time=6:00:00
#SBATCH --output=../../snakemake/slurm-logs/s05_vrf_tbam-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gthomas@g.harvard.edu

# Regenerate the mason -oa TRUTH BAM (with oH haplotype tags) for the s05 mapping-vs-calling
# decomposition. s05 did not keep the truth alignment, but mason is deterministic: the same
# reference + truth VCF + seed reproduces the identical reads and their true positions. This
# is the truth-centric input variant_read_fate.py needs to know which reads TRULY cover each
# divergent site. Built for the reference characterization config (L20000-S1.0, mean D=0.06),
# whose per-site local divergence spans the full 0->~0.6 range where the caller walls appear.
set -eo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims
cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations

FIELD=L20000-S1.0; DIV=0.06; HET=0.005; SEED=42
READLEN=150; FRAGMEAN=300; COV=30
D=data/s05-heterogeneous-divergence
REF=$D/ref-regions/mm39-19.fa
VCF=$D/sim/${COV}X/19/$DIV/$HET/$FIELD/$SEED/mm39-truth.vcf.gz
OUT=$D/varfate/$FIELD
mkdir -p "$OUT"

[ -s "$REF" ]  || { echo "ERROR: reference missing: $REF"; exit 1; }
[ -s "$VCF" ]  || { echo "ERROR: truth VCF missing: $VCF"; exit 1; }
if [ -s "$OUT/truth.oa.bam" ]; then echo "already have $OUT/truth.oa.bam"; exit 0; fi

# same NFRAG formula the workflow's simulate_reads rule uses -> identical read set
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
