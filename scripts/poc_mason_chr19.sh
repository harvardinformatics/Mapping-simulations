#!/bin/bash
#SBATCH --job-name=mason_poc_chr19
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=24g
#SBATCH --time=8:00:00
#SBATCH --output=snakemake/slurm-logs/mason_poc_chr19-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gthomas@g.harvard.edu

# Phase 1 proof-of-concept for the mason2-based reference-bias pipeline (v2).
# Single chromosome (chr19), single divergence level (0.04), 30X, diploid (het via -n 2).
# Validates: (1) realized SNP rate ~= target, (2) mason truth BAM is in REFERENCE
# coordinates so no CrossMap/liftover is needed, (3) the constant "missing reads"
# artifact from the v1 liftover is gone.
# NOTE: standalone POC; snakemake-ification is Phase 2. GATK/variant-call comparison
# is deferred (it is the slow part and not needed to confirm the core mechanism).

set -euo pipefail

PROJ=/n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations
REGION_FA=$PROJ/data/reference-genomes/mm39/mm39-19.fa   # extracted chr19
OUTDIR=$PROJ/data/mm39-mason-poc
THREADS=8
SEED=42
COV=30
READLEN=150
FRAGMEAN=300
SNP_RATE=0.04          # target divergence
SMALL_INDEL_RATE=0.0004
DIV_LABEL=0.04

source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims

mkdir -p "$OUTDIR"
cd "$OUTDIR"

# ---- reference for this run (copy chr19 in, index) -------------------------
cp -f "$REGION_FA" ref.fa
samtools faidx ref.fa
bwa index ref.fa

LEN=$(cut -f2 ref.fa.fai | head -1)
NFRAG=$(( COV * LEN / (2 * READLEN) ))
echo "# region length = $LEN ; num fragments for ${COV}X = $NFRAG"

# ---- 1. simulate variants (divergence + het) -> truth VCF (ref coords) -----
# -n 2 => diploid (heterozygosity); SVs disabled for the POC.
mason_variator \
    -ir ref.fa \
    -ov truth.vcf \
    -of truth_haplotypes.fa \
    --out-breakpoints truth_breakpoints.tsv \
    -n 2 \
    --snp-rate "$SNP_RATE" \
    --small-indel-rate "$SMALL_INDEL_RATE" \
    --sv-indel-rate 0 --sv-inversion-rate 0 --sv-translocation-rate 0 --sv-duplication-rate 0 \
    -s "$SEED" \
    --version-check 0

# ---- 2. simulate reads + TRUTH alignment BAM (against ref coords) ----------
mason_simulator \
    -ir ref.fa \
    -iv truth.vcf \
    -n "$NFRAG" \
    --illumina-read-length "$READLEN" \
    --fragment-mean-size "$FRAGMEAN" \
    --num-threads "$THREADS" \
    --seed "$SEED" \
    -o reads_R1.fq.gz \
    -or reads_R2.fq.gz \
    -oa truth.unsorted.bam \
    --version-check 0

samtools sort -@ "$THREADS" -o truth.bam truth.unsorted.bam
samtools index truth.bam
rm -f truth.unsorted.bam

# ---- 3. map reads back to the reference with bwa ---------------------------
bwa mem -t "$THREADS" ref.fa reads_R1.fq.gz reads_R2.fq.gz \
    | samtools sort -@ "$THREADS" -o mapped.bam -
samtools index mapped.bam

# ---- 4. validation gates ---------------------------------------------------
{
  echo "############ mason POC validation (chr19, div=$DIV_LABEL, ${COV}X) ############"
  echo "## region length: $LEN"
  echo
  echo "== GATE 1: realized variant rates vs target (target SNP rate = $SNP_RATE) =="
  NSNP=$(bcftools view -H -v snps truth.vcf | wc -l)
  NIND=$(bcftools view -H -v indels truth.vcf | wc -l)
  echo "truth SNPs:   $NSNP  => realized SNP rate = $(echo "scale=5; $NSNP/$LEN" | bc)"
  echo "truth indels: $NIND  => realized indel rate = $(echo "scale=6; $NIND/$LEN" | bc)"
  echo
  echo "== GATE 3: no 'missing reads' artifact (truth vs mapped read sets) =="
  samtools view -f 0x40 truth.bam  | cut -f1 | sort -u > truth_r1.qnames
  samtools view -f 0x40 mapped.bam | cut -f1 | sort -u > mapped_r1.qnames
  NT=$(wc -l < truth_r1.qnames)
  NM=$(wc -l < mapped_r1.qnames)
  NMISS=$(comm -23 truth_r1.qnames mapped_r1.qnames | wc -l)
  echo "R1 fragments in truth BAM:  $NT"
  echo "R1 fragments in mapped BAM: $NM"
  echo "in truth but MISSING from mapped: $NMISS  (v1 had a constant ~840k-1.4M here)"
  echo
  echo "== flagstat: truth.bam =="; samtools flagstat truth.bam
  echo "== flagstat: mapped.bam =="; samtools flagstat mapped.bam
  echo "############ done ############"
} > poc_validation.txt 2>&1

echo "POC complete. See $OUTDIR/poc_validation.txt"
