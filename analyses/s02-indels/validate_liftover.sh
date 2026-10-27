#!/bin/bash
#SBATCH --job-name=s02_liftval
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=24g
#SBATCH --time=4:00:00
#SBATCH --output=../../snakemake/slurm-logs/s02_liftval-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gthomas@g.harvard.edu

# s02 liftover-fidelity PROBE (go/no-go gate before building the indel sweep).
#
# Does CrossMap + a chainSwap'd bcftools-consensus chain lift indel-containing calls AND
# read alignments from consensus coords back to reference coords with ZERO spurious loss
# (the v1 "missing reads" artifact)?
#
# Design: bake the truth's homozygous SNP+indel variants into the reference to build a
# PERFECT indel consensus (simulates perfect iter1 recovery), map reads to it, call, keep
# the BAM in consensus coords, then lift both back to reference coords and check against
# reference-coordinate truth. Perfect consensus => reads map ~exactly, so a liftover
# artifact shows up starkly.
#
# GENOTYPE NOTE: make_truth_vcf writes mason-format truth (FORMAT='.', genotype in the
# sample col, no declared GT). Both bcftools consensus AND compare_vcfs need a real GT, so
# we tag it via scripts/mason2gt.sh (same idea as s01's normalize_truth awk $9="GT"). The
# real s02 workflow builds consensus from CALLER VCFs (which have GT) and normalizes truth
# exactly like s01 -- this probe just has to reproduce both.
#
# Direction: bcftools writes a ref->cons chain; CrossMap lifts target->query, so to lift
# consensus-coord data to reference we chainSwap it (target=cons). (Confirmed: the swapped
# cons2ref chain lifted 307,079 het SNPs to correct ref positions; the unswapped direction
# produced out-of-bounds coords, so we only use cons2ref now.)
#
# Reuses every valid intermediate from a prior run (sim, consensus, mapping, calls) so a
# rerun is fast. Self-contained sbatch (survives disconnect).

set -euo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims

BASE=/n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations
REF=$BASE/data/s01-baseline-iterative/ref-regions/mm39-19.fa
LIB=$BASE/analyses/lib
M2GT=$BASE/analyses/s02-indels/scripts/mason2gt.sh
OUT=$BASE/data/s02-indels/validation
mkdir -p $OUT; cd $OUT

D=0.06; DI=0.0072; H=0.005; SEED=42; COV=30; RL=150; FM=300   # Di = 0.12 * D

echo "### 1. truth VCF + reads + mason -oa alignment (reuse sim if present) ###"
if [ ! -s truth.vcf.gz ]; then
  python $LIB/make_truth_vcf.py --reference $REF --divergence $D --heterozygosity $H \
      --indel-divergence $DI --indel-heterozygosity 0 --min-indel-size 1 --max-indel-size 10 \
      --seed $SEED --out truth.vcf
  bgzip -f truth.vcf; tabix -f -p vcf truth.vcf.gz
fi
if [ ! -s truth.oa.bam ] || [ ! -s r1.fq.gz ] || [ ! -s r2.fq.gz ]; then
  LEN=$(cut -f2 $REF.fai | head -1); NFRAG=$(( COV * LEN / (2 * RL) ))
  mason_simulator -ir $REF -iv truth.vcf.gz -n $NFRAG --illumina-read-length $RL \
      --fragment-mean-size $FM --num-threads 8 --seed $SEED \
      -o r1.fq.gz -or r2.fq.gz -oa truth.oa.sam --version-check 0
  samtools view -b truth.oa.sam -o truth.oa.bam; rm -f truth.oa.sam
else
  echo "reusing existing reads + truth.oa.bam"
fi

echo "### 1b. GT-tagged, normalized truth for scoring (the bug fix) ###"
# plain 'bcftools norm' on mason-format truth DROPS the genotype; tag GT first (like s01).
bash $M2GT truth.vcf.gz truth.gt.vcf.gz
bcftools norm -f $REF -m- --atomize -Oz -o truth.norm.vcf.gz truth.gt.vcf.gz
tabix -f -p vcf truth.norm.vcf.gz
echo "truth GT check (want GT / 1|1): $(zcat truth.norm.vcf.gz | grep -v '^#' | head -1 | cut -f9,10)"

echo "### 2. PERFECT indel consensus from GT homozygous truth + chain (reuse) ###"
if [ ! -s cons.fa ] || [ ! -s cons2ref.chain ]; then
  bash $M2GT truth.vcf.gz truth.hom.gt.vcf.gz hom
  bcftools consensus -f $REF -c ref2cons.chain truth.hom.gt.vcf.gz > cons.fa
  samtools faidx cons.fa; bwa index cons.fa 2>/dev/null
  conda run -n chainswap chainSwap ref2cons.chain cons2ref.chain
fi
echo "ref length=$(cut -f2 $REF.fai|head -1)  cons length=$(cut -f2 cons.fa.fai|head -1)"

echo "### 3. map reads to the indel consensus (reuse) ###"
if [ ! -s map.cons.bam ]; then
  bwa mem -t 8 cons.fa r1.fq.gz r2.fq.gz 2>/dev/null | samtools sort -@ 8 -o map.cons.bam -
  samtools index map.cons.bam
fi

echo "### 4. call on the consensus (reuse) ###"
if [ ! -s calls.cons.norm.vcf.gz ]; then
  bcftools mpileup -f cons.fa -q 20 map.cons.bam 2>/dev/null | bcftools call -mv -Oz -o calls.cons.vcf.gz
  bcftools norm -f cons.fa -m- --atomize -Oz -o calls.cons.norm.vcf.gz calls.cons.vcf.gz
  tabix -f -p vcf calls.cons.norm.vcf.gz
fi

echo "### 5. VCF liftover cons->ref (chainSwap'd chain) + score vs truth ###"
rm -f calls.lifted.vcf calls.lifted.vcf.gz
CrossMap vcf cons2ref.chain calls.cons.norm.vcf.gz $REF calls.lifted.vcf > crossmap.vcf.log 2>&1
bgzip -f calls.lifted.vcf
bcftools norm -f $REF -m- --atomize -Oz -o calls.lifted.norm.vcf.gz calls.lifted.vcf.gz
tabix -f -p vcf calls.lifted.norm.vcf.gz
python $LIB/compare_vcfs.py --truth-vcf truth.norm.vcf.gz --called-vcf calls.lifted.norm.vcf.gz \
    --region 19 --coverage 30 --divergence $D --heterozygosity $H --seed $SEED \
    --caller bcftools --iteration 2 --out-prefix vcfscore

echo "### 6. BAM liftover cons->ref + read classification vs truth ###"
rm -f map.lifted*.bam
CrossMap bam cons2ref.chain map.cons.bam map.lifted.bam > crossmap.bam.log 2>&1
LB=$(ls -1 map.lifted*.bam 2>/dev/null | grep -v unmap | head -1)
samtools sort -@ 8 -o map.lifted.sort.bam "$LB"; samtools index map.lifted.sort.bam
python $LIB/compare_bams.py --truth-bam truth.oa.bam --mapped-bam map.lifted.sort.bam \
    --mason-vcf truth.vcf.gz --reference $REF \
    --region 19 --coverage 30 --divergence $D --heterozygosity $H --seed $SEED \
    --caller bcftools --iteration 2 --out-prefix bamscore

echo
echo "############## VALIDATION SUMMARY ##############"
echo "### perfect consensus => het SNPs should recover ~fully at correct REF coords,"
echo "### hom div (baked, not re-called) shows as recovered-in-consensus not here,"
echo "### and BAM not_in_truth/missing ~ 0 with high exact = liftover loses nothing."
echo "--- VCF liftover (compare_vcfs, lifted calls vs GT-proper truth) ---"
column -s, -t vcfscore-vcf-summary.csv
echo "--- BAM liftover (compare_bams, lifted reads vs truth) ---"
column -s, -t bamscore-bam-summary.csv
echo "############## END #### "
