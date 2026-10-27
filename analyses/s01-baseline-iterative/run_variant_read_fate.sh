#!/bin/bash
#SBATCH --job-name=s01_vrf
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=32g
#SBATCH --time=12:00:00
#SBATCH --output=../../snakemake/slurm-logs/s01_vrf-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gthomas@g.harvard.edu

# s01 mechanism dive: per-variant CALL fate x READ-placement fate across iterations.
# For one divergence: regenerate the mason -oa truth BAM (temp() in the workflow, so gone),
# then for each caller join each true divergent SNP's recovered? status with the placement
# classes of the reads that truly cover it, at iters 1/2/3 (variant_read_fate.py). Shows the
# causal chain: SNPs get recovered once their reads stop being mismapped/clipped.
# Usage:  sbatch run_variant_read_fate.sh <DIV>
set -eo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims
cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations

D=$1; SEED=42; READLEN=150; FRAG=300; SAMPLE=8000
BASE=data/s01-baseline-iterative
REF=$BASE/ref-regions/mm39-19.fa
SIM=$BASE/sim/30X/19/$D/0.005/$SEED
MTRUTH=$SIM/mm39-truth.vcf.gz            # mason-format (for mason -iv, and the true read spans)
NTRUTH=$SIM/mm39-truth.norm.vcf.gz       # normalized (for SNP loading)
WORK=data/varfate/$D; mkdir -p $WORK
TB=$WORK/truth.oa.bam
OUT=analyses/s01-baseline-iterative/results/s01-varfate-$D.csv

# --- regenerate the mason -oa truth BAM (same seed -> identical reads) ---
if [ ! -s "$TB" ]; then
  LEN=$(cut -f2 $REF.fai | head -1); NFRAG=$(( 30 * LEN / (2 * READLEN) ))
  TMP=$(mktemp -d)
  mason_simulator -ir $REF -iv $MTRUTH -n $NFRAG \
      --illumina-read-length $READLEN --fragment-mean-size $FRAG \
      --num-threads 8 --seed $SEED \
      -o $TMP/r1.fq.gz -or $TMP/r2.fq.gz -oa $TMP/truth.sam --version-check 0
  samtools view -b $TMP/truth.sam -o $TB; rm -rf $TMP
fi

# --- per-variant fate for each caller (3 iterations in one invocation) ---
first=1
for C in bcftools gatk freebayes; do
  it_args=""
  for k in 1 2 3; do
    d=$BASE/iter/$C/30X/19/$D/0.005/$SEED/iter$k
    it_args="$it_args --iter $k:$d/mm39_mapped.bam:$d/mm39_recovered.vcf.gz"
  done
  python analyses/lib/variant_read_fate.py --truth-bam $TB --truth-vcf $NTRUTH \
      --mason-vcf $MTRUTH --reference $REF \
      --region 19 --divergence $D --caller $C --sample $SAMPLE $it_args --out $WORK/vrf-$C.csv
  if [ $first -eq 1 ]; then cp $WORK/vrf-$C.csv $OUT; first=0
  else tail -n +2 $WORK/vrf-$C.csv >> $OUT; fi
done

# data dictionary for the combined table: this script's provenance block + the per-SNP ledger
# definitions variant_read_fate.py wrote for the per-caller tables (identical for every caller)
sed -n '/^Sites      :/,$p' $WORK/vrf-bcftools.csv.dict.txt | \
  python analyses/lib/provenance.py --table $OUT --script analyses/s01-baseline-iterative/run_variant_read_fate.sh \
      --command "cd analyses/s01-baseline-iterative && sbatch run_variant_read_fate.sh $D" \
      --note "rows of the three per-caller tables $WORK/vrf-{bcftools,gatk,freebayes}.csv concatenated, in that order"
echo "WROTE $OUT ($(($(wc -l < $OUT)-1)) rows)"
