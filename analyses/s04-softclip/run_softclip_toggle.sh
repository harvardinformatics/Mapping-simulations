#!/bin/bash
#SBATCH --job-name=s04_sc
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=48g
#SBATCH --time=16:00:00
#SBATCH --output=../../snakemake/slurm-logs/s04_sc-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gthomas@g.harvard.edu

# s04 soft-clip test: does gatk's high-divergence SNP loss depend on soft-clipped bases?
# Re-run HaplotypeCaller on the SAME existing iter1 BAM twice -- default (uses soft-clipped
# bases) vs --dont-use-soft-clipped-bases -- and score SNP recall vs truth. Paired on one BAM,
# so any difference is the flag alone. Prediction (assembler-not-clipping): recall barely moves.
# Usage:  sbatch run_softclip_toggle.sh <DIV>
set -eo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims
cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations

D=$1
SEED=42; CONTIG=19; N=8
BAM=data/s01-baseline-iterative/iter/gatk/30X/19/$D/0.005/$SEED/iter1/mm39_mapped.bam
REF=data/s01-baseline-iterative/ref-regions/mm39-19.fa
TRUTH=data/s01-baseline-iterative/sim/30X/19/$D/0.005/$SEED/mm39-truth.norm.vcf.gz
WORK=data/s04-softclip/$D
OUT=analyses/s04-softclip/results/softclip-toggle-$D.csv
LEN=$(cut -f2 $REF.fai | head -1); CHUNK=$(( (LEN + N - 1) / N ))

call_gatk() {  # $1 = mode label ; $2 = extra HC flags
  local mode=$1 flag=$2 d=$WORK/$1
  rm -rf $d; mkdir -p $d/scatter
  local pids=""
  for i in $(seq 1 $N); do
    local s=$(( (i-1)*CHUNK + 1 )) e=$(( i*CHUNK )); [ $e -gt $LEN ] && e=$LEN
    gatk HaplotypeCaller -R $REF -I $BAM -ERC GVCF -L $CONTIG:$s-$e \
        -O $d/scatter/g_$i.vcf.gz --java-options -Xmx4g $flag > $d/hc_$i.log 2>&1 &
    pids="$pids $!"
  done
  local fail=0; for p in $pids; do wait $p || fail=1; done; [ $fail -eq 0 ]
  local gargs=""; for i in $(seq 1 $N); do gargs="$gargs -I $d/scatter/g_$i.vcf.gz"; done
  gatk GatherVcfs $gargs -O $d/g.vcf.gz > $d/gather.log 2>&1; tabix -p vcf $d/g.vcf.gz
  gatk GenotypeGVCFs -R $REF -V $d/g.vcf.gz -O $d/raw.vcf.gz > $d/geno.log 2>&1
  bcftools norm -f $REF -m- --atomize -Oz -o $d/calls.norm.vcf.gz $d/raw.vcf.gz 2> $d/norm.log
  python analyses/lib/compare_vcfs.py --truth-vcf $TRUTH --called-vcf $d/calls.norm.vcf.gz \
      --region $CONTIG --coverage 30 --divergence $D --heterozygosity 0.005 --seed $SEED \
      --caller gatk_$mode --iteration 1 --out-prefix $d/score
  rm -rf $d/scatter
}

echo "### D=$D : default (uses soft-clipped bases) ###"
call_gatk default ""
echo "### D=$D : --dont-use-soft-clipped-bases ###"
call_gatk nosoftclip "--dont-use-soft-clipped-bases"

# combine the two score rows into one comparison CSV
head -1 $WORK/default/score-vcf-summary.csv > $OUT
tail -n +2 $WORK/default/score-vcf-summary.csv    >> $OUT
tail -n +2 $WORK/nosoftclip/score-vcf-summary.csv >> $OUT
echo "WROTE $OUT"; column -t -s, $OUT
