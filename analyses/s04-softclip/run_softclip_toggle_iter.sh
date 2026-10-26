#!/bin/bash
#SBATCH --job-name=s04_sci
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=48g
#SBATCH --time=16:00:00
#SBATCH --output=../../snakemake/slurm-logs/s04_sci-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gthomas@g.harvard.edu

# s04 soft-clip test, ITERATION >= 2. Does the soft-clip effect (or lack of it) hold once iteration
# has converged the consensus? At iter k the gatk BAM is mapped to consensus_{k-1}, so we call
# HaplotypeCaller against THAT consensus (default vs --dont-use-soft-clipped-bases) and score the
# CUMULATIVE recovered set the way the s01 pipeline defines it:
#   recovered_k = (already-recovered set baked into consensus_{k-1}) UNION (new hom-alt SNP calls)
# The prior set (iter_{k-1}/mm39_recovered.vcf.gz) is identical for both arms -- only iter-k calling
# is toggled -- so any difference in cumulative div recall is the flag alone. s01 is SNP-only, so the
# consensus keeps reference (identity) coordinates and no liftover is needed.
# Usage:  cd analyses/s04-softclip && sbatch run_softclip_toggle_iter.sh <DIV> <ITER>   (ITER >= 2;
#         submit_softclip_toggle_iter.sh submits D = 0.06, 0.10 x ITER = 2, 3)
#         bash analyses/s04-softclip/run_softclip_toggle_iter.sh <DIV> <ITER> --dict-only   (from the
#         project root: only (re)write the data dictionary of an existing output, marked as a backfill)
set -eo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims
cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations

D=$1; K=$2; SEED=42; CONTIG=19; N=8
BASE=data/s01-baseline-iterative
ITERD=$BASE/iter/gatk/30X/19/$D/0.005/$SEED
BAM=$ITERD/iter$K/mm39_mapped.bam                       # reads mapped to consensus_{K-1}
CONSREF=$ITERD/iter$((K - 1))/mm39-consensus.fa         # that consensus == the calling reference
PRIOR=$ITERD/iter$((K - 1))/mm39_recovered.vcf.gz       # already recovered (baked into consensus_{K-1})
TRUTH=$BASE/sim/30X/19/$D/0.005/$SEED/mm39-truth.norm.vcf.gz
WORK=data/s04-softclip/$D/iter$K
OUT=analyses/s04-softclip/results/softclip-toggle-$D-iter$K.csv

# data dictionary for $OUT: provenance block + what the two rows are + compare_vcfs.py's columns.
# Optional $1 = --backfill (the table already exists and is not being regenerated).
write_dict() {
    python analyses/lib/provenance.py --table "$OUT" --script analyses/s04-softclip/run_softclip_toggle_iter.sh \
        --command "cd analyses/s04-softclip && sbatch run_softclip_toggle_iter.sh $D $K" \
        --append-definitions compare_vcfs $1 <<EOF
s04 soft-clip toggle at iteration $K, divergence D = $D (s01 SNP-only simulation, seed $SEED).
GATK HaplotypeCaller was re-run twice on the SAME existing s01 iter$K gatk BAM (reads mapped to
consensus_$((K - 1))), calling against that consensus, once per arm:
  gatk_default    : HaplotypeCaller defaults (soft-clipped bases are used)
  gatk_nosoftclip : HaplotypeCaller --dont-use-soft-clipped-bases
Each arm's new hom-alt SNP calls were merged with the already-recovered set baked into
consensus_$((K - 1)) (identical for both arms), and that CUMULATIVE recovered set was scored
against the truth with compare_vcfs.py -- the s01 pipeline's definition of recovery at iteration $K.
One row per arm (caller column). Het columns are 0 by construction: het calls are never baked into
the haploid consensus, so div recall (div_tp / n_truth_div) is the metric for this test.
Inputs: BAM     $BAM
        CONSREF $CONSREF
        PRIOR   $PRIOR
        TRUTH   $TRUTH
Per-arm work files (git-ignored): $WORK/{default,nosoftclip}/
EOF
}
if [ "$3" = "--dict-only" ]; then write_dict --backfill --replace; exit 0; fi

for f in "$BAM" "$CONSREF" "$PRIOR" "$TRUTH"; do
    [ -s "$f" ] || { echo "ERROR: missing input $f"; exit 1; }
done
LEN=$(cut -f2 "$CONSREF".fai | head -1); CHUNK=$(( (LEN + N - 1) / N ))

# call gatk on the iter-K BAM (against consensus_{K-1}), build cumulative recovered, score vs truth
call_and_score() {   # $1 = mode label ; $2 = extra HC flag
    local mode=$1 flag=$2 d=$WORK/$1
    rm -rf "$d"; mkdir -p "$d"/scatter                   # $d is this run's own scratch dir
    local pids=""
    for i in $(seq 1 $N); do
        local s=$(( (i - 1) * CHUNK + 1 )) e=$(( i * CHUNK )); [ $e -gt $LEN ] && e=$LEN
        gatk HaplotypeCaller -R "$CONSREF" -I "$BAM" -ERC GVCF -L $CONTIG:$s-$e \
            -O "$d"/scatter/g_$i.vcf.gz --java-options -Xmx4g $flag > "$d"/hc_$i.log 2>&1 &
        pids="$pids $!"
    done
    local fail=0; for p in $pids; do wait $p || fail=1; done; [ $fail -eq 0 ]
    local gargs=""; for i in $(seq 1 $N); do gargs="$gargs -I $d/scatter/g_$i.vcf.gz"; done
    gatk GatherVcfs $gargs -O "$d"/g.vcf.gz > "$d"/gather.log 2>&1; tabix -p vcf "$d"/g.vcf.gz
    gatk GenotypeGVCFs -R "$CONSREF" -V "$d"/g.vcf.gz -O "$d"/raw.vcf.gz > "$d"/geno.log 2>&1
    bcftools norm -f "$CONSREF" -m- --atomize -Oz -o "$d"/calls.norm.vcf.gz "$d"/raw.vcf.gz 2> "$d"/norm.log
    tabix -p vcf "$d"/calls.norm.vcf.gz

    # new hom-alt SNP calls (exactly the pipeline's extraction), reheadered to a common sample so
    # they can be concatenated with the prior recovered set into one cumulative recovered VCF.
    echo rec > "$d"/sn.txt
    bcftools view -v snps -g hom -Oz -o "$d"/homsnps.raw.vcf.gz "$d"/calls.norm.vcf.gz
    bcftools reheader -s "$d"/sn.txt "$d"/homsnps.raw.vcf.gz -o "$d"/homsnps.vcf.gz; tabix -p vcf "$d"/homsnps.vcf.gz
    bcftools reheader -s "$d"/sn.txt "$PRIOR" -o "$d"/prior.vcf.gz; tabix -p vcf "$d"/prior.vcf.gz
    bcftools concat -a "$d"/prior.vcf.gz "$d"/homsnps.vcf.gz 2> "$d"/concat.log \
        | bcftools norm -d all -Oz -o "$d"/recovered.vcf.gz 2>> "$d"/concat.log
    tabix -p vcf "$d"/recovered.vcf.gz

    # score cumulative recovered vs truth -> div recall (het is 0 here: het is never baked into the
    # haploid consensus, so div recall is the metric for this test, matching the pipeline).
    python analyses/lib/compare_vcfs.py --truth-vcf "$TRUTH" --called-vcf "$d"/recovered.vcf.gz \
        --region $CONTIG --coverage 30 --divergence "$D" --heterozygosity 0.005 --seed $SEED \
        --caller gatk_$mode --iteration $K --out-prefix "$d"/score
    rm -rf "$d"/scatter
}

echo "### D=$D iter=$K : default (uses soft-clipped bases) ###"
call_and_score default ""
echo "### D=$D iter=$K : --dont-use-soft-clipped-bases ###"
call_and_score nosoftclip "--dont-use-soft-clipped-bases"

# combine the two score rows into one comparison CSV
head -1 "$WORK"/default/score-vcf-summary.csv > "$OUT"
tail -n +2 "$WORK"/default/score-vcf-summary.csv    >> "$OUT"
tail -n +2 "$WORK"/nosoftclip/score-vcf-summary.csv >> "$OUT"
write_dict
echo "WROTE $OUT"; column -t -s, "$OUT"
