#!/bin/bash
# Build the EFFECTIVE truth (truth minus phantom SNPs mason never simulated) for each s02 divergence.
# Same tool as s01; s02 truth VCFs also contain indels, which effective_truth.py passes through to the
# effective VCF unchanged (only SNP read-level presence is checked). One truth-BAM scan per divergence.
# Run:  bash run_effective_truth.sh
set -eo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims
cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations

SEED=42; HET=0.005
BASE=data/s02-indels
OUT=analyses/s02-indels/results/effective-truth
mkdir -p "$OUT"

for D in 0.04 0.06 0.08 0.10; do
    S=$BASE/sim/30X/19/$D/$HET/$SEED
    TB=$S/mm39-truth.oa.bam
    NT=$S/mm39-truth.norm.vcf.gz
    if [ ! -s "$TB" ] || [ ! -s "$NT" ]; then echo "SKIP D=$D (missing truth BAM/VCF)"; continue; fi
    echo "### s02 effective truth D=$D ###"
    python analyses/lib/effective_truth.py --truth-bam "$TB" --truth-vcf "$NT" \
        --region 19 --divergence "$D" --seed $SEED --out-prefix "$OUT/eff-$D"
done

head -1 "$OUT"/eff-0.10-summary.csv > "$OUT/phantom-by-divergence.csv"
for D in 0.04 0.06 0.08 0.10; do
    [ -s "$OUT/eff-$D-summary.csv" ] && tail -n +2 "$OUT/eff-$D-summary.csv" >> "$OUT/phantom-by-divergence.csv"
done
echo "### combined -> $OUT/phantom-by-divergence.csv ###"; column -t -s, "$OUT/phantom-by-divergence.csv"
