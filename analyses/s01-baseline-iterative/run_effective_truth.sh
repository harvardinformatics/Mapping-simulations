#!/bin/bash
# Build the EFFECTIVE truth (truth minus phantom SNPs mason never simulated) for each s01 divergence,
# so recall can be scored against variants that actually exist in the reads. One truth-BAM scan per
# divergence (seed 42). Outputs -> results/effective-truth/. Run:  bash run_effective_truth.sh
set -eo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims
cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations

SEED=42
BASE=data/s01-baseline-iterative
OUT=analyses/s01-baseline-iterative/results/effective-truth
mkdir -p "$OUT"

for D in 0.04 0.06 0.08 0.10; do
    TB=data/varfate/$D/truth.oa.bam
    NT=$BASE/sim/30X/19/$D/0.005/$SEED/mm39-truth.norm.vcf.gz
    if [ ! -s "$TB" ] || [ ! -s "$NT" ]; then echo "SKIP D=$D (missing truth BAM/VCF)"; continue; fi
    echo "### effective truth D=$D ###"
    python analyses/lib/effective_truth.py --truth-bam "$TB" --truth-vcf "$NT" \
        --region 19 --divergence "$D" --seed $SEED --out-prefix "$OUT/eff-$D"
done

# combine the per-divergence summaries into one table
head -1 "$OUT"/eff-0.10-summary.csv > "$OUT/phantom-by-divergence.csv"
for D in 0.04 0.06 0.08 0.10; do
    [ -s "$OUT/eff-$D-summary.csv" ] && tail -n +2 "$OUT/eff-$D-summary.csv" >> "$OUT/phantom-by-divergence.csv"
done
echo "### combined -> $OUT/phantom-by-divergence.csv ###"
column -t -s, "$OUT/phantom-by-divergence.csv"
