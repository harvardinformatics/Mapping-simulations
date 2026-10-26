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

# data dictionary for the combined phantom-by-divergence.csv (the per-divergence outputs get theirs
# from effective_truth.py). Optional $1 = --backfill (table exists, not being regenerated).
# `bash analyses/s02-indels/run_effective_truth.sh --dict-only` writes only this dictionary.
write_dict() {
    python analyses/lib/provenance.py --table "$OUT/phantom-by-divergence.csv" \
        --script analyses/s02-indels/run_effective_truth.sh \
        --command "bash analyses/s02-indels/run_effective_truth.sh" \
        --note "$(python -c 'import sys; sys.path.insert(0, "analyses/lib"); import effective_truth; print(effective_truth.ARTIFACT_NOTE)')" \
        "$@" <<EOT
The per-divergence eff-<D>-summary.csv tables (effective_truth.py, D = 0.04 0.06 0.08 0.10, seed
$SEED) concatenated under one header, one row per divergence that had a truth BAM and truth VCF.

Fields: identical to eff-<D>-summary.csv; see eff-<D>-summary.csv.dict.txt next to this file.
EOT
}
if [ "${1:-}" = "--dict-only" ]; then write_dict --backfill --replace; exit 0; fi

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
write_dict
echo "### combined -> $OUT/phantom-by-divergence.csv ###"; column -t -s, "$OUT/phantom-by-divergence.csv"
