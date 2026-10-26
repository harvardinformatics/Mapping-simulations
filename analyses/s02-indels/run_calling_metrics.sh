#!/bin/bash
#SBATCH --job-name=s02_calling_metrics
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=16g
#SBATCH --time=2:00:00
#SBATCH --output=slurm-logs/s02_calling_metrics-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gthomas@g.harvard.edu

# Per-SNP CALLING-side covariates for s02 (indels present). Runs analyses/lib/calling_metrics.py
# over each divergence's truth VCF and concatenates into one table. These are STATIC site
# properties (window +/-100 bp): homopolymer_frac, indel_density, dist_to_indel, snp_density,
# het_density -- the local sequence/variant context that gates gatk's local reassembler.
# Downstream they are compared BY FATE (see scripts/aggregate_calling_metrics.R).
set -eo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims
cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations

REF=data/s02-indels/ref-regions/mm39-19.fa
OUTDIR=analyses/s02-indels/results/calling-metrics
mkdir -p "$OUTDIR"
COMBINED=$OUTDIR/calling-metrics.csv

# per-divergence truth VCFs (mason-format IDs: div_snp / div_indel / het_snp)
DIVS="0.06 0.08 0.10"

# data dictionary for the combined table (each per-divergence part gets its own from
# calling_metrics.py). Optional $1 = --backfill (table exists, not being regenerated).
# `bash analyses/s02-indels/run_calling_metrics.sh --dict-only` writes only this dictionary.
write_dict() {
    python analyses/lib/provenance.py --table "$COMBINED" --script analyses/s02-indels/run_calling_metrics.sh \
        --command "cd analyses/s02-indels && sbatch run_calling_metrics.sh" "$@" <<EOF
s02 per-SNP calling-side covariates, all divergences in one table: the per-divergence tables
calling-metrics-<D>.csv (D = $DIVS; analyses/lib/calling_metrics.py, window +/-100 bp, homopolymer
run >= 5 bp, on data/s02-indels/sim/30X/19/<D>/0.005/42/mm39-truth.vcf.gz) concatenated in that
order under one header. One row per divergent (1|1) truth SNP per divergence.

Fields: identical to the per-divergence tables; see calling-metrics-<D>.csv.dict.txt next to
this file (committed to git, while the tables themselves are git-ignored for size).
EOF
}
if [ "${1:-}" = "--dict-only" ]; then
    shift
    write_dict --backfill --replace "$@"
    exit 0
fi

first=1
for D in $DIVS; do
    VCF=data/s02-indels/sim/30X/19/$D/0.005/42/mm39-truth.vcf.gz
    [ -s "$VCF" ] || { echo "ERROR: missing truth VCF $VCF"; exit 1; }
    PART=$OUTDIR/calling-metrics-$D.csv
    echo "### calling_metrics.py D=$D ###"
    python analyses/lib/calling_metrics.py \
        --truth-vcf "$VCF" --reference "$REF" --divergence "$D" \
        --window 100 --hp-min 5 --out "$PART"

    # concatenate: keep the header from the first divergence only
    if [ "$first" -eq 1 ]; then
        cp "$PART" "$COMBINED"; first=0
    else
        tail -n +2 "$PART" >> "$COMBINED"
    fi
done

write_dict
echo "### wrote $COMBINED ($(wc -l < "$COMBINED") lines incl. header) ###"
