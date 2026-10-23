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

echo "### wrote $COMBINED ($(wc -l < "$COMBINED") lines incl. header) ###"
