#!/bin/bash
#SBATCH --job-name=site_linkage
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=32g
#SBATCH --time=6:00:00
#SBATCH --output=snakemake/slurm-logs/site_linkage-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gthomas@g.harvard.edu

# Run read_site_linkage.py across all divergence levels with divergence sites
# (D>0), then concatenate the per-condition per-site tables into one combined CSV
# in summary-data/ for the notebook. Runs as its own SLURM job so it survives a
# disconnected session.

set -euo pipefail

source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims

PROJ=/n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations
cd "$PROJ"

REGION=19
COV=30
HET=0.005
DIVS=("0.02" "0.04" "0.06" "0.08" "0.10")   # 0.00 has no divergence sites
COMBINED=summary-data/mm39-mason-chr19-site-linkage.csv

per_div_files=()
for div in "${DIVS[@]}"; do
    base=data/mm39-mason-chr19
    out=$base/summary-files/${COV}X/${REGION}/${div}/${HET}/mm39-site-linkage.csv
    echo "### divergence $div -> $out"
    python scripts/read_site_linkage.py \
        --truth-vcf  $base/truth/${REGION}/${div}/${HET}/mm39-${REGION}-${div}d-${HET}h-truth.vcf.gz \
        --called-vcf $base/called-variants/${COV}X/${REGION}/${div}/${HET}/mm39_called-snps.vcf.gz \
        --truth-bam  $base/simulated-reads/${COV}X/${REGION}/${div}/${HET}/mm39_truth.bam \
        --mapped-bam $base/mapped-reads/${COV}X/${REGION}/${div}/${HET}/mm39_mapped.bam \
        --region $REGION --coverage $COV --divergence $div --heterozygosity $HET \
        --n-per-class 10000 --seed 42 \
        --out "$out"
    per_div_files+=("$out")
done

# concatenate: header from the first file, data rows from all
python scripts/combine_summaries.py --inputs "${per_div_files[@]}" --out "$COMBINED" --kind site-linkage
echo "DONE -> $COMBINED"
