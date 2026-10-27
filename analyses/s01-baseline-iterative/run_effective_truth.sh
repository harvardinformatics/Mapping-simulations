#!/bin/bash
#SBATCH --job-name=s01_truthsites
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=16g
#SBATCH --time=6:00:00
#SBATCH --output=../../snakemake/slurm-logs/s01_truthsites-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gthomas@g.harvard.edu

# Truth-site check for every s01 divergence: are all truth SNPs actually in the simulated reads?
# analyses/lib/effective_truth.py rebuilds mason's two haplotypes from the mason-format truth VCF and
# reads each truth read's base at every truth SNP from its haplotype position (oP), not the -oa CIGAR.
# One truth-BAM scan per divergence (seed 42; truth BAMs under data/varfate/<D>/), then the summaries
# are combined into phantom-by-divergence.csv. Outputs -> results/effective-truth/.
# Run:  cd analyses/s01-baseline-iterative && sbatch run_effective_truth.sh
#       bash analyses/s01-baseline-iterative/run_effective_truth.sh --dict-only   (combined dictionary only)
set -eo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims
cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations

SEED=42
DIVS="0.00 0.001 0.005 0.02 0.04 0.06 0.08 0.10"
BASE=data/s01-baseline-iterative
REF=$BASE/ref-regions/mm39-19.fa
OUT=analyses/s01-baseline-iterative/results/effective-truth
mkdir -p "$OUT"

# data dictionary for the combined phantom-by-divergence.csv (the per-divergence outputs get theirs
# from effective_truth.py). Optional extra provenance.py flags (e.g. --backfill --replace).
write_dict() {
    python analyses/lib/provenance.py --table "$OUT/phantom-by-divergence.csv" \
        --script analyses/s01-baseline-iterative/run_effective_truth.sh \
        --command "cd analyses/s01-baseline-iterative && sbatch run_effective_truth.sh" \
        --note "$(python -c 'import sys; sys.path.insert(0, "analyses/lib"); import effective_truth; print(effective_truth.CORRECTION_NOTE)')" \
        "$@" <<EOT
The per-divergence eff-<D>-summary.csv tables (effective_truth.py, D = $DIVS, seed $SEED)
concatenated under one header, one row per divergence, in that order.

Fields: identical to eff-<D>-summary.csv; see eff-<D>-summary.csv.dict.txt next to this file.
EOT
}
if [ "${1:-}" = "--dict-only" ]; then write_dict --backfill --replace; exit 0; fi

# every input must exist: a missing divergence is an error, not a silent skip
for D in $DIVS; do
    for f in data/varfate/$D/truth.oa.bam $BASE/sim/30X/19/$D/0.005/$SEED/mm39-truth.vcf.gz; do
        [ -s "$f" ] || { echo "ERROR: missing input $f"; exit 1; }
    done
done

for D in $DIVS; do
    echo "### truth-site check D=$D ###"
    python analyses/lib/effective_truth.py --mason-vcf $BASE/sim/30X/19/$D/0.005/$SEED/mm39-truth.vcf.gz \
        --reference $REF --truth-bam data/varfate/$D/truth.oa.bam \
        --region 19 --divergence "$D" --seed $SEED --out-prefix "$OUT/eff-$D"
done

# combine the per-divergence summaries into one table (header from the first)
first=1
for D in $DIVS; do
    if [ $first -eq 1 ]; then head -1 "$OUT/eff-$D-summary.csv" > "$OUT/phantom-by-divergence.csv"; first=0; fi
    tail -n +2 "$OUT/eff-$D-summary.csv" >> "$OUT/phantom-by-divergence.csv"
done
write_dict
echo "### combined -> $OUT/phantom-by-divergence.csv ###"
column -t -s, "$OUT/phantom-by-divergence.csv"
