#!/bin/bash
#SBATCH --job-name=s03_vrf_pile
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=32g
#SBATCH --time=8:00:00
#SBATCH --output=../../snakemake/slurm-logs/s03_vrf_pile-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gthomas@g.harvard.edu

# s03 companion decomposition: per-variant read fate + pileup at p for HETEROZYGOUS sites (--het).
# Runs on the s01 (SNP-only) data so coordinates are identity across iterations (no chain). n_true_alt
# = the alt-HAPLOTYPE reads (~half); "recovered" = the caller called the site correctly AS A HET
# (GT 0/1) in its per-iteration VCF (mm39_calls.norm.vcf.gz). NB for het, ref reads at p are REAL
# (the ref haplotype), so the contamination axis to watch is mismapped-IN ref, not raw ref_p.
# One condition per job:  sbatch run_vrf_pileup.sh <CALLER> <DIV>
set -eo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims
cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations

C=$1; D=$2; SEED=42; SAMPLE=30000
TB=data/varfate/$D/truth.oa.bam                          # cached mason -oa truth BAM (has oH)
[ -s "$TB" ] || { echo "ERROR: truth BAM missing: $TB"; exit 1; }
BASE=data/s01-baseline-iterative
NTRUTH=$BASE/sim/30X/19/$D/0.005/$SEED/mm39-truth.norm.vcf.gz
MTRUTH=$BASE/sim/30X/19/$D/0.005/$SEED/mm39-truth.vcf.gz      # mason-format (true read spans)
REF=$BASE/ref-regions/mm39-19.fa
OUT=analyses/s03-heterozygosity/results/vrf-pileup
mkdir -p "$OUT"

IT=""
for k in 1 2 3; do
    d=$BASE/iter/$C/30X/19/$D/0.005/$SEED/iter$k
    IT="$IT --iter $k:$d/mm39_mapped.bam:$d/mm39_calls.norm.vcf.gz"
done

echo "### s03 per-variant HET read-fate + pileup: caller=$C div=$D (sample $SAMPLE) ###"
python analyses/lib/variant_read_fate.py --het --truth-bam "$TB" --truth-vcf "$NTRUTH" \
    --mason-vcf "$MTRUTH" --reference "$REF" \
    --region 19 --divergence "$D" --caller "$C" --sample $SAMPLE --pileup $IT \
    --out "$OUT/vrf-pileup-$C-$D.csv"
echo "### DONE caller=$C div=$D ###"
