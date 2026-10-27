#!/bin/bash
#SBATCH --job-name=s05_vrf_pile_cap
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=32g
#SBATCH --time=10:00:00
#SBATCH --output=../../snakemake/slurm-logs/s05_vrf_pile_cap-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gthomas@g.harvard.edu

# mapping-vs-calling ledger on the REALISTIC CAPPED config (L20000-S1.0-C0.20, mean D=0.03,
# local divergence capped at 0.20). One caller per job: sbatch run_vrf_pileup_capped.sh <caller>
set -eo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims
cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations

C=$1
[ -n "$C" ] || { echo "usage: sbatch run_vrf_pileup_capped.sh <caller>"; exit 1; }
FIELD=L20000-S1.0-C0.20; DIV=0.03; HET=0.005; SEED=42; SAMPLE=50000
D=data/s05-heterogeneous-divergence
TB=$D/varfate/$FIELD/truth.oa.bam
NTRUTH=$D/sim/30X/19/$DIV/$HET/$FIELD/$SEED/mm39-truth.norm.vcf.gz
MTRUTH=$D/sim/30X/19/$DIV/$HET/$FIELD/$SEED/mm39-truth.vcf.gz   # mason-format (true read spans)
REF=$D/ref-regions/mm39-19.fa
OUT=analyses/s05-heterogeneous-divergence/results/vrf-pileup-capped
mkdir -p "$OUT"

[ -s "$TB" ]     || { echo "ERROR: truth BAM missing: $TB (run run_vrf_truthbam_capped.sh first)"; exit 1; }
[ -s "$NTRUTH" ] || { echo "ERROR: normalized truth VCF missing: $NTRUTH"; exit 1; }

IT=""
for k in 1 2 3; do
    d=$D/iter/$C/30X/19/$DIV/$HET/$FIELD/$SEED/iter$k
    IT="$IT --iter $k:$d/mm39_mapped.bam:$d/mm39_recovered.vcf.gz"
done

echo "### s05 CAPPED read-fate + pileup: caller=$C field=$FIELD D=$DIV (sample $SAMPLE) ###"
python analyses/lib/variant_read_fate.py --truth-bam "$TB" --truth-vcf "$NTRUTH" \
    --mason-vcf "$MTRUTH" --reference "$REF" \
    --region 19 --divergence "$DIV" --caller "$C" --sample $SAMPLE --pileup $IT \
    --out "$OUT/vrf-pileup-$C.csv"
echo "### DONE caller=$C ###"
