#!/bin/bash
#SBATCH --job-name=s05_vrf_pile
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=32g
#SBATCH --time=10:00:00
#SBATCH --output=../../snakemake/slurm-logs/s05_vrf_pile-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gthomas@g.harvard.edu

# s05 mapping-vs-calling decomposition: per-divergent-SNP read-fate ledger + position-centric
# pileup at p, for the reference characterization config (L20000-S1.0, mean D=0.06). For every
# sampled truth divergent site the ledger records (a) whether it was recovered, (b) alt_p = alt-
# supporting reads present in the caller's pileup at p, and (c) the truth-centric read fate
# (where the alt-carrying reads actually went). Downstream this splits each unrecovered site into
# MAPPING-limited (alt evidence absent at p) vs CALLING-limited (alt present but uncalled), as a
# function of local divergence.
#
# One caller per job:  sbatch run_vrf_pileup.sh <gatk|bcftools|freebayes>
set -eo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims
cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations

C=$1
[ -n "$C" ] || { echo "usage: sbatch run_vrf_pileup.sh <caller>"; exit 1; }
FIELD=L20000-S1.0; DIV=0.06; HET=0.005; SEED=42; SAMPLE=50000
D=data/s05-heterogeneous-divergence
TB=$D/varfate/$FIELD/truth.oa.bam                         # mason -oa truth BAM (has oH)
NTRUTH=$D/sim/30X/19/$DIV/$HET/$FIELD/$SEED/mm39-truth.norm.vcf.gz
OUT=analyses/s05-heterogeneous-divergence/results/vrf-pileup
mkdir -p "$OUT"

[ -s "$TB" ]     || { echo "ERROR: truth BAM missing: $TB (run run_vrf_truthbam.sh first)"; exit 1; }
[ -s "$NTRUTH" ] || { echo "ERROR: normalized truth VCF missing: $NTRUTH"; exit 1; }

# recovered set = the cumulative consensus-vs-reference SNP diff (mm39_recovered.vcf.gz), same
# as the recovery scoring; one --iter spec per iteration.
IT=""
for k in 1 2 3; do
    d=$D/iter/$C/30X/19/$DIV/$HET/$FIELD/$SEED/iter$k
    IT="$IT --iter $k:$d/mm39_mapped.bam:$d/mm39_recovered.vcf.gz"
done

echo "### s05 read-fate + pileup: caller=$C field=$FIELD D=$DIV (sample $SAMPLE) ###"
python analyses/lib/variant_read_fate.py --truth-bam "$TB" --truth-vcf "$NTRUTH" \
    --region 19 --divergence "$DIV" --caller "$C" --sample $SAMPLE --pileup $IT \
    --out "$OUT/vrf-pileup-$C.csv"
echo "### DONE caller=$C ###"
