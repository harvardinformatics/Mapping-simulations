#!/bin/bash
#SBATCH --job-name=s02_vrf_agg
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=40g
#SBATCH --time=10:00:00
#SBATCH --output=../../snakemake/slurm-logs/s02_vrf_agg-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gthomas@g.harvard.edu

# Full-chromosome per-SNP read-fate ledger for s02 (SNPs IN THE PRESENCE OF INDELS), aggregate mode.
# Same analysis as s01, but s02's consensus is length-shifted by recovered indels, so at iter>1 the
# reads are in CONSENSUS coordinates -- we pass the previous iteration's ref2cons.chain and the tool
# lifts reference coords (each SNP position, each read's true start) into consensus coords before
# classifying at p. Recovered set = the accumulated reference-coord calls (mm39_recovered.acc.vcf.gz).
# One condition per job:  sbatch run_vrf_aggregate.sh <CALLER> <DIV>
set -eo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims
cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations

C=$1; D=$2; SEED=42; HET=0.005
BASE=data/s02-indels
S=$BASE/sim/30X/19/$D/$HET/$SEED
TB=$S/mm39-truth.oa.bam
[ -s "$TB" ] || { echo "ERROR: truth BAM missing: $TB"; exit 1; }
NTRUTH=$S/mm39-truth.norm.vcf.gz
OUT=analyses/s02-indels/results/vrf-agg
mkdir -p "$OUT"

# build --iter specs: iter1 identity coords (no chain); iter>1 gets the PREVIOUS iter's ref2cons.chain
IT=""
for k in 1 2 3; do
    d=$BASE/iter/$C/30X/19/$D/$HET/$SEED/iter$k
    spec="$k:$d/mm39_mapped.bam:$d/mm39_recovered.acc.vcf.gz"
    if [ "$k" -gt 1 ]; then
        prev=$BASE/iter/$C/30X/19/$D/$HET/$SEED/iter$((k - 1))
        spec="$spec:$prev/ref2cons.chain"
    fi
    IT="$IT --iter $spec"
done

echo "### s02 aggregate SNP read-fate: caller=$C div=$D (full chromosome, chain-lifted iter>1) ###"
python analyses/lib/variant_read_fate.py --truth-bam "$TB" --truth-vcf "$NTRUTH" \
    --region 19 --divergence "$D" --caller "$C" --sample 0 --aggregate $IT \
    --out "$OUT/vrf-agg-$C-$D"
echo "### DONE caller=$C div=$D ###"
