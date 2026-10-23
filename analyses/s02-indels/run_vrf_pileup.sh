#!/bin/bash
#SBATCH --job-name=s02_vrf_pile
#SBATCH --partition=shared
#SBATCH --cpus-per-task=2
#SBATCH --mem=32g
#SBATCH --time=8:00:00
#SBATCH --output=../../snakemake/slurm-logs/s02_vrf_pile-%j.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gthomas@g.harvard.edu

# s02 companion decomposition: per-SNP read fate + position-centric pileup at p (--pileup), same as s01
# but for SNPs IN THE PRESENCE OF INDELS. At iter>1 the consensus is length-shifted, so we pass the
# previous iteration's ref2cons.chain and the tool lifts reference coords -> consensus coords before
# both the ledger classification AND the pileup fetch; the mismapped-in check stays in reference coords
# (truth is reference), so contamination is measured correctly. Sampled (~30k SNPs).
# One condition per job:  sbatch run_vrf_pileup.sh <CALLER> <DIV>
set -eo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims
cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations

C=$1; D=$2; SEED=42; HET=0.005; SAMPLE=30000
BASE=data/s02-indels
S=$BASE/sim/30X/19/$D/$HET/$SEED
TB=$S/mm39-truth.oa.bam
[ -s "$TB" ] || { echo "ERROR: truth BAM missing: $TB"; exit 1; }
NTRUTH=$S/mm39-truth.norm.vcf.gz
OUT=analyses/s02-indels/results/vrf-pileup
mkdir -p "$OUT"

# iter1 identity coords (no chain); iter>1 gets the PREVIOUS iter's ref2cons.chain
IT=""
for k in 1 2 3; do
    d=$BASE/iter/$C/30X/19/$D/$HET/$SEED/iter$k
    spec="$k:$d/mm39_mapped.bam:$d/mm39_recovered.acc.vcf.gz"
    if [ "$k" -gt 1 ]; then
        spec="$spec:$BASE/iter/$C/30X/19/$D/$HET/$SEED/iter$((k - 1))/ref2cons.chain"
    fi
    IT="$IT --iter $spec"
done

echo "### s02 per-SNP read-fate + pileup: caller=$C div=$D (sample $SAMPLE, chain-lifted iter>1) ###"
python analyses/lib/variant_read_fate.py --truth-bam "$TB" --truth-vcf "$NTRUTH" \
    --region 19 --divergence "$D" --caller "$C" --sample $SAMPLE --pileup $IT \
    --out "$OUT/vrf-pileup-$C-$D.csv"
echo "### DONE caller=$C div=$D ###"
