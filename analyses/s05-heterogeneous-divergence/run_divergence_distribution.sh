#!/bin/bash
# QC figure for the s05 REALISTIC CAPPED config (L20000-S1.0-C0.20, mean D=0.03): per-read and
# per-site divergence distributions that justify the mapping-vs-calling plot's x-axis (drawn to the
# 99.9th percentile). Extracts the two quantities from the actual simulated data, then plots them.
# Run:  bash analyses/s05-heterogeneous-divergence/run_divergence_distribution.sh
set -eo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims
cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations

FIELD=L20000-S1.0-C0.20; DIV=0.03; HET=0.005; SEED=42
D=data/s05-heterogeneous-divergence
TB=$D/varfate/$FIELD/truth.oa.bam
LD=$D/sim/30X/19/$DIV/$HET/$FIELD/$SEED/mm39-truth-locald.tsv
OUT=analyses/s05-heterogeneous-divergence/results/capped-divergence-distribution.png

[ -s "$TB" ] || { echo "ERROR: truth BAM missing: $TB"; exit 1; }
[ -s "$LD" ] || { echo "ERROR: local-divergence tsv missing: $LD"; exit 1; }

TMP=$(mktemp -d)
# per-READ divergence = NM edit distance / read length, for a deterministic 15% sample of primary
# alignments (samtools -s uses a fixed hash seed, so this is reproducible). Grid = 1/readlen.
echo "### extracting per-read divergence (NM/readlen) ###"
samtools view -s 0.15 -F 0x900 "$TB" \
  | awk 'BEGIN{OFS="\t"}
      { nm=-1;
        for(i=12;i<=NF;i++){ if($i ~ /^NM:i:/){ split($i,a,":"); nm=a[3]; break } }
        if(nm<0) next;
        rl=length($10);
        if(rl>0) print nm/rl }' > "$TMP/readdiv.txt"

# per-SITE read-scale local divergence at divergent SNPs (local_d_150 column). Grid = 1/300.
echo "### extracting per-site local_d_150 ###"
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="local_d_150") c=i; next}{print $c}' "$LD" > "$TMP/locald150.txt"

echo "### plotting (bins aligned to the discrete value grid) ###"
conda run -n renv Rscript analyses/s05-heterogeneous-divergence/scripts/plot_divergence_distribution.R \
    "$TMP/readdiv.txt" "$TMP/locald150.txt" "$OUT"
rm -rf "$TMP"
echo "wrote $OUT"
