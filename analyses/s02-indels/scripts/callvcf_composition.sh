#!/bin/bash
# callvcf_composition.sh -- composition of each caller's normalized call VCF in s02 (SNPs + indels)
# at D = 0.10, every iteration: how many records are SNPs, indels and spanning-deletion '*' alleles,
# and how many extra records share a position (split multiallelics). Shows freebayes' record
# inflation and how it shrinks as iteration improves mapping. Used by s02-results.rmd.
#
# Codifies the ad hoc loop that first made the table (2026-08-02); the counting rules and output
# format are unchanged.
#
# Usage (from the project root):
#   bash analyses/s02-indels/scripts/callvcf_composition.sh [OUT]
#   OUT defaults to analyses/s02-indels/results/s02-indels-callvcf-composition.csv
#   bash analyses/s02-indels/scripts/callvcf_composition.sh --dict-only [OUT]
#       only (re)write the data dictionary of an existing table, as a backfill
set -eo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims          # (not `conda run`: it does not pass the heredoc on stdin)
set -u
cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations

# --- settings ------------------------------------------------------------------------------------
SCRIPT=analyses/s02-indels/scripts/callvcf_composition.sh
D=0.10; HET=0.005; SEED=42
CALLERS="gatk bcftools freebayes"; ITERS="1 2 3"
DICT_ONLY=0
if [ "${1:-}" = "--dict-only" ]; then DICT_ONLY=1; shift; fi
OUT=${1:-analyses/s02-indels/results/s02-indels-callvcf-composition.csv}

# --- data dictionary -----------------------------------------------------------------------------
write_dict() {   # optional extra provenance.py flags (e.g. --backfill --note ...)
    python analyses/lib/provenance.py --table "$OUT" --script "$SCRIPT" \
        --command "bash $SCRIPT" "$@" <<EOF
Record composition of the normalized caller VCFs
data/s02-indels/iter/<caller>/30X/19/$D/$HET/$SEED/iter<k>/mm39_calls.norm.vcf.gz
(D = $D; all records, no filtering). One row per (caller, iteration).

Fields:
  caller             : variant caller (gatk | bcftools | freebayes)
  iteration          : mapping iteration (1 = single-shot)
  total              : all records (= snp + indel + star)
  snp                : records with 1-bp REF and 1-bp ALT
  indel              : other records with a non-'*' ALT (REF and/or ALT longer than 1 bp)
  star               : records whose ALT is '*' (spanning deletion)
  multiallelic_extra : records beyond the first at the same CHROM:POS (sum over positions of
                       records - 1), i.e. how many records come from split multiallelic sites
EOF
}
if [ "$DICT_ONLY" -eq 1 ]; then
    write_dict --backfill --replace --note "first made 2026-08-02 by an equivalent ad hoc loop (same awk counting and output format); this script reproduces it byte for byte (checked 2026-10-26)"
    exit 0
fi

# --- count ---------------------------------------------------------------------------------------
# one awk pass per VCF; C and IT label the row
AWK='
!/^#/ {
  tot++; pos[$1"_"$2]++;
  if ($5=="*") star++;
  else if (length($4)==1 && length($5)==1) snp++;
  else indel++;
}
END {
  extra=0; for (p in pos) if (pos[p]>1) extra+=pos[p]-1;
  printf "%s,%s,%d,%d,%d,%d,%d\n", C, IT, tot+0, snp+0, indel+0, star+0, extra+0;
}'
echo "caller,iteration,total,snp,indel,star,multiallelic_extra" > "$OUT"
for c in $CALLERS; do
    for k in $ITERS; do
        f=data/s02-indels/iter/$c/30X/19/$D/$HET/$SEED/iter$k/mm39_calls.norm.vcf.gz
        [ -s "$f" ] || { echo "ERROR: missing call VCF $f" >&2; exit 1; }
        zcat "$f" | awk -v C="$c" -v IT="$k" "$AWK" >> "$OUT"
    done
done
write_dict
echo "wrote $OUT"
