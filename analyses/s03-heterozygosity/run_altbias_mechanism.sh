#!/bin/bash
# s03 mechanism dive: alt-read deficit vs divergence, from the s01 (SNP-only, identity-coord)
# BAMs. iter1 is single-shot to the reference (caller-independent -> one caller suffices);
# iter3 is measured per caller to show the correction (and its caller-dependence). Emits a
# tidy CSV consumed by the "mechanism" subsection of s03-results.rmd. Existing BAMs only.
set -uo pipefail
cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations
LIB=analyses/lib/altbias_mechanism.py
DATA=data/s01-baseline-iterative
OUT=analyses/s03-heterozygosity/results/s03-altbias-mechanism.csv
STRIDE=25
DIVS="0.00 0.02 0.04 0.06 0.08 0.10"

# data dictionary for $OUT. Optional $1 = --backfill (table exists, not being regenerated).
# Run `bash analyses/s03-heterozygosity/run_altbias_mechanism.sh --dict-only` to write only this.
write_dict() {
    python analyses/lib/provenance.py --table "$OUT" --script analyses/s03-heterozygosity/run_altbias_mechanism.sh \
        --command "bash analyses/s03-heterozygosity/run_altbias_mechanism.sh" ${1:-} <<EOF
s03 mechanism: is reference bias at het sites alt-read EXCLUSION or degradation? At every
${STRIDE}th true heterozygous SNP (stride $STRIDE over the het sites of each s01 truth VCF,
$DATA/sim/30X/19/<D>/0.005/42/mm39-truth.vcf.gz), the primary reads of that condition's mapped BAM
(iter<k>/mm39_mapped.bam) with an ALIGNED base at the het position are split by the base they carry
(ref vs alt). Reads with the het position inside a deletion or soft-clip have no aligned base and are
not counted. Secondary / supplementary / duplicate / unmapped reads are skipped. NO MAPQ filter.
Rows: iteration 1 (single-shot to the reference; caller-independent, so bcftools stands in) for each
D, then iteration 3 for each caller x D. Produced by analyses/lib/altbias_mechanism.py --csv, one
call per row (a missing BAM is reported on stderr and its row is absent).

Fields:
  divergence       : simulated hom-div divergence D (s01 SNP-only)
  iteration        : mapping iteration (1 = single-shot, 3 = after two refinement rounds)
  caller           : caller whose consensus the iteration-3 reads were mapped to
  n_sites          : het sites sampled
  n_ref, n_alt     : reads presenting the ref / alt base at the het position
  alt_ref_ratio    : n_alt / n_ref (1.0 = no loss; < 1 = alt-bearing reads missing)
  mapq_ref/alt     : mean MAPQ of those ref / alt reads
  pct_lowmapq_ref/alt : % of ref / alt reads with MAPQ < 20 (would be dropped by a MAPQ-20 filter)
  pct_clip_ref/alt : % of ref / alt reads with any soft-clip in the CIGAR
EOF
}
if [ "${1:-}" = "--dict-only" ]; then write_dict --backfill --replace; exit 0; fi

echo "divergence,iteration,caller,n_sites,n_ref,n_alt,alt_ref_ratio,mapq_ref,mapq_alt,pct_lowmapq_ref,pct_lowmapq_alt,pct_clip_ref,pct_clip_alt" > "$OUT"
run() { # div iter caller
  local D=$1 IT=$2 C=$3
  local B=$DATA/iter/$C/30X/19/$D/0.005/42/iter$IT/mm39_mapped.bam
  local T=$DATA/sim/30X/19/$D/0.005/42/mm39-truth.vcf.gz
  [ -f "$B" ] || { echo "  MISSING $B" >&2; return; }
  python "$LIB" "$B" "$T" 19 --stride $STRIDE --label "$D,$IT,$C" --csv >> "$OUT"
  echo "  done D=$D iter$IT $C" >&2
}
for D in $DIVS; do run "$D" 1 bcftools; done           # single-shot: caller-independent
for D in $DIVS; do for C in bcftools gatk freebayes; do run "$D" 3 "$C"; done; done
write_dict
echo "WROTE $OUT ($(($(wc -l < "$OUT")-1)) rows)" >&2
