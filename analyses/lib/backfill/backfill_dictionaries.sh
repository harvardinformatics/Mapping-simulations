#!/bin/bash
# backfill_dictionaries.sh -- the record of how the data dictionaries of OLDER outputs were written.
#
# Tables made before provenance tracking (2026-10-26) had no dictionary, or one without its producing
# command. This script writes those dictionaries WITHOUT regenerating any table: each producing tool
# or run script is called in its dictionary-only mode with the command that ORIGINALLY made the table
# (recovered from the run scripts and the session log), so each dictionary starts with a "backfill"
# provenance block. It is kept so the dictionaries' own origin is reproducible; re-running it rewrites
# the same dictionaries. New runs of the tools write their dictionaries themselves.
#
# Verifications noted below ("reproduces ... byte for byte") were done by re-running the current code
# into a scratch directory and comparing with cmp; the tables in the repo were not touched.
#
# Run from the project root:  bash analyses/lib/backfill/backfill_dictionaries.sh
set -eo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh
conda activate mapping-sims
set -u
cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations
n=0
step() { n=$((n + 1)); }

# ---------------------------------------------------------------------------------------------------
# 1. variant_read_fate.py outputs (s01-s03 trajectory + per-site pileup tables; s05 pileup tables)
# ---------------------------------------------------------------------------------------------------
V=analyses/lib/variant_read_fate.py
for s in s01-baseline-iterative s02-indels s03-heterozygosity; do
    het=""; [ $s = s03-heterozygosity ] && het="--het"
    # aggregate trajectory tables (full chromosome, --sample 0), one job per caller x D
    for g in analyses/$s/results/vrf-agg/*-gtraj.csv; do
        b=$(basename "$g" -gtraj.csv); C=$(echo "$b" | cut -d- -f3); D=$(echo "$b" | cut -d- -f4)
        python $V --dicts-only "cd analyses/$s && sbatch run_vrf_aggregate.sh $C $D" \
            --aggregate $het --sample 0 --divergence "$D" --caller "$C" \
            --out "analyses/$s/results/vrf-agg/vrf-agg-$C-$D"; step
    done
    # per-site pileup ledgers (30,000 sampled SNPs), one job per caller x D
    for f in analyses/$s/results/vrf-pileup/vrf-pileup-*.csv; do
        b=$(basename "$f" .csv); C=$(echo "$b" | cut -d- -f3); D=$(echo "$b" | cut -d- -f4)
        python $V --dicts-only "cd analyses/$s && sbatch run_vrf_pileup.sh $C $D" \
            --pileup $het --sample 30000 --divergence "$D" --caller "$C" --out "$f"; step
    done
done
# s05: one job per caller; field config and mean D are fixed inside each run script
S5=analyses/s05-heterogeneous-divergence
for C in gatk bcftools freebayes; do
    python $V --dicts-only "cd $S5 && sbatch run_vrf_pileup.sh $C" --pileup --sample 50000 \
        --divergence 0.06 --caller $C --out $S5/results/vrf-pileup/vrf-pileup-$C.csv; step
    python $V --dicts-only "cd $S5 && sbatch run_vrf_pileup_capped.sh $C" --pileup --sample 50000 \
        --divergence 0.03 --caller $C --out $S5/results/vrf-pileup-capped/vrf-pileup-$C.csv; step
done

# ---------------------------------------------------------------------------------------------------
# 2. neighborhood_assembly.py (s01, s03: made 2026-08-06 by an earlier version; s02: current version)
# ---------------------------------------------------------------------------------------------------
NB=analyses/lib/neighborhood_assembly.py
OLDNB="made 2026-08-06 by an earlier version of neighborhood_assembly.py that did not yet write the n_neighbor_truth_indel / n_unrecovered_indel columns; the current version reproduces neighborhood-summary.csv and every other per-SNP column byte for byte (checked 2026-10-26)"
nb_dicts() {   # section, truth-data section, extra args...
    local s=$1 t=$2; shift 2
    python $NB --dicts-only "bash analyses/$s/run_neighborhood.sh" \
        --truth "data/$t/sim/30X/19/{D}/0.005/42/mm39-truth.norm.vcf.gz" \
        --recovered "data/$t/iter/{C}/30X/19/{D}/0.005/42/iter{k}/$REC" \
        --pile "analyses/$s/results/vrf-pileup/vrf-pileup-{C}-{D}.csv" \
        --out "analyses/$s/results/neighborhood" --divs 0.06,0.08,0.10 "$@"; step
}
REC=mm39_recovered.vcf.gz
nb_dicts s01-baseline-iterative s01-baseline-iterative --dict-note "$OLDNB"
nb_dicts s03-heterozygosity s01-baseline-iterative --dict-note "$OLDNB"
# s02: keep the existing hand-written per-SNP dictionary (it carries an interpretation section)
REC=mm39_recovered.acc.vcf.gz
nb_dicts s02-indels s02-indels --dicts-only-table summary

# ---------------------------------------------------------------------------------------------------
# 3. tables assembled by run scripts (dictionary-only modes of the scripts themselves)
# ---------------------------------------------------------------------------------------------------
bash analyses/s03-heterozygosity/run_altbias_mechanism.sh --dict-only; step
for D in 0.06 0.10; do
    bash analyses/s04-softclip/run_softclip_toggle.sh $D --dict-only; step
    for K in 2 3; do bash analyses/s04-softclip/run_softclip_toggle_iter.sh $D $K --dict-only; step; done
done
bash analyses/s02-indels/run_calling_metrics.sh --dict-only \
    --note "made 2026-08-13 by an equivalent ad hoc loop over the same three calling_metrics.py commands, ~2 h before this script was written; the current calling_metrics.py reproduces calling-metrics-0.06.csv byte for byte (checked 2026-10-26)"; step

# ---------------------------------------------------------------------------------------------------
# 4. section scripts run with ad hoc commands (commands recovered from the session log)
# ---------------------------------------------------------------------------------------------------
LL=$S5/scripts/local_landscape.py
REF5=data/s05-heterogeneous-divergence/ref-regions/mm39-19.fa
CHK5="the current local_landscape.py reproduces both tables byte for byte (checked 2026-10-26)"
python $LL --data-dir data/s05-heterogeneous-divergence --reference $REF5 \
    --configs L20000-S0.0 L2000-S1.0 L20000-S1.0 L200000-S1.0 L20000-S2.0 --out-prefix $S5/results/s05 \
    --dicts-only "python $LL --data-dir data/s05-heterogeneous-divergence --reference $REF5 --configs L20000-S0.0 L2000-S1.0 L20000-S1.0 L200000-S1.0 L20000-S2.0 --out-prefix $S5/results/s05" \
    --dict-note "ad hoc command (2026-08-08), made by a version of the script edited afterwards; $CHK5"; step
python $LL --data-dir data/s05-heterogeneous-divergence --reference $REF5 --div 0.03 \
    --configs L20000-S0.0 L20000-S1.0-C0.20 --out-prefix $S5/results/s05-realistic \
    --dicts-only "python $LL --data-dir data/s05-heterogeneous-divergence --reference $REF5 --div 0.03 --configs L20000-S0.0 L20000-S1.0-C0.20 --out-prefix $S5/results/s05-realistic" \
    --dict-note "ad hoc command (2026-08-09); $CHK5"; step
SV=analyses/s07-svs/scripts/sv_confounder.py
python $SV --data-dir data/s07-svs --div 0.02 --out analyses/s07-svs/results/s07-sv-confounder.csv \
    --dicts-only "python $SV --data-dir data/s07-svs --div 0.02 --out analyses/s07-svs/results/s07-sv-confounder.csv"; step

# ---------------------------------------------------------------------------------------------------
# 5. (removed) phantom-site outputs: backfilled here on 2026-10-26, then regenerated the same day by
#    the corrected effective_truth.py and run_effective_truth.sh, which write their own dictionaries
#    (re-running the old steps would overwrite them with the old method's text).

# ---------------------------------------------------------------------------------------------------
# 6. orphans: tables whose producing code is no longer (or was never) in the repo as it ran
# ---------------------------------------------------------------------------------------------------
# 6a. codified ad hoc tables: the new scripts reproduce them byte for byte
conda run -n renv Rscript analyses/synthesis/scripts/s01_mapping_vs_calling.R --dict-only; step
bash analyses/s02-indels/scripts/callvcf_composition.sh --dict-only; step

# 6b. (removed) s01-varfate-<D>.csv: backfilled here on 2026-10-26, then regenerated the same day
#     with the current variant_read_fate.py by run_variant_read_fate.sh, which writes its own
#     dictionary (re-running the old step would overwrite it with the old column layout).

# 6c. (removed) pileup_decomp_{summary,bybin}.csv: their dictionaries were backfilled here on
#     2026-10-26; the tables were then deleted the same day (superseded, nothing read them).

echo "wrote dictionaries in $n steps"
