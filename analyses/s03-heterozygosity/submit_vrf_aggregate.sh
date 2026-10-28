#!/bin/bash
# Submit the s03 full-chromosome aggregate HET read-fate jobs (3 callers x 6 divergences) to SLURM.
# Uses the s01 (SNP-only) data + cached oH truth BAMs, matching the s01/s02 divergence set.
# Launch:  bash submit_vrf_aggregate.sh
set -euo pipefail
cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations/analyses/s03-heterozygosity

for C in gatk bcftools freebayes; do
    for D in 0.005 0.02 0.04 0.06 0.08 0.10; do      # 0.005/0.02 were an ad hoc loop until 2026-10-27
        TB=../../data/varfate/$D/truth.oa.bam
        if [ ! -s "$TB" ]; then echo "ERROR $C $D: truth BAM missing ($TB)"; exit 1; fi
        jid=$(sbatch --parsable run_vrf_aggregate.sh "$C" "$D")
        echo "submitted caller=$C div=$D -> job $jid"
    done
done
