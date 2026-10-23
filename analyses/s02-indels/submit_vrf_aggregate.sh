#!/bin/bash
# Submit the s02 full-chromosome aggregate SNP read-fate jobs (3 callers x 4 divergences) to SLURM.
# Matches the s01 divergence set (0.04-0.10) for a direct s01-vs-s02 comparison of how the indel-rich
# context changes SNP read fate. Launch:  bash submit_vrf_aggregate.sh
set -euo pipefail
cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations/analyses/s02-indels

for C in gatk bcftools freebayes; do
    for D in 0.04 0.06 0.08 0.10; do
        TB=../../data/s02-indels/sim/30X/19/$D/0.005/42/mm39-truth.oa.bam
        if [ ! -s "$TB" ]; then echo "SKIP $C $D: truth BAM missing"; continue; fi
        jid=$(sbatch --parsable run_vrf_aggregate.sh "$C" "$D")
        echo "submitted caller=$C div=$D -> job $jid"
    done
done
