#!/bin/bash
# Submit the 12 full-chromosome aggregate-ledger jobs (3 callers x 4 divergences) to SLURM.
# Each condition is its own `sbatch run_vrf_aggregate.sh <CALLER> <DIV>` job, so the heavy compute
# runs on the cluster, not the login node. Launch from this directory:  bash submit_vrf_aggregate.sh
set -euo pipefail
cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations/analyses/s01-baseline-iterative

for C in gatk bcftools freebayes; do
    for D in 0.04 0.06 0.08 0.10; do
        TB=../../data/varfate/$D/truth.oa.bam
        if [ ! -s "$TB" ]; then echo "SKIP $C $D: truth BAM missing ($TB)"; continue; fi
        jid=$(sbatch --parsable run_vrf_aggregate.sh "$C" "$D")
        echo "submitted caller=$C div=$D -> job $jid"
    done
done
