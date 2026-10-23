#!/bin/bash
# Submit the s02 companion-decomposition jobs (per-SNP read fate + pileup), 3 callers x 3 divergences
# (0.06/0.08/0.10 -- matching the s01 decomposition). Launch:  bash submit_vrf_pileup.sh
set -euo pipefail
cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations/analyses/s02-indels

for C in gatk bcftools freebayes; do
    for D in 0.06 0.08 0.10; do
        TB=../../data/s02-indels/sim/30X/19/$D/0.005/42/mm39-truth.oa.bam
        if [ ! -s "$TB" ]; then echo "SKIP $C $D: truth BAM missing"; continue; fi
        jid=$(sbatch --parsable run_vrf_pileup.sh "$C" "$D")
        echo "submitted caller=$C div=$D -> job $jid"
    done
done
