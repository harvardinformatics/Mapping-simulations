#!/bin/bash
# Submit the s01 companion-decomposition jobs (per-SNP read fate + position-centric pileup).
# 3 callers x 3 divergences (0.06/0.08/0.10 -- the range with meaningful not-called->called flips).
# Launch:  bash submit_vrf_pileup.sh
set -euo pipefail
cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations/analyses/s01-baseline-iterative

for C in gatk bcftools freebayes; do
    for D in 0.06 0.08 0.10; do
        TB=../../data/varfate/$D/truth.oa.bam
        if [ ! -s "$TB" ]; then echo "SKIP $C $D: truth BAM missing"; continue; fi
        jid=$(sbatch --parsable run_vrf_pileup.sh "$C" "$D")
        echo "submitted caller=$C div=$D -> job $jid"
    done
done
