#!/bin/bash
# Submit the s04 soft-clip toggle for ITERATIONS 2 and 3 (iteration 1 is done by run_softclip_toggle.sh).
# Confirms empirically that the iter-1 result (soft-clipping is not gatk's limiter) holds as iteration
# converges the consensus and clipping falls. Launch:  bash submit_softclip_toggle_iter.sh
set -euo pipefail
cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations/analyses/s04-softclip

for D in 0.06 0.10; do
    for K in 2 3; do
        BAM=../../data/s01-baseline-iterative/iter/gatk/30X/19/$D/0.005/42/iter$K/mm39_mapped.bam
        if [ ! -s "$BAM" ]; then echo "SKIP D=$D iter=$K: BAM missing"; continue; fi
        jid=$(sbatch --parsable run_softclip_toggle_iter.sh "$D" "$K")
        echo "submitted D=$D iter=$K -> job $jid"
    done
done
