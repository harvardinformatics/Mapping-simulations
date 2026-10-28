#!/bin/bash
# Submit the s03 HET companion-decomposition jobs (per-variant read fate + pileup), 3 callers x 3
# divergences (0.06/0.08/0.10). Launch:  bash submit_vrf_pileup.sh
set -euo pipefail
cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations/analyses/s03-heterozygosity

for C in gatk bcftools freebayes; do
    for D in 0.06 0.08 0.10; do
        TB=../../data/varfate/$D/truth.oa.bam
        if [ ! -s "$TB" ]; then echo "ERROR $C $D: truth BAM missing ($TB)"; exit 1; fi
        jid=$(sbatch --parsable run_vrf_pileup.sh "$C" "$D")
        echo "submitted caller=$C div=$D -> job $jid"
    done
done
