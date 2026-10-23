#!/bin/bash
# Neighborhood assembly analysis for s03 (het query sites; hom-div neighbor load from s01 data).
# Run:  bash run_neighborhood.sh
set -eo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh; conda activate mapping-sims
cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations
python analyses/lib/neighborhood_assembly.py \
    --truth     "data/s01-baseline-iterative/sim/30X/19/{D}/0.005/42/mm39-truth.norm.vcf.gz" \
    --recovered "data/s01-baseline-iterative/iter/{C}/30X/19/{D}/0.005/42/iter{k}/mm39_recovered.vcf.gz" \
    --pile      "analyses/s03-heterozygosity/results/vrf-pileup/vrf-pileup-{C}-{D}.csv" \
    --out       analyses/s03-heterozygosity/results/neighborhood --divs 0.06,0.08,0.10
