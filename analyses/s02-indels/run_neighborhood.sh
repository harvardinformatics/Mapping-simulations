#!/bin/bash
# Neighborhood assembly analysis for s02 (SNPs in indel context). Run:  bash run_neighborhood.sh
set -eo pipefail
source /n/home07/gthomas/miniconda3/etc/profile.d/conda.sh; conda activate mapping-sims
cd /n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations
python analyses/lib/neighborhood_assembly.py \
    --truth     "data/s02-indels/sim/30X/19/{D}/0.005/42/mm39-truth.norm.vcf.gz" \
    --recovered "data/s02-indels/iter/{C}/30X/19/{D}/0.005/42/iter{k}/mm39_recovered.acc.vcf.gz" \
    --pile      "analyses/s02-indels/results/vrf-pileup/vrf-pileup-{C}-{D}.csv" \
    --out       analyses/s02-indels/results/neighborhood --divs 0.06,0.08,0.10
