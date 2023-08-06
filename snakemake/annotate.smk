#############################################################################
# Pipeline for read mapping simulations with varying divergence
#############################################################################

import os
import re

#############################################################################
# Example cmd for mouse genome

# snakemake -p -s annotate.smk --configfile ../simulation-configs/mm39-iterative.yaml --profile profiles/slurm_profile/ --dryrun

# To generate rulegraph image:
# snakemake -p -s annotate.smk --configfile ../simulation-configs/mm39-iterative.yaml --profile profiles/slurm_profile/ --dryrun --rulegraph | dot -Tpng > annotate-dag.png

#############################################################################
# Reference file and path info

indir = config["sim_indir"];

ref_str = config["ref_str"];
# The subfolder and abbreviation for the current reference genome

REF = config["ref_file"];
# The full reference genome file

REF_INDEX = config["ref_index"];

outdir = config["sim_outdir"];
# The directory for all output for the current reference genome

num_iters = config["iterations"];

FILTER_STR = config["filter"]

REF_REPEAT_BED = config["ref_repeat_bed"]
REF_GFF = config["ref_gff"]

#############################################################################
# Simulation parameters

regions = config["regions"];
# Regions to split by for read simulation

covs = config["covs"];
# Average coverage to simulate

divs = config["divs"];
# Average divergence to simulate

hets = config["hets"];
# Average heterozygosity to simulate

var_types = ["tps", "fns"];
map_types = ["exact", "unmapped"];

#############################################################################
# Python checks

if not os.path.isfile(REF):
    raise FileNotFoundError ("\n * ERROR: Cannot find reference genome: " + REF + "\n");
# Check for the reference genome at the provided path

ref_indexes = { "samtools faidx" : [REF + ".fai"], 
                "bwa index" : multiext(REF, ".amb", ".ann", ".bwt", ".pac", ".sa"), 
                "picard CreateSequenceDictionary" : [REF.replace(".fa", ".dict")] };
for index_type in ref_indexes:
    for index_file in ref_indexes[index_type]:
        if not os.path.isfile(index_file):
            raise FileNotFoundError ("\n * ERROR: You are missing a reference genome file to be generated with " + index_type + "\n");
# A check for all the expected reference file indexes

#############################################################################
# Final rule - rule that depends on final expected output file and initiates all
# the other rules

localrules: all

rule all:
    input:
        expand(os.path.join(outdir, "summary-files", "{cov}X", "{div}", "regions", ref_str + "-iter{n}-{region}-{cov}X-{div}d-{het}h-compare-vcf-{var_type}-repeats.bed"), cov=covs, region=regions, div=divs, n=list(range(1,num_iters+1)), het=hets, var_type=var_types),
        expand(os.path.join(outdir, "summary-files", "{cov}X", "{div}", "regions", ref_str + "-iter{n}-{region}-{cov}X-{div}d-{het}h-compare-vcf-{var_type}-genes.bed"), cov=covs, region=regions, div=divs, n=list(range(1,num_iters+1)), het=hets, var_type=var_types),
        
        expand(os.path.join(outdir, "summary-files", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h-iter{n}-compare-bam-{map_type}-repeats.bed"), cov=covs, div=divs, n=list(range(1,num_iters+1)), het=hets, map_type=map_types),
        expand(os.path.join(outdir, "summary-files", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h-iter{n}-compare-bam-{map_type}-genes.bed"), cov=covs, div=divs, n=list(range(1,num_iters+1)), het=hets, map_type=map_types)

#############################################################################
# Pipeline functions

rule variant_overlaps:
    input:
        variant_file = os.path.join(outdir, "summary-files", "{cov}X", "{div}", "regions", ref_str + "-iter{n}-{region}-{cov}X-{div}d-{het}h-compare-vcf-{var_type}.bed"),
        repeat_bed = REF_REPEAT_BED,
        gff = REF_GFF
    output:
        repeat_overlaps = os.path.join(outdir, "summary-files", "{cov}X", "{div}", "regions", ref_str + "-iter{n}-{region}-{cov}X-{div}d-{het}h-compare-vcf-{var_type}-repeats.bed"),
        gene_overlaps = os.path.join(outdir, "summary-files", "{cov}X", "{div}", "regions", ref_str + "-iter{n}-{region}-{cov}X-{div}d-{het}h-compare-vcf-{var_type}-genes.bed")
    shell:
        """
        bedtools intersect -a {input.variant_file} -b {input.repeat_bed} -c > {output.repeat_overlaps}
        bedtools intersect -a {input.variant_file} -b {input.gff} -c > {output.gene_overlaps}
        """

#################

rule read_overlaps:
    input:
        bam_file = os.path.join(outdir, "summary-files", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h-iter{n}-compare-bam-{map_type}.bam"),
        repeat_bed = REF_REPEAT_BED,
        gff = REF_GFF
    output:
        repeat_overlaps = os.path.join(outdir, "summary-files", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h-iter{n}-compare-bam-{map_type}-repeats.bed"),
        gene_overlaps = os.path.join(outdir, "summary-files", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h-iter{n}-compare-bam-{map_type}-genes.bed")
    shell:
        """
        bedtools intersect -a {input.bam_file} -b {input.repeat_bed} -c -bed > {output.repeat_overlaps}
        bedtools intersect -a {input.bam_file} -b {input.gff} -c -bed > {output.gene_overlaps}
        """

#############################################################################