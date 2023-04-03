#############################################################################
# Pipeline for read mapping simulations with varying divergence
#############################################################################

import os
import re

#############################################################################
# Example cmd for mouse genome

# snakemake -p -s iterative_divergence.smk --configfile ../simulation-configs/mm39-iterative.yaml --profile profiles/slurm_profile/ --dryrun

# To generate rulegraph image:
# snakemake -p -s iterative_divergence.smk --configfile ../simulation-configs/mm39-iterative.yaml --profile profiles/slurm_profile/ --dryrun --rulegraph | dot -Tpng > iterative-divergence-dag.png

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
iter_list = [ str(i) for i in range(1,num_iters+1) ];

FILTER_STR = config["filter"];

window_size = 10000;
window_str = "10kb";

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
        expand(os.path.join(outdir, "summary-files", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-golden-snps-" + window_str + "-intersect.bed"), cov=covs, div=divs),
        expand(os.path.join(outdir, "summary-files", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-golden-snps-" + window_str + "-intersect-annotated.tab"), cov=covs, div=divs),
        expand(os.path.join(outdir, "summary-files", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h-iter{n}-snps-" + window_str + "-intersect.bed"), cov=covs, div=divs, het=hets, n=iter_list),
        expand(os.path.join(outdir, "summary-files", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h-iter{n}-snps-" + window_str + "-intersect-annotated.tab"), cov=covs, div=divs, het=hets, n=iter_list),
        expand(os.path.join(outdir, "summary-files", "{cov}X", "window-snps-{cov}X-summary.tab"), cov=covs)
        # Expected output from combine summaries 

## The final expected outputs should be listed in this rule. Only necessary to list final output from final rule, but I found it useful to list them 
## for all rules for debugging (can comment out outputs for rules you don't want to run), though there's also probably a better way to do this

#############################################################################

rule windows:
    input:
        ref_ind = REF_INDEX
    output:
        bed_windows = os.path.join(os.path.dirname(os.path.realpath(REF_INDEX)), ref_str + "-" + window_str + ".bed")
    params:
        winsize = window_size
    shell:
        """
        bedtools makewindows -g {input.ref_ind} -w {params.winsize} > {output.bed_windows}
        """

rule intersect_ref:
    input:
        bed_windows = os.path.join(os.path.dirname(os.path.realpath(REF_INDEX)), ref_str + "-" + window_str + ".bed"),
        sim_vcf = os.path.join(indir, "simulated-reads", "{cov}X", "{div}", "divergent", ref_str + "_golden.vcf.gz"),
    output:
        ref_intersect = os.path.join(outdir, "summary-files", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-golden-snps-" + window_str + "-intersect.bed"),
        ref_intersect_ann = os.path.join(outdir, "summary-files", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-golden-snps-" + window_str + "-intersect-annotated.tab")
    params:
        cov = "{cov}",
        div = "{div}"
    shell:
        """
        bedtools intersect -c -a {input.bed_windows} -b {input.sim_vcf} > {output.ref_intersect}
        awk 'OFS="\\t" {{print $0, {params.cov}, {params.div}, "NA", 0}}' {output.ref_intersect} > {output.ref_intersect_ann}
        """

rule intersect_iters:
    input:
        bed_windows = os.path.join(os.path.dirname(os.path.realpath(REF_INDEX)), ref_str + "-" + window_str + ".bed"),
        iter_vcf = os.path.join(outdir, "consensus", "gatk", "{cov}X", "{div}", "iter{n}", ref_str + "-{cov}X-{div}d-{het}h-snps-consensus-to-ref.vcf")
    output:
        iter_intersect = os.path.join(outdir, "summary-files", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h-iter{n}-snps-" + window_str + "-intersect.bed"),
        iter_intersect_ann = os.path.join(outdir, "summary-files", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h-iter{n}-snps-" + window_str + "-intersect-annotated.tab")
    params:
        cov = "{cov}",
        div = "{div}",
        het = "{het}",
        n = "{n}"
    shell:
        """
        bedtools intersect -c -a {input.bed_windows} -b {input.iter_vcf} > {output.iter_intersect}
        awk 'OFS="\\t" {{print $0, {params.cov}, {params.div}, {params.het}, {params.n}}}' {output.iter_intersect} > {output.iter_intersect_ann}
        """

rule combine_intersects:
    input:
        ref_intersect_ann = expand(os.path.join(outdir, "summary-files", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-golden-snps-" + window_str + "-intersect-annotated.tab"), cov=covs, div=divs),
        iter_intersect_ann = expand(os.path.join(outdir, "summary-files", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h-iter{n}-snps-" + window_str + "-intersect-annotated.tab"), cov=covs, div=divs, het=hets, n=iter_list)
    output:
        combined_intersect = os.path.join(outdir, "{cov}X", "summary-files", "window-snps-{cov}X-summary.tab")
    shell:
        """
        cat {input.ref_intersect_ann} {input.iter_intersect_ann} > {output}
        """