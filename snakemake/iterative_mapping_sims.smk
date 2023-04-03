#############################################################################
# Pipeline for read mapping simulations with varying divergence
#############################################################################

import os
import re

#############################################################################
# Example cmd for mouse genome

# snakemake -p -s iterative_mapping_sims.smk --configfile ../simulation-configs/mm39-iterative.yaml --profile profiles/slurm_profile/

# To generate rulegraph image:
# snakemake -p -s iterative_mapping_sims.smk --configfile ../simulation-configs/mm39-iterative.yaml --profile profiles/slurm_profile/ --rulegraph | dot -Tpng > iterative-mapping-dag.png

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
        expand(os.path.join(outdir, "consensus", "gatk", "{cov}X", "{div}", "iter" + str(num_iters), ref_str + "-{cov}X-{div}d-{het}h-snps-consensus.fa"), cov=covs, div=divs, het=hets),
        expand(os.path.join(outdir, "consensus", "gatk", "{cov}X", "{div}", "iter" + str(num_iters), ref_str + "-{cov}X-{div}d-{het}h-snps-consensus.chain"), cov=covs, div=divs, het=hets),
        # Expected output from generate_consensus

        expand(os.path.join(outdir, "consensus", "gatk", "{cov}X", "{div}", "iter{n}", ref_str + "-{cov}X-{div}d-{het}h-snps-consensus-to-ref.vcf"), cov=covs, div=divs, n=list(range(1,num_iters+1)), het=hets),
        # Expected output from align_to_vcf

        expand(os.path.join(outdir, "summary-files", "{cov}X", ref_str + "-{cov}X-{het}h-" + str(num_iters) + "i-vcf-summary.csv"), cov=covs, het=hets),
        expand(os.path.join(outdir, "summary-files", "{cov}X", ref_str + "-{cov}X-{het}h-" + str(num_iters) + "i-snps.csv.gz"), cov=covs, het=hets),
        expand(os.path.join(outdir, "summary-files", "{cov}X", ref_str + "-{cov}X-{het}h-" + str(num_iters) + "i-bam-summary.csv"), cov=covs, het=hets),
        # Expected output from combine summaries 

## The final expected outputs should be listed in this rule. Only necessary to list final output from final rule, but I found it useful to list them 
## for all rules for debugging (can comment out outputs for rules you don't want to run), though there's also probably a better way to do this

#############################################################################
# Pipeline rules

rule extract_region_sim:
    input:
        ref = os.path.join(indir, "simulated-genomes", "{cov}X", "{div}", ref_str + "-consensus.fa")
    output:
        os.path.join(indir, "simulated-genomes", "{cov}X", "{div}", "regions", ref_str + "-{region}.fa")
    params:
        region = "{region}"
    resources:
        time = "2:00:00",
        mem = "2g"
    shell:
        """
        samtools faidx {input.ref} {params.region} > {output}
        """
# Extracting regions to simulate reads in

#################

# rule index_region_sim:
#     input:
#         os.path.join(indir, "simulated-genomes", "{cov}X", "{div}", "regions", ref_str + "-{region}.fa")
#     output:
#         #bwa_indices = multiext(os.path.join(outdir, "simulated-genomes", "{cov}X", "{div}", ref_str + "-consensus.fa"), ".amb", ".ann", ".bwt", ".pac", ".sa"),
#         #samtools_index = os.path.join(outdir, "simulated-genomes", "{cov}X", "{div}", ref_str + "-consensus.fa.fai"),
#         picard_dict = os.path.join(indir, "simulated-genomes", "{cov}X", "{div}", "regions", ref_str + "-{region}.dict")
#     log:
#         os.path.join(outdir, "simulated-genomes", "{cov}X", "{div}", "regions", "logs", ref_str + "{region}-dict.log")
#     shell:
#         """
#         rm -f {output.picard_dict}
#         picard CreateSequenceDictionary R={input} O={output.picard_dict}
#         """
# Index the simulated genome
## NOT USED

#################

def get_runtime_neat(wildcards, attempt):
    max_time = 168;
    base_time = 6;
    multiplier = 100;
    runtime = base_time + ( float(wildcards.div) * attempt * multiplier );
    if runtime > max_time:
        runtime = max_time;
    return str(round(runtime)) + ":00:00";
# Sets the run time based on the number of times the job has been restarted

rule simulate_reads:
    input:
        os.path.join(indir, "simulated-genomes", "{cov}X", "{div}", "regions", ref_str + "-{region}.fa")
    output:
        read1 = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", "regions", ref_str + "-{region}_read1.fq.gz"),
        read2 = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", "regions", ref_str + "-{region}_read2.fq.gz"),
        bam = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", "regions", ref_str + "-{region}_golden.bam"),
        vcf = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", "regions", ref_str + "-{region}_golden.vcf.gz"),
    params:
        prefix = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", "regions", ref_str + "-{region}"),
        div = "{div}",
        cov = "{cov}"
    log:
        os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", "logs", ref_str + "-{region}.log")
    resources:
        mem = "12g",
        time = get_runtime_neat
    retries: 3
    shell:
        """
        python ../pkgs/NEAT/gen_reads.py -r {input} -o {params.prefix} --bam --vcf -R 150 --pe 300 30 -c {params.cov} -M {params.div} &> {log}
        """
# Simulate reads per chromosome with varying levels of divergence (and possibly varying coverage) with NEAT
# This generates read pairs, a golden VCF and a golden BAM for each chromosome
# https://github.com/ncsa/NEAT

#################

def get_runtime_merge_fq(wildcards, attempt):
    max_time = 168;
    base_time = 6;
    multiplier = 100;
    runtime = base_time + ( float(wildcards.div) * attempt * multiplier );
    if runtime > max_time:
        runtime = max_time;
    return str(round(runtime)) + ":00:00";
# Sets the run time based on the number of times the job has been restarted

rule merge_simulated_fastqs:
    input:
        read1 = expand(os.path.join(outdir, "simulated-reads", "{{cov}}X", "{{div}}", "heterozygous", "{{het}}", "regions", ref_str + "-{region}_read1.fq.gz"), region=regions),
        read2 = expand(os.path.join(outdir, "simulated-reads", "{{cov}}X", "{{div}}", "heterozygous", "{{het}}", "regions", ref_str + "-{region}_read2.fq.gz"), region=regions)
    output:
        read1 = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", ref_str + "_read1.fq.gz"),
        read2 = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", ref_str + "_read2.fq.gz"),
    params:
        read1 = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", ref_str + "_read1.fq"),
        read2 = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", ref_str + "_read2.fq"),
    resources:
        mem = "4g",
        time = get_runtime_merge_fq
    retries: 3    
    shell:
        """
        zcat {input.read1} > {params.read1}
        gzip {params.read1}
        zcat {input.read2} > {params.read2}
        gzip {params.read2}
        """
# Merge simulated FASTQ files by chromosome and read pair

#################

rule merge_golden_bams:
    input:
        expand(os.path.join(outdir, "simulated-reads", "{{cov}}X", "{{div}}", "heterozygous", "{{het}}", "regions", ref_str + "-{region}_golden.bam"), region=regions)
    output:
        os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", ref_str + "_golden.bam"),
    log:
        os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", "logs", ref_str + "_golden-bam-merge.log"),
    resources:
        cpus = 8,
        time = "2:00:00",
        mem = "4g"
    shell:
        """
        samtools merge -@ {resources.cpus} -o {output} {input} &> {log}
        """
# Merge simulated BAM files for all chromosomes

#################

rule index_golden_bams:
    input:
        os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", ref_str + "_golden.bam")
    output:
        os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", ref_str + "_golden.bam.bai")
    log:
        os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", "logs", ref_str + "-{cov}X-{div}-index.log")
    resources:
        time = "4:00:00",
        mem = "4g"
    shell:
        """
        samtools index {input} 2> {log}
        """
# Index the mapped reads

#################

# rule index_golden_vcfs_region:
#     input:
#         os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "regions", ref_str + "-{region}_golden.vcf.gz")
#     output:
#         os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "regions", ref_str + "-{region}_golden.vcf.gz.tbi")
#     log:
#         os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "logs", ref_str + "-{region}-tabix.log")
#     resources:
#         mem = "2g",
#         time = "2:00:00"
#     shell:
#         """
#         tabix {input} &> {log}
#         """
# # Index the golden VCFs, otherwise bcftools concat errors

# #################

# rule merge_golden_vcfs:
#     input:
#         vcf = expand(os.path.join(outdir, "simulated-reads", "{{cov}}X", "{{div}}", "regions", ref_str + "-{region}_golden.vcf.gz"), region=regions),
#         index = expand(os.path.join(outdir, "simulated-reads", "{{cov}}X", "{{div}}", "regions", ref_str + "-{region}_golden.vcf.gz.tbi"), region=regions)
#     output:
#         os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", ref_str + "_golden.vcf.gz"),
#     log:
#         os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "logs", ref_str + "_golden-vcf-merge.log"),
#     resources:
#         cpus = 8,
#         time = "2:00:00",
#         mem = "4g"
#     shell:
#         """
#         bcftools concat --threads {resources.cpus} -Oz -o {output} {input.vcf} &> {log}
#         """
# # Merge simulated VCF files for all chromosomes

# #################

# rule index_golden_vcfs_merged:
#     input:
#         os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", ref_str + "_golden.vcf.gz"),
#     output:
#         os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", ref_str + "_golden.vcf.gz.tbi"),
#     log:
#         os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "logs", ref_str + "tabix.log")
#     resources:
#         mem = "2g",
#         time = "2:00:00"
#     shell:
#         """
#         tabix {input} &> {log}
#         """
# # Index the merged golden VCFs

#################

rule trim_simulated_reads:
    input:
        read1 = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", ref_str + "_read1.fq.gz"),
        read2 = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", ref_str + "_read2.fq.gz")
    output:
        trim_read1 = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", ref_str + "_read1.fastp.fq.gz"),
        trim_read2 = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", ref_str + "_read2.fastp.fq.gz"),
        html = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", ref_str + "_fastp.html"),
        json = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", ref_str + "_fastp.json"),
    resources:
        cpus = 4,
        time = "4:00:00"
    log:
        os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", "logs", ref_str + "-fastp.log")
    shell:
        """
        fastp -i {input.read1} -I {input.read2} --out1 {output.trim_read1} --out2 {output.trim_read2} --detect_adapter_for_pe --cut_front --cut_front_window_size 5 --cut_front_mean_quality 20 -l 25 -j {output.json} -h {output.html} -w {resources.cpus} 2> {log}
        """

# QC on simulated reads
#################

def map_iters_ref(wcs):
    n = int(wcs.n)
    if n == 1:
        r = REF
    elif n > 1:
        r = os.path.join(outdir, "consensus", "gatk", "{cov}X", "{div}", "iter" + str(n-1), ref_str + "-{cov}X-{div}d-{het}h-snps-consensus.fa")
    else:
        raise ValueError("loop numbers must be 1 or greater: received %s" % wcs.n)
    #print(r);
    return r

def map_iters_ind(wcs):
    n = int(wcs.n)
    if n == 1:
        r = multiext(REF, ".amb", ".ann", ".bwt", ".pac", ".sa")
    elif n > 1:
        r = multiext(os.path.join(outdir, "consensus", "gatk", "{cov}X", "{div}", "iter" + str(n-1), ref_str + "-{cov}X-{div}d-{het}h-snps-consensus.fa"), ".amb", ".ann", ".bwt", ".pac", ".sa")
    else:
        raise ValueError("loop numbers must be 1 or greater: received %s" % wcs.n)
    #print(r);
    return r

rule map_simulated_reads:
    input:
        read1 = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", ref_str + "_read1.fastp.fq.gz"),
        read2 = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", ref_str + "_read2.fastp.fq.gz"),
        ref = map_iters_ref,
        ref_index = map_iters_ind
    output:
        os.path.join(outdir, "mapped-reads", "{cov}X", "{div}", "{het}", "iter{n}", ref_str + "-{cov}X-{div}d-{het}h.bam")
    params:
        read_group = "@RG\\tID:" + ref_str + "-{cov}X-{div}d-{het}h\\tLB:NEAT\\tPL:ILLUMINA" + "\\tSM:" + ref_str + "-{cov}X-{div}d-{het}h"
    log:
        os.path.join(outdir, "mapped-reads", "{cov}X", "{div}", "{het}", "iter{n}", "logs", ref_str + "-{cov}X-{div}d-{het}h-iter{n}.log")
    wildcard_constraints:
        n="[0-9]+"
    resources:
        mem= "24g",
        cpus = 12,
        time = "16:00:00"
    shell:
        """
        bwa mem -t {resources.cpus} -R '{params.read_group}' {input.ref} {input.read1} {input.read2} 2> {log} | samtools view -b - 2>> {log} | samtools sort - -o {output} 2>> {log}
        """
# Map the simulated reads with BWA to the full reference genome

#################

rule index_bams:
    input:
        os.path.join(outdir, "mapped-reads", "{cov}X", "{div}", "{het}", "iter{n}", ref_str + "-{cov}X-{div}d-{het}h.bam")
    output:
        os.path.join(outdir, "mapped-reads", "{cov}X", "{div}", "{het}", "iter{n}", ref_str + "-{cov}X-{div}d-{het}h.bam.bai")
    log:
        os.path.join(outdir, "mapped-reads", "{cov}X", "{div}", "{het}", "iter{n}", "logs", ref_str + "-{cov}X-{div}d-{het}h-iter{n}-index.log")
    resources:
        time = "4:00:00",
        mem = "4g"
    shell:
        """
        samtools index {input} 2> {log}
        """
# Index the mapped reads

#################

def get_runtime_gatk_hc(wildcards, attempt):
    max_time = 168;
    base_time = 160;
    multiplier = 200;
    runtime = base_time + ( float(wildcards.div) * attempt * multiplier );
    if runtime > max_time:
        runtime = max_time;
    return str(round(runtime)) + ":00:00";
# Sets the run time based on the number of times the job has been restarted

rule gatk_haplotypecaller: 
    input:
        bam = os.path.join(outdir, "mapped-reads", "{cov}X", "{div}", "{het}", "iter{n}", ref_str + "-{cov}X-{div}d-{het}h.bam"),
        bam_index = os.path.join(outdir, "mapped-reads", "{cov}X", "{div}", "{het}", "iter{n}", ref_str + "-{cov}X-{div}d-{het}h.bam.bai"),
        ref = REF
    output:
        os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "iter{n}", "regions", ref_str + "-{region}-{cov}X-{div}d-{het}h.gvcf.gz")
    params:
        region = "{region}"
    log:
        os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "logs", ref_str + "-haplotypecaller-{region}-{cov}X-{div}d-{het}h-iter{n}.log"),
    resources:
        cpus = 4,
        mem="12g",
        time = get_runtime_gatk_hc
    retries: 3
    shell:
        """
        gatk HaplotypeCaller -R {input.ref} -I {input.bam} -stand-call-conf 30 --native-pair-hmm-threads {resources.cpus} -L {params.region} -ERC GVCF -O {output} &> {log}
        """
# Call variants with GATK HaplotypeCaller by input region

#################

def get_runtime_gatk_gt(wildcards, attempt):
    max_time = 168;
    base_time = 6;
    multiplier = 100;
    runtime = base_time + ( float(wildcards.div) * attempt * multiplier );
    if runtime > max_time:
        runtime = max_time;
    return str(round(runtime)) + ":00:00";
# Sets the run time based on the number of times the job has been restarted

rule gatk_genotypegvcfs: 
    input:
        gvcf = os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "iter{n}", "regions", ref_str + "-{region}-{cov}X-{div}d-{het}h.gvcf.gz"),
        ref = REF
    output:
        os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "iter{n}", "regions", ref_str + "-{region}-{cov}X-{div}d-{het}h.vcf.gz")
    params:
        region = "{region}"
    log:
        os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "logs", ref_str + "-genotypegvcfs-{region}-{cov}X-{div}d-{het}h-iter{n}.log"),
    resources:
        mem="4g",
        time = get_runtime_gatk_gt
    retries: 3
    shell:
        """
        gatk GenotypeGVCFs -R {input.ref} -V {input.gvcf} -O {output} --include-non-variant-sites &> {log}
        """
# Genotype all sites with GATK GenotypeGVCF by input region

#################

rule index_vcfs_regions:
    input:
        os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "iter{n}", "regions", ref_str + "-{region}-{cov}X-{div}d-{het}h.vcf.gz")
    output:
        os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "iter{n}", "regions", ref_str + "-{region}-{cov}X-{div}d-{het}h.vcf.gz.tbi")
    log:
        os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "logs", ref_str + "-{region}-{cov}X-{div}d-{het}h-iter{n}-tabix.log")
    resources:
        mem = "2g",
        time = "2:00:00"
    shell:
        """
        tabix {input} &> {log}
        """

#################

rule merge_vcfs:
    input:
        vcf = expand(os.path.join(outdir, "called-variants", "gatk", "{{cov}}X", "{{div}}", "iter{{n}}", "regions", ref_str + "-{region}-{{cov}}X-{{div}}d-{{het}}h.vcf.gz"), region=regions),
        index = expand(os.path.join(outdir, "called-variants", "gatk", "{{cov}}X", "{{div}}", "iter{{n}}", "regions", ref_str + "-{region}-{{cov}}X-{{div}}d-{{het}}h.vcf.gz.tbi"), region=regions)
    output:
        os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "iter{n}", ref_str + "-{cov}X-{div}d-{het}h.vcf.gz")
    log:
        os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "logs", ref_str + "-{cov}X-{div}d-{het}h-iter{n}-vcf-merge.log"),
    resources:
        cpus = 8,
        time = "2:00:00",
        mem = "4g"
    shell:
        """
        bcftools concat --threads {resources.cpus} -Oz -o {output} {input.vcf} &> {log}
        """
# Merge called VCF files for all chromosomes

#################

rule index_vcfs_merged:
    input:
        os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "iter{n}", ref_str + "-{cov}X-{div}d-{het}h.vcf.gz")
    output:
        os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "iter{n}", ref_str + "-{cov}X-{div}d-{het}h.vcf.gz.tbi")
    log:
        os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "logs", ref_str + "-{cov}X-{div}d-{het}h-iter{n}-tabix.log")
    resources:
        mem = "2g",
        time = "2:00:00"
    shell:
        """
        tabix {input} &> {log}
        """
# Index the merged VCFs

#################

rule compare_vcfs:
    input:
        golden = os.path.join(indir, "simulated-reads", "{cov}X", "{div}", "divergent", "regions", ref_str + "-{region}_golden.vcf.gz"),
        golden_index = os.path.join(indir, "simulated-reads", "{cov}X", "{div}", "divergent", "regions", ref_str + "-{region}_golden.vcf.gz.tbi"), 
        query = os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "iter{n}", "regions", ref_str + "-{region}-{cov}X-{div}d-{het}h.vcf.gz"),
        query_index = os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "iter{n}", "regions", ref_str + "-{region}-{cov}X-{div}d-{het}h.vcf.gz.tbi")
    params:
        outdir = os.path.join(outdir, "summary-files", "{cov}X", "{div}", "regions"),
        div = "{div}",
        het = "{het}",
        cov = "{cov}",
        iteration = "{n}"
    output:
        summary = os.path.join(outdir, "summary-files", "{cov}X", "{div}", "regions", ref_str + "-iter{n}-{region}-{cov}X-{div}d-{het}h-compare-vcf-summary.csv"),
        snps = os.path.join(outdir, "summary-files", "{cov}X", "{div}", "regions", ref_str + "-iter{n}-{region}-{cov}X-{div}d-{het}h-compare-vcf-snps.csv")
    resources:
        mem = "4g",
        time = "1:00:00"
    shell:
        """
        python /n/home07/gthomas/projects/Mapping-simulations/scripts/compare_vcfs_2.py {params.cov} {params.div} {params.het} {params.iteration} {input.golden} {input.query} {params.outdir} {output.summary} {output.snps}
        """
# Run the compare_vcfs script to get number of variants compared to golden file
# Use to combine files:
# for f in *.csv; do (cat "${f}"; echo); done | grep -v "#" | sort -r | uniq

#################

rule liftover_bams:
    input:
        bam = os.path.join(outdir, "mapped-reads", "{cov}X", "{div}","{het}", "iter{n}", ref_str + "-{cov}X-{div}d-{het}h.bam"),
        chain = os.path.join(indir, "simulated-genomes", "{cov}X", "{div}", ref_str + "-consensus.chain")
    output:
        bam = os.path.join(outdir, "mapped-reads", "{cov}X", "{div}", "{het}", "iter{n}", ref_str + "-{cov}X-{div}d-{het}h-crossmap.bam"),
        bam_sorted = os.path.join(outdir, "mapped-reads", "{cov}X", "{div}", "{het}", "iter{n}", ref_str + "-{cov}X-{div}d-{het}h-crossmap.sorted.bam"),
        bam_index = os.path.join(outdir, "mapped-reads", "{cov}X", "{div}", "{het}", "iter{n}", ref_str + "-{cov}X-{div}d-{het}h-crossmap.sorted.bam.bai")
    params:
        out_prefix = os.path.join(outdir, "mapped-reads", "{cov}X", "{div}","{het}", "iter{n}", ref_str + "-{cov}X-{div}d-{het}h-crossmap"),
    log:
        os.path.join(outdir, "mapped-reads", "{cov}X", "{div}", "{het}", "iter{n}", "logs", ref_str + "-{cov}X-{div}d-{het}h-crossmap.log")
    resources:
        mem = "24g"
    shell:
        """
        CrossMap.py bam -a {input.chain} {input.bam} {params.out_prefix} &> {log}
        """
# Liftover the coordinates of the mapped reads from the mouse reference to the simulated divergent genome (same as golden bam here)
# https://crossmap.sourceforge.net/#convert-bam-cram-sam-format-files 

#################

rule compare_bams:
    input:
        golden = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", ref_str + "_golden.bam"),
        golden_index = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", ref_str + "_golden.bam.bai"),
        query = os.path.join(outdir, "mapped-reads", "{cov}X", "{div}", "{het}", "iter{n}", ref_str + "-{cov}X-{div}d-{het}h-crossmap.sorted.bam"),
        query_index = os.path.join(outdir, "mapped-reads", "{cov}X", "{div}", "{het}", "iter{n}", ref_str + "-{cov}X-{div}d-{het}h-crossmap.sorted.bam.bai"),
        ref_index = REF_INDEX
    params:
       div = "{div}",
       het = "{het}",
       cov = "{cov}",
       reg = ",".join(regions),
       iteration = "{n}"
    output:
        summary = os.path.join(outdir, "summary-files", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h-iter{n}-compare-bam-summary.csv"),
    resources:
        mem = "24g",
        time = "8:00:00"
    shell:
        """
        python /n/home07/gthomas/projects/Mapping-simulations/scripts/compare_bams.py {params.reg} {params.cov} {params.div} {params.het} {params.iteration} {input.ref_index} {input.golden} {input.query} {output.summary}
        """
# Run the compare_vcfs script to get number of variants compared to golden file
# Use to combine files:
# for f in *.csv; do (cat "${f}"; echo); done | grep -v "#" | sort -r | uniq

#################

# def get_runtime_qualimap(wildcards, attempt):
#     times = ["8:00:00", "12:00:00", "16:00:00", "20:00:00", "32:00:00"];
#     return times[attempt];
# # Sets the run time based on the number of times the job has been restarted
# rule qualimap_bams:
#     input:
#         golden_bam = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", ref_str + "_golden.bam"),
#         mapped_bam = os.path.join(outdir, "mapped-reads", "{cov}X", "{div}", "iter{n}", ref_str + "-{cov}X-{div}.bam")
#     output:
#         golden_qualimap = os.path.join(outdir, "qualimap", "{cov}X", "{div}", "iter{n}", "golden", "qualimapReport.html"),
#         mapped_qualimap = os.path.join(outdir, "qualimap", "{cov}X", "{div}", "iter{n}", "mapped", "qualimapReport.html")
#     params:
#         golden_outdir = os.path.join(outdir, "qualimap", "{cov}X", "{div}", "iter{n}", "golden"),
#         mapped_outdir = os.path.join(outdir, "qualimap", "{cov}X", "{div}", "iter{n}", "mapped")
#     log:
#         golden_log = os.path.join(outdir, "qualimap", "{cov}X", "{div}", "logs", ref_str + "-iter{n}-golden.log"),
#         mapped_log = os.path.join(outdir, "qualimap", "{cov}X", "{div}", "logs", ref_str + "-iter{n}-mapped.log")
#     resources:
#         cpus = 32,
#         time = get_runtime_qualimap
#     shell:
#         """
#         qualimap bamqc -bam {input.golden_bam} -nt {resources.cpus} -outdir {params.golden_outdir} -outformat html --java-mem-size=12G &> {log.golden_log}
#         qualimap bamqc -bam {input.mapped_bam} -nt {resources.cpus} -outdir {params.mapped_outdir} -outformat html --java-mem-size=12G &> {log.mapped_log}
#         """        

#################

rule filter_vcf:
    input:
        vcf = os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "iter{n}", ref_str + "-{cov}X-{div}d-{het}h.vcf.gz"),
        vcf_index = os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "iter{n}", ref_str + "-{cov}X-{div}d-{het}h.vcf.gz.tbi")
    output:
        os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "iter{n}", ref_str + "-{cov}X-{div}d-{het}h-filtered.vcf.gz")
    params:
        filt = FILTER_STR
    log:
        os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "logs", ref_str + "-iter{n}-{cov}X-{div}d-{het}h-vcf-filter.log")
    shell:
        """
        bcftools filter -m+ -e {params.filt} -s pseudoit --IndelGap 5 -Oz -o {output} {input.vcf} &> {log}
        """
#rule filter_snps
# NOTE: bcftools leaves an empty file if it errors, but snakemake thinks this has completed when it tries again....

#################

rule index_vcfs_filtered:
    input:
        os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "iter{n}", ref_str + "-{cov}X-{div}d-{het}h-filtered.vcf.gz")
    output:
        os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "iter{n}", ref_str + "-{cov}X-{div}d-{het}h-filtered.vcf.gz.tbi")
    log:
        os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "logs", ref_str + "-iter{n}-{cov}X-{div}d-{het}h-filter-tabix.log")
    resources:
        mem = "2g",
        time = "2:00:00"
    shell:
        """
        tabix {input} &> {log}
        """
# Index the merged VCFs

#################

rule select_snps:
    input:
        vcf = os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "iter{n}", ref_str + "-{cov}X-{div}d-{het}h-filtered.vcf.gz"),
        vcf_index = os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "iter{n}", ref_str + "-{cov}X-{div}d-{het}h-filtered.vcf.gz.tbi")
    output:
        os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "iter{n}", ref_str + "-{cov}X-{div}d-{het}h-filtered-snps.vcf.gz"),
    log:
        os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "logs", ref_str + "-iter{n}-{cov}X-{div}d-{het}h-vcf-select-snps.log")
    shell:
        """
        gatk SelectVariants -V {input.vcf} -O {output} -select-type SNP -xl-select-type INDEL -xl-select-type MIXED -xl-select-type SYMBOLIC &> {log}
        """
# Select SNPs only

#################

rule index_vcfs_snps:
    input:
        os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "iter{n}", ref_str + "-{cov}X-{div}d-{het}h-filtered-snps.vcf.gz")
    output:
        os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "iter{n}", ref_str + "-{cov}X-{div}d-{het}h-filtered-snps.vcf.gz.tbi")
    log:
        os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "logs", ref_str + "-iter{n}-{cov}X-{div}d-{het}h-snps-tabix.log")
    resources:
        mem = "2g",
        time = "2:00:00"
    shell:
        """
        tabix {input} &> {log}
        """
# Index the SNP VCFs
#################

# rule get_mask_sites:
#     input:
#         os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "iter{n}", ref_str + "-{cov}X-{div}d-{het}h-filtered-snps.vcf.gz")
#     output:
#         os.path.join(outdir, "consensus", "gatk", "{cov}X", "{div}", "iter" + str(num_iters), ref_str + "-{cov}X-{div}d-{het}h-masksites.bed")
#     shell:
#         """
#         zgrep \"\./\.\" {input} | awk '{{OFS=\"\t\"; if ($0 !~ /\#/); print $1, $2-1, $2}}' | bedtools merge -i - > {output}
#         """

# rule mask_fa:
#     input:
#         ref = os.path.join(outdir, "consensus", "gatk", "{cov}X", "{div}", "{het}", "iter" + str(num_iters-1), ref_str + "-snps-consensus.fa"),
#         bed = os.path.join(outdir, "consensus", "gatk", "{cov}X", "{div}", "iter" + str(num_iters), ref_str + "-{cov}X-{div}d-{het}h-masksites.bed")
#     output:
#         os.path.join(outdir, "consensus", "gatk", "{cov}X", "{div}", "iter" + str(num_iters), ref_str + "-snps-consensus-iter" + str(num_iters-1) + "-{cov}X-{div}d-{het}h-masked.fa")
#     shell:
#         """
#         bedtools maskfasta -fi {input.ref} -bed {input.bed} -soft -fo {output}
#         """

#################

rule generate_consensus:
    input:
        vcf = os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "iter{n}", ref_str + "-{cov}X-{div}d-{het}h-filtered-snps.vcf.gz"),
        vcf_index = os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "iter{n}", ref_str + "-{cov}X-{div}d-{het}h-filtered-snps.vcf.gz.tbi"),
        ref = REF
    output:
        fasta = os.path.join(outdir, "consensus", "gatk", "{cov}X", "{div}", "iter{n}", ref_str + "-{cov}X-{div}d-{het}h-snps-consensus.fa"),
        chain = os.path.join(outdir, "consensus", "gatk", "{cov}X", "{div}", "iter{n}", ref_str + "-{cov}X-{div}d-{het}h-snps-consensus.chain")
    log:
        os.path.join(outdir, "consensus", "gatk", "{cov}X", "{div}", "logs", ref_str + "-iter{n}-{cov}X-{div}d-{het}h-consensus.log")
    shell:
        """
        bcftools consensus -f {input.ref} -o {output.fasta} -c {output.chain} -e "FILTER='pseudoit' || FILTER='IndelGap'" {input.vcf} &> {log}
        """   

#################

rule index_consensus:
    input:
        os.path.join(outdir, "consensus", "gatk", "{cov}X", "{div}", "iter{n}", ref_str + "-{cov}X-{div}d-{het}h-snps-consensus.fa")
    output:
        samtools_index = os.path.join(outdir, "consensus", "gatk", "{cov}X", "{div}", "iter{n}", ref_str + "-{cov}X-{div}d-{het}h-snps-consensus.fa.fai"),
        bwa_index = multiext(os.path.join(outdir, "consensus", "gatk", "{cov}X", "{div}", "iter{n}", ref_str + "-{cov}X-{div}d-{het}h-snps-consensus.fa"), ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        os.path.join(outdir, "consensus", "gatk", "{cov}X", "{div}", "logs", ref_str + "-iter{n}-{cov}X-{div}d-{het}h-consensus-index.log")
    shell:
        """
        samtools faidx {input} 2> {log}
        bwa index {input} 2>> {log}
        """   

#################

rule align_consensus:
    input:
        ref = REF,
        consensus = os.path.join(outdir, "consensus", "gatk", "{cov}X", "{div}", "iter{n}", ref_str + "-{cov}X-{div}d-{het}h-snps-consensus.fa"),
        consensus_index = os.path.join(outdir, "consensus", "gatk", "{cov}X", "{div}", "iter{n}", ref_str + "-{cov}X-{div}d-{het}h-snps-consensus.fa.fai")
    output:
        paf = os.path.join(outdir, "consensus", "gatk", "{cov}X", "{div}", "iter{n}", ref_str + "-{cov}X-{div}d-{het}h-snps-consensus-to-ref.paf")
    log:
        os.path.join(outdir, "consensus", "gatk", "{cov}X", "{div}", "logs", ref_str + "-iter{n}-{cov}X-{div}d-{het}h-minimap.log")
    resources:
        cpus = 16,
        mem = "100g"
    shell:
        """
        minimap2 -x asm20 --cs -K 200M -I 20G -t {resources.cpus} {input.ref} {input.consensus} -o {output.paf} > {log}
        """

#################

rule paf_to_vcf:
    input:
        ref = REF,
        paf = os.path.join(outdir, "consensus", "gatk", "{cov}X", "{div}", "iter{n}", ref_str + "-{cov}X-{div}d-{het}h-snps-consensus-to-ref.paf")
    output:
        os.path.join(outdir, "consensus", "gatk", "{cov}X", "{div}", "iter{n}", ref_str + "-{cov}X-{div}d-{het}h-snps-consensus-to-ref.vcf")
    resources:
        time='96:00:00'
    shell:
        """
        sort -k6,6 -k8,8n {input.paf} | paftools.js call -f {input.ref} - > {output}
        """

#################

rule combine_summaries:
    input:
        vcf_summary = expand(os.path.join(outdir, "summary-files", "{cov}X", "{div}", "regions", ref_str + "-iter{n}-{region}-{cov}X-{div}d-{het}h-compare-vcf-summary.csv"), cov=covs, div=divs, region=regions, het=hets, n=range(1, num_iters+1)),
        snp_summary = expand(os.path.join(outdir, "summary-files", "{cov}X", "{div}", "regions", ref_str + "-iter{n}-{region}-{cov}X-{div}d-{het}h-compare-vcf-snps.csv"), cov=covs, div=divs, region=regions, het=hets, n=range(1, num_iters+1)),
        bam_summary = expand(os.path.join(outdir, "summary-files", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h-iter{n}-compare-bam-summary.csv"), cov=covs, div=divs, het=hets, n=range(1, num_iters+1)) 
    output:
        vcf_out = os.path.join(outdir, "summary-files", "{cov}X", ref_str + "-{cov}X-{het}h-" + str(num_iters) + "i-vcf-summary.csv"),
        snp_out = os.path.join(outdir, "summary-files", "{cov}X", ref_str + "-{cov}X-{het}h-" + str(num_iters) + "i-snps.csv.gz"),
        bam_out = os.path.join(outdir, "summary-files", "{cov}X", ref_str + "-{cov}X-{het}h-" + str(num_iters) + "i-bam-summary.csv")
    params:
        indir = os.path.join(outdir, "summary-files", "{cov}X"),
        snp_int = os.path.join(outdir, "summary-files", "{cov}X", ref_str + "-{cov}X-{het}h-" + str(num_iters) + "i-snps.csv")
    resources:
        mem = "4g",
        time = "1:00:00"
    shell:
        """
        find {params.indir} -type f -name *-compare-vcf-summary.csv -exec awk 1 {{}} \; | grep -v "#" | sort -r | uniq > {output.vcf_out}
        find {params.indir} -type f -name *-compare-vcf-snps.csv -exec awk 1 {{}} \; > {params.snp_int}
        gzip {params.snp_int}
        find {params.indir} -type f -name *-compare-bam-summary.csv -exec awk 1 {{}} \; | grep -v "#" | sort -r | uniq > {output.bam_out}
        """
# Combines the output files from compare_vcfs and compare_bams
# Originally ran separately as:
# find . -type f -name *-vcf-summary.csv -not -name mm39-20X-vcf-summary.csv -exec awk 1 {} \; | grep -v "#" | sort -r | uniq > ../mm39-20X-vcf-summary.csv
# find . -type f -name *-vcf-snps.csv -not -name mm39-20X-vcf-snps.csv -exec awk 1 {} \; > ../mm39-20X-vcf-snps.csv
# find . -type f -name *-bam-summary.csv -exec awk 1 {} \; | grep -v "#" | sort -r | uniq > ../mm39-20X-bam-summary.csv


#############################################################################
