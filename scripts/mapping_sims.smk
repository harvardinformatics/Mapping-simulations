#############################################################################
# Pipeline for read mapping simulations with varying divergence
#############################################################################

import os
import re

#############################################################################
# Example cmd

# snakemake -p -s mapping_sims.smk --dryrun --configfile mm39-config.yaml --profile profiles/slurm_profile/

#############################################################################
# Reference file and path info

data_dir = config["data_dir"];
# The base datadir where reference genomes and output data are stored

base_ref_dir = os.path.join(data_dir, "reference-genomes");
# The location of all reference genomes for the project

ref_str = config["ref_str"];
# The subfolder and abbreviation for the current reference genome

ref_dir = os.path.join(base_ref_dir, ref_str);
# The current reference genome directory

REF = os.path.join(ref_dir, config["ref_file"]);
# The full reference genome file

#ref_files = [ f for f in os.listdir(ref_dir) if re.findall(ref_str + "-chr(.*).fa", f) != [] ];
# Get all the individiual chromosome FASTA files for the current reference genome
# Assumes that files are named as <ref_str>-chr<chrome name>.fa

# chromes = [ f.replace(ref_str + "-", "").replace(".fa", "") for f in ref_files ];
# print(chromes);
# Get the chromosome names from the chromosome FASTA files

#ref_files = [ os.path.join(ref_dir, f) for f in ref_files ];
# Prepend the reference directory to each chromosome FASTA file
## This was all for getting the individual chromosome/region fasta files from the reference
## directory -- now handled as regions in config file and split out in the first rule if the
## files don't exist

#############################################################################
# Simulation parameters

regions = config["regions"];
# Regions to split by for read simulation

covs = config["covs"];
# Average coverage to simulate

#divs = ["0.00", "0.02", "0.04", "0.06", "0.08", "0.10", "0.12", "0.14", "0.16", "0.18", "0.20"];
divs = config["divs"];
# Average divergence to simulate

#############################################################################
# Python checks

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
        expand(os.path.join(ref_dir, ref_str + "-{region}.fa"), region=regions),
        # Expected outputs from extract_region

        expand(os.path.join(data_dir, "simulated-reads", ref_str, "{cov}X", "{div}", "regions", ref_str + "-{region}_read1.fq.gz"), cov=covs, div=divs, region=regions),
        expand(os.path.join(data_dir, "simulated-reads", ref_str, "{cov}X", "{div}", "regions", ref_str + "-{region}_read2.fq.gz"), cov=covs, div=divs, region=regions),
        # Expected outputs from simulate_reads

        expand(os.path.join(data_dir, "simulated-reads", ref_str, "{cov}X", "{div}", ref_str + "_read1.fq.gz"), cov=covs, div=divs),
        expand(os.path.join(data_dir, "simulated-reads", ref_str, "{cov}X", "{div}", ref_str + "_read2.fq.gz"), cov=covs, div=divs),
        # Expected outputs from merge_simulated_reads

        expand(os.path.join(data_dir, "simulated-reads", ref_str, "{cov}X", "{div}", ref_str + "_golden.bam"), cov=covs, div=divs),
        # Expected output from merge_golden_bams

        expand(os.path.join(data_dir, "simulated-reads", ref_str, "{cov}X", "{div}", ref_str + "_read1.fastp.fq.gz"), cov=covs, div=divs),
        expand(os.path.join(data_dir, "simulated-reads", ref_str, "{cov}X", "{div}", ref_str + "_read2.fastp.fq.gz"), cov=covs, div=divs),
        # Expected outputs from trim_simulated_reads

        expand(os.path.join(data_dir, "mapped-reads", ref_str, "{cov}X", "{div}", ref_str + "-{cov}X-{div}.bam"), cov=covs, div=divs),
        # Expected outputs from map_simulated_reads

        expand(os.path.join(data_dir, "mapped-reads", ref_str, "{cov}X", "{div}", ref_str + "-{cov}X-{div}.bam.bai"), cov=covs, div=divs),
        # Expected outputs from index_bams

        expand(os.path.join(data_dir, "called-variants", ref_str, "gatk", "{cov}X", "{div}", ref_str + "-{region}-{cov}X-{div}.vcf.gz"), cov=covs, div=divs, region=regions),
        # Expected outputs from gatk_haplotypecaller

        expand(os.path.join(data_dir, "called-variants", ref_str, "gatk", "{cov}X", "{div}", ref_str + "-{region}-{cov}X-{div}.vcf.gz.tbi"), cov=covs, div=divs, region=regions),
        # Expected output from index_vcfs

        expand(os.path.join(data_dir, "qualimap", ref_str, "{cov}X", "{div}", "golden", "qualimapReport.html"), cov=covs, div=divs),
        expand(os.path.join(data_dir, "qualimap", ref_str, "{cov}X", "{div}", "mapped", "qualimapReport.html"), cov=covs, div=divs)
        # Expected output from qualimap_bams

## The final expected outputs should be listed in this rule. Only necessary to list final output from final rule, but I found it useful to list them 
## for all rules for debugging (can comment out outputs for rules you don't want to run), though there's also probably a better way to do this

#############################################################################
# Pipeline rules

rule extract_region:
    input:
        ref = REF
    output:
        os.path.join(ref_dir, ref_str + "-{region}.fa")
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

rule index_fasta:
    input:
        os.path.join(ref_dir, ref_str + "-{region}.fa")
    output:
        os.path.join(ref_dir, ref_str + "-{region}.fa.fai")
    shell:
        """
        samtools faidx {input}
        """
# Indexing the chromosome FASTAs... not needed if merging before mapping.

#################

rule simulate_reads:
    input:
        os.path.join(ref_dir, ref_str + "-{region}.fa")
    output:
        read1 = os.path.join(data_dir, "simulated-reads", ref_str, "{cov}X", "{div}", "regions", ref_str + "-{region}_read1.fq.gz"),
        read2 = os.path.join(data_dir, "simulated-reads", ref_str, "{cov}X", "{div}", "regions", ref_str + "-{region}_read2.fq.gz"),
        bam = os.path.join(data_dir, "simulated-reads", ref_str, "{cov}X", "{div}", "regions", ref_str + "-{region}_golden.bam"),
        vcf = os.path.join(data_dir, "simulated-reads", ref_str, "{cov}X", "{div}", "regions", ref_str + "-{region}_golden.vcf.gz"),
    params:
        prefix = os.path.join(data_dir, "simulated-reads", ref_str, "{cov}X", "{div}", "regions", ref_str + "-{region}"),
        div = "{div}",
        cov = "{cov}"
    log:
        os.path.join(data_dir, "simulated-reads", ref_str, "{cov}X", "{div}", "logs", ref_str + "-{region}.log")
    shell:
        """
        python ../pkgs/NEAT/gen_reads.py -r {input} -o {params.prefix} --bam --vcf -R 150 --pe 300 30 -c {params.cov} -M {params.div} &> {log}
        """
# Simulate reads per chromosome with varying levels of divergence (and possibly varying coverage) with NEAT
# This generates read pairs, a golden VCF and a golden BAM for each chromosome
# https://github.com/ncsa/NEAT

#################

rule merge_simulated_fastqs:
    input:
        read1 = expand(os.path.join(data_dir, "simulated-reads", ref_str, "{{cov}}X", "{{div}}", "regions", ref_str + "-{region}_read1.fq.gz"), region=regions),
        read2 = expand(os.path.join(data_dir, "simulated-reads", ref_str, "{{cov}}X", "{{div}}", "regions", ref_str + "-{region}_read2.fq.gz"), region=regions)
    output:
        read1 = os.path.join(data_dir, "simulated-reads", ref_str, "{cov}X", "{div}", ref_str + "_read1.fq.gz"),
        read2 = os.path.join(data_dir, "simulated-reads", ref_str, "{cov}X", "{div}", ref_str + "_read2.fq.gz"),
    params:
        read1 = os.path.join(data_dir, "simulated-reads", ref_str, "{cov}X", "{div}", ref_str + "_read1.fq"),
        read2 = os.path.join(data_dir, "simulated-reads", ref_str, "{cov}X", "{div}", ref_str + "_read2.fq"),
    resources:
        time = "2:00:00",
        mem = "4g"
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
        expand(os.path.join(data_dir, "simulated-reads", ref_str, "{{cov}}X", "{{div}}", "regions", ref_str + "-{region}_golden.bam"), region=regions)
    output:
        os.path.join(data_dir, "simulated-reads", ref_str, "{cov}X", "{div}", ref_str + "_golden.bam"),
    log:
        os.path.join(data_dir, "simulated-reads", ref_str, "{cov}X", "{div}", "logs", ref_str + "_golden-bam-merge.log"),
    resources:
        cpus = 8,
        time = "2:00:00",
        mem = "4g"
    shell:
        """
        samtools merge -@ {resources.cpus} -o {output} {input} &> {log}
        """
# Merge simulated FASTQ files by chromosome and read pair

#################

rule trim_simulated_reads:
    input:
        read1 = os.path.join(data_dir, "simulated-reads", ref_str, "{cov}X", "{div}", ref_str + "_read1.fq.gz"),
        read2 = os.path.join(data_dir, "simulated-reads", ref_str, "{cov}X", "{div}", ref_str + "_read2.fq.gz")
    output:
        trim_read1 = os.path.join(data_dir, "simulated-reads", ref_str, "{cov}X", "{div}", ref_str + "_read1.fastp.fq.gz"),
        trim_read2 = os.path.join(data_dir, "simulated-reads", ref_str, "{cov}X", "{div}", ref_str + "_read2.fastp.fq.gz"),
        html = os.path.join(data_dir, "simulated-reads", ref_str, "{cov}X", "{div}", ref_str + "_fastp.html"),
        json = os.path.join(data_dir, "simulated-reads", ref_str, "{cov}X", "{div}", ref_str + "_fastp.json"),
    resources:
        cpus = 4,
        time = "2:00:00"
    log:
        os.path.join(data_dir, "simulated-reads", ref_str, "{cov}X", "{div}", "logs", ref_str + "-fastp.log")
    shell:
        """
        fastp -i {input.read1} -I {input.read2} --out1 {output.trim_read1} --out2 {output.trim_read2} --detect_adapter_for_pe --cut_front --cut_front_window_size 5 --cut_front_mean_quality 20 -l 25 -j {output.json} -h {output.html} -w {resources.cpus} 2> {log}
        """

# QC on simulated reads
#################

rule map_simulated_reads:
    input:
        read1 = os.path.join(data_dir, "simulated-reads", ref_str, "{cov}X", "{div}", ref_str + "_read1.fastp.fq.gz"),
        read2 = os.path.join(data_dir, "simulated-reads", ref_str, "{cov}X", "{div}", ref_str + "_read2.fastp.fq.gz"),
        ref = REF,
        ref_index = multiext(REF, ".amb", ".ann", ".bwt", ".pac", ".sa")
    output:
        os.path.join(data_dir, "mapped-reads", ref_str, "{cov}X", "{div}", ref_str + "-{cov}X-{div}.bam")
    params:
        read_group = "@RG\\tID:" + ref_str + "-{cov}X-{div}\\tLB:NEAT\\tPL:ILLUMINA" + "\\tSM:" + ref_str + "-{cov}X-{div}"
    log:
        os.path.join(data_dir, "mapped-reads", ref_str, "{cov}X", "{div}", "logs", ref_str + "-{cov}X-{div}.log")
    resources:
        cpus = 12,
        time = "4:00:00"
    shell:
        """
        bwa mem -t {resources.cpus} -R '{params.read_group}' {input.ref} {input.read1} {input.read2} 2> {log} | samtools view -b - 2>> {log} | samtools sort - -o {output} 2>> {log}
        """
# Map the simulated reads with BWA to the full reference genome

#################

rule index_bams:
    input:
        os.path.join(data_dir, "mapped-reads", ref_str, "{cov}X", "{div}", ref_str + "-{cov}X-{div}.bam")
    output:
        os.path.join(data_dir, "mapped-reads", ref_str, "{cov}X", "{div}", ref_str + "-{cov}X-{div}.bam.bai")
    log:
        os.path.join(data_dir, "mapped-reads", ref_str, "{cov}X", "{div}", "logs", ref_str + "-{cov}X-{div}-index.log")
    resources:
        time = "2:00:00",
        mem = "4g"
    shell:
        """
        samtools index {input} 2> {log}
        """
# Index the mapped reads

#################

def get_runtime(wildcards, attempt):
    times = ["24:00:00", "32:00:00", "40:00:00", "48:00:00", "56:00:00"];
    return times[attempt];
# Sets the run time based on the number of times the job has been restarted

rule gatk_haplotypecaller:
    input:
        bam = os.path.join(data_dir, "mapped-reads", ref_str, "{cov}X", "{div}", ref_str + "-{cov}X-{div}.bam"),
        bam_index = os.path.join(data_dir, "mapped-reads", ref_str, "{cov}X", "{div}", ref_str + "-{cov}X-{div}.bam.bai"),
        ref = REF
    output:
        os.path.join(data_dir, "called-variants", ref_str, "gatk", "{cov}X", "{div}", ref_str + "-{region}-{cov}X-{div}.vcf.gz")
    params:
        region = "{region}"
    log:
        os.path.join(data_dir, "called-variants", ref_str, "gatk", "{cov}X", "{div}", "logs", ref_str + "-{region}-{cov}X-{div}.log"),
    resources:
        cpus = 16,
        time = get_runtime
    shell:
        """
        gatk HaplotypeCaller -R {input.ref} -I {input.bam} -stand-call-conf 30 --native-pair-hmm-threads 8 -L {params.region} -O {output} &> {log}
        """
# Call variants with GATK HaplotypeCaller by input region

#################

rule index_vcfs:
    input:
        os.path.join(data_dir, "called-variants", ref_str, "gatk", "{cov}X", "{div}", ref_str + "-{region}-{cov}X-{div}.vcf.gz")
    output:
        os.path.join(data_dir, "called-variants", ref_str, "gatk", "{cov}X", "{div}", ref_str + "-{region}-{cov}X-{div}.vcf.gz.tbi")
    log:
        os.path.join(data_dir, "called-variants", ref_str, "gatk", "{cov}X", "{div}", "logs", ref_str + "-{region}-{cov}X-{div}-tabix.log")
    resources:
        mem = "2g",
        time = "2:00:00"
    shell:
        """
        tabix {input} &> {log}
        """

#################

rule qualimap_bams:
    input:
        golden_bam = os.path.join(data_dir, "simulated-reads", ref_str, "{cov}X", "{div}", ref_str + "_golden.bam"),
        mapped_bam = os.path.join(data_dir, "mapped-reads", ref_str, "{cov}X", "{div}", ref_str + "-{cov}X-{div}.bam")
    output:
        golden_qualimap = os.path.join(data_dir, "qualimap", ref_str, "{cov}X", "{div}", "golden", "qualimapReport.html"),
        mapped_qualimap = os.path.join(data_dir, "qualimap", ref_str, "{cov}X", "{div}", "mapped", "qualimapReport.html")
    params:
        golden_outdir = os.path.join(data_dir, "qualimap", ref_str, "{cov}X", "{div}", "golden"),
        mapped_outdir = os.path.join(data_dir, "qualimap", ref_str, "{cov}X", "{div}", "mapped")
    log:
        golden_log = os.path.join(data_dir, "qualimap", ref_str, "{cov}X", "{div}", "logs", ref_str + "_golden.log"),
        mapped_log = os.path.join(data_dir, "qualimap", ref_str, "{cov}X", "{div}", "logs", ref_str + "_mapped.log")
    resources:
        cpus = 8,
        time = "2:00:00"
    shell:
        """
        qualimap bamqc -bam {input.golden_bam} -nt {resources.cpus} -outdir {params.golden_outdir} -outformat html --java-mem-size=4G &> {log.golden_log}
        qualimap bamqc -bam {input.mapped_bam} -nt {resources.cpus} -outdir {params.mapped_outdir} -outformat html --java-mem-size=4G &> {log.mapped_log}
        """        

#################
#############################################################################