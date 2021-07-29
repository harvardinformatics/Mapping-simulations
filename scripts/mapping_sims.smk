#############################################################################
# Pipeline for read mapping simulations with increasing divergence
#############################################################################

import sys
import os
import re

#############################################################################
# Example cmd

# snakemake -p -s mapping_sims.smk --dryrun --configfile mm39-config.yaml --profile profiles/slurm_profile/

#############################################################################
# Setup - These variables will need to change for each user/genome

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

base_sim_dir = os.path.join(data_dir, "simulated-reads");

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
# Final rule - rule that depends on final expected output file and initiates all
# the other rules

localrules: all

rule all:
    input:
        expand(os.path.join(ref_dir, ref_str + "-{region}.fa"), region=regions),
        # Expected outputs from extract_region

        expand(os.path.join(data_dir, "simulated-reads", ref_str, "{cov}X", "{div}", ref_str + "-{region}_read1.fq.gz"), cov=covs, div=divs, region=regions),
        expand(os.path.join(data_dir, "simulated-reads", ref_str, "{cov}X", "{div}", ref_str + "-{region}_read2.fq.gz"), cov=covs, div=divs, region=regions),
        # Expected outputs from simulate_reads

        expand(os.path.join(data_dir, "simulated-reads", ref_str, "{cov}X", "{div}", ref_str + "_read1.fq.gz"), cov=covs, div=divs),
        expand(os.path.join(data_dir, "simulated-reads", ref_str, "{cov}X", "{div}", ref_str + "_read2.fq.gz"), cov=covs, div=divs),
        # Extracted outputs from merge_simulated_reads

        expand(os.path.join(data_dir, "mapped-reads", ref_str, "{cov}X", "{div}", ref_str + "-{cov}X-{div}.bam"), cov=covs, div=divs)
        # Extracted outputs from map_simulated_reads

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
        time = "1:00:00",
        mem = "2g"
    shell:
        """
        samtools faidx {input.ref} {params.region} > {output}
        """
# Indexing the chromosome FASTAs... not needed if merging before mapping.

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

print(expand(os.path.join(ref_dir, ref_str + "-{region}.fa"), region=regions))

rule simulate_reads:
    input:
        os.path.join(ref_dir, ref_str + "-{region}.fa")
    output:
        read1 = os.path.join(data_dir, "simulated-reads", ref_str, "{cov}X", "{div}", ref_str + "-{region}_read1.fq.gz"),
        read2 = os.path.join(data_dir, "simulated-reads", ref_str, "{cov}X", "{div}", ref_str + "-{region}_read2.fq.gz")
    params:
        prefix = os.path.join(data_dir, "simulated-reads", ref_str, "{cov}X", "{div}", ref_str + "-{region}"),
        div = "{div}",
        cov = "{cov}"
    log:
        os.path.join(data_dir, "simulated-reads", ref_str, "{cov}X", "{div}", "logs", ref_str + "-{region}.log")
    resources:
        time = "6:00:00"
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
        read1 = expand(os.path.join(data_dir, "simulated-reads", ref_str, "{{cov}}X", "{{div}}", ref_str + "-{region}_read1.fq.gz"), region=regions),
        read2 = expand(os.path.join(data_dir, "simulated-reads", ref_str, "{{cov}}X", "{{div}}", ref_str + "-{region}_read2.fq.gz"), region=regions)
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

rule map_simulated_reads:
    input:
        read1 = os.path.join(data_dir, "simulated-reads", ref_str, "{cov}X", "{div}", ref_str + "_read1.fq.gz"),
        read2 = os.path.join(data_dir, "simulated-reads", ref_str, "{cov}X", "{div}", ref_str + "_read2.fq.gz"),
        ref = REF,
        ref_index = REF + ".fai"
    output:
        os.path.join(data_dir, "mapped-reads", ref_str, "{cov}X", "{div}", ref_str + "-{cov}X-{div}.bam")
    log:
        os.path.join(data_dir, "mapped-reads", ref_str, "{cov}X", "{div}", "logs", ref_str + "-{cov}X-{div}.log")
    resources:
        cpus = 8
    shell:
        """
        bwa mem -t {resources.cpus} {input.ref} {input.read1} {input.read2} 2> {log} | samtools view -b - 2>> {log} | samtools sort - -o {output} 2>> {log}
        """
# Map the simulated reads with BWA to the full reference genome

#################

## TODO: 
# 1. pre-mapping qc?
# 2. variant calling
# 3. fix name in profile

#############################################################################