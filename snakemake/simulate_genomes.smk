#############################################################################
# Pipeline for generating consensus genomes from reads simulated
# with varying levels of divergence
#############################################################################

import os
import re

#############################################################################
# Example cmd for mouse genome

# snakemake -p -s mapping_sims.smk --configfile simulation-configs/mm39.yaml --profile profiles/slurm_profile/ --cluster-status scripts/slurm_status.py --dryrun

# To generate rulegraph image:
# snakemake -p -s mapping_sims.smk --configfile simulation-configs/mm39.yaml --profile profiles/slurm_profile/ --dryrun --rulegraph | dot -Tpng > dag.png

#############################################################################
# Reference file and path info

ref_str = config["ref_str"];
# The subfolder and abbreviation for the current reference genome

REF = config["ref_file"];
# The full reference genome file

REF_INDEX = config["ref_index"];

outdir = config["sim_outdir"]
# The directory for all output for the current reference genome

#############################################################################
# Simulation parameters

regions = config["regions"];
# Regions to split by for read simulation

covs = config["covs"];
# Average coverage to simulate

divs = config["divs"];
# Average divergence to simulate

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
        expand(os.path.join(outdir, "simulated-genomes", "{cov}X", "{div}", ref_str + "-consensus.fa.amb"), cov=covs, div=divs),
        expand(os.path.join(outdir, "simulated-genomes", "{cov}X", "{div}", ref_str + "-consensus.fa.fai"), cov=covs, div=divs)
        #expand(multiext(os.path.join(outdir, "simulated-genomes", "{cov}X", "{div}", ref_str + "-consensus.fa"), ".amb", ".ann", ".bwt", ".pac", ".sa"), cov=covs, div=divs),
#############################################################################
# Pipeline rules

rule extract_region_ref:
    input:
        ref = REF
    output:
        os.path.join(outdir, "ref-regions", ref_str + "-{region}.fa")
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

def get_runtime_neat(wildcards, attempt):
    times = ["8:00:00", "8:00:00", "16:00:00", "20:00:00", "24:00:00"];
    return times[attempt];
# Sets the run time based on the number of times the job has been restarted

rule simulate_divergent_reads:
    input:
        os.path.join(outdir, "ref-regions", ref_str + "-{region}.fa")
    output:
        read1 = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "divergent", "regions", ref_str + "-{region}_read1.fq.gz"),
        read2 = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "divergent", "regions", ref_str + "-{region}_read2.fq.gz"),
        bam = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "divergent", "regions", ref_str + "-{region}_golden.bam"),
        vcf = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "divergent", "regions", ref_str + "-{region}_golden.vcf.gz"),
    params:
        prefix = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "divergent", "regions", ref_str + "-{region}"),
        div = "{div}",
        cov = "{cov}"
    log:
        os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "divergent", "logs", ref_str + "-{region}.log")
    resources:
        mem = "12g",
        time = get_runtime_neat
    shell:
        """
        python /n/home07/gthomas/projects/Mapping-simulations/pkgs/NEAT/gen_reads.py -r {input} -o {params.prefix} --bam --vcf -R 150 --pe 300 30 -c {params.cov} -M {params.div} &> {log}
        """
# Simulate reads per chromosome with varying levels of divergence (and possibly varying coverage) with NEAT
# This generates read pairs, a golden VCF and a golden BAM for each chromosome
# https://github.com/ncsa/NEAT

#################

def get_runtime_merge_fq(wildcards, attempt):
    times = ["8:00:00", "8:00:00", "12:00:00", "16:00:00", "24:00:00"];
    return times[attempt];
# Sets the run time based on the number of times the job has been restarted

rule merge_simulated_divergent_fastqs:
    input:
        read1 = expand(os.path.join(outdir, "simulated-reads", "{{cov}}X", "{{div}}", "divergent", "regions", ref_str + "-{region}_read1.fq.gz"), region=regions),
        read2 = expand(os.path.join(outdir, "simulated-reads", "{{cov}}X", "{{div}}", "divergent", "regions", ref_str + "-{region}_read2.fq.gz"), region=regions)
    output:
        read1 = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "divergent", ref_str + "_read1.fq.gz"),
        read2 = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "divergent", ref_str + "_read2.fq.gz"),
    params:
        read1 = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "divergent", ref_str + "_read1.fq"),
        read2 = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "divergent", ref_str + "_read2.fq"),
    resources:
        mem = "4g",
        time = get_runtime_merge_fq
    shell:
        """
        zcat {input.read1} > {params.read1}
        gzip {params.read1}
        zcat {input.read2} > {params.read2}
        gzip {params.read2}
        """
# Merge simulated FASTQ files by chromosome and read pair

#################

rule merge_divergent_golden_bams:
    input:
        expand(os.path.join(outdir, "simulated-reads", "{{cov}}X", "{{div}}", "divergent", "regions", ref_str + "-{region}_golden.bam"), region=regions)
    output:
        os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "divergent", ref_str + "_golden.bam"),
    log:
        os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "divergent", "logs", ref_str + "_golden-bam-merge.log"),
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

rule index_divergent_golden_bams:
    input:
        os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "divergent", ref_str + "_golden.bam")
    output:
        os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "divergent", ref_str + "_golden.bam.bai")
    log:
        os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "divergent", "logs", ref_str + "-{cov}X-{div}-index.log")
    resources:
        time = "4:00:00",
        mem = "4g"
    shell:
        """
        samtools index {input} 2> {log}
        """
# Index the mapped reads

#################

rule index_divergent_golden_vcfs_region:
    input:
        os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "divergent", "regions", ref_str + "-{region}_golden.vcf.gz")
    output:
        os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "divergent", "regions", ref_str + "-{region}_golden.vcf.gz.tbi")
    log:
        os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "divergent", "logs", ref_str + "-{region}-tabix.log")
    resources:
        mem = "2g",
        time = "2:00:00"
    shell:
        """
        tabix {input} &> {log}
        """
# Index the golden VCFs, otherwise bcftools concat errors

#################

rule merge_divergent_golden_vcfs:
    input:
        vcf = expand(os.path.join(outdir, "simulated-reads", "{{cov}}X", "{{div}}", "divergent", "regions", ref_str + "-{region}_golden.vcf.gz"), region=regions),
        index = expand(os.path.join(outdir, "simulated-reads", "{{cov}}X", "{{div}}", "divergent", "regions", ref_str + "-{region}_golden.vcf.gz.tbi"), region=regions)
    output:
        os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "divergent", ref_str + "_golden.vcf.gz"),
    log:
        os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "divergent", "logs", ref_str + "_golden-vcf-merge.log"),
    resources:
        cpus = 8,
        time = "2:00:00",
        mem = "4g"
    shell:
        """
        bcftools concat --threads {resources.cpus} -Oz -o {output} {input.vcf} &> {log}
        """
# Merge simulated VCF files for all chromosomes

#################

rule index_divergent_golden_vcfs_merged:
    input:
        os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "divergent", ref_str + "_golden.vcf.gz"),
    output:
        os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "divergent", ref_str + "_golden.vcf.gz.tbi"),
    log:
        os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "divergent", "logs", ref_str + "tabix.log")
    resources:
        mem = "2g",
        time = "2:00:00"
    shell:
        """
        tabix {input} &> {log}
        """
# Index the merged golden VCFs

#################

rule generate_divergent_consensus:
    input:
        vcf = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "divergent", ref_str + "_golden.vcf.gz"),
        vcf_index = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "divergent", ref_str + "_golden.vcf.gz.tbi"),
        ref = REF
    output:
        fasta = os.path.join(outdir, "simulated-genomes", "{cov}X", "{div}", ref_str + "-consensus.fa"),
        chain = os.path.join(outdir, "simulated-genomes", "{cov}X", "{div}", ref_str + "-consensus.chain")
    log:
        os.path.join(outdir, "simulated-genomes", "{cov}X", "{div}", "logs", ref_str + "-consensus.log")
    shell:
        """
        bcftools consensus -f {input.ref} -o {output.fasta} -c {output.chain} {input.vcf} &> {log}
        """
# Generate a consensus by injecting all simulated variants into the reference genome to simulate a diverged reference
#################

rule index_divergent_consensus:
    input:
        os.path.join(outdir, "simulated-genomes", "{cov}X", "{div}", ref_str + "-consensus.fa")
    output:
        bwa_indices = multiext(os.path.join(outdir, "simulated-genomes", "{cov}X", "{div}", ref_str + "-consensus.fa"), ".amb", ".ann", ".bwt", ".pac", ".sa"),
        samtools_index = os.path.join(outdir, "simulated-genomes", "{cov}X", "{div}", ref_str + "-consensus.fa.fai"),
        picard_dict = os.path.join(outdir, "simulated-genomes", "{cov}X", "{div}", ref_str + "-consensus.dict")
    log:
        os.path.join(outdir, "simulated-genomes", "{cov}X", "{div}", "logs", ref_str + "-index.log")
    shell:
        """
        bwa index {input}
        samtools faidx {input}
        rm -f {output.picard_dict}
        picard CreateSequenceDictionary R={input} O={output.picard_dict}
        """
# Index the simulated genome

#############################################################################