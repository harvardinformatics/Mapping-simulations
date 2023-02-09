#############################################################################
# Pipeline for read mapping simulations with varying divergence
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
# Average divergence simulated in genome

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
        expand(os.path.join(outdir, "summary-files", "{cov}X", "{div}", "regions", ref_str + "-{region}-{cov}X-{div}d-{het}h-compare-vcf-summary.csv"), cov=covs, div=divs, het=hets, region=regions),
        expand(os.path.join(outdir, "summary-files", "{cov}X", "{div}", "regions", ref_str + "-{region}-{cov}X-{div}d-{het}h-compare-vcf-snps.csv"), cov=covs, div=divs, het=hets, region=regions),
        # Expected output from compare_vcfs

        expand(os.path.join(outdir, "summary-files", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h-compare-bam-summary.csv"), cov=covs, div=divs, het=hets),
        # Expected output from compare_bams

        #expand(os.path.join(outdir, "qualimap", "{cov}X", "{div}", "golden", "qualimapReport.html"), cov=covs, div=divs),
        #expand(os.path.join(outdir, "qualimap", "{cov}X", "{div}", "mapped", "qualimapReport.html"), cov=covs, div=divs)
        # Expected output from qualimap_bams

        #expand(os.path.join(outdir, "ref-regions", ref_str + "-{region}.fa"), region=regions),
        # Expected outputs from extract_region

        #expand(os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "regions", ref_str + "-{region}_read1.fq.gz"), cov=covs, div=divs, region=regions),
        #expand(os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "regions", ref_str + "-{region}_read2.fq.gz"), cov=covs, div=divs, region=regions),
        # Expected outputs from simulate_reads

        #expand(os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", ref_str + "_read1.fq.gz"), cov=covs, div=divs),
        #expand(os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", ref_str + "_read2.fq.gz"), cov=covs, div=divs),
        # Expected outputs from merge_simulated_reads

        #expand(os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", ref_str + "_golden.bam"), cov=covs, div=divs),
        # Expected output from merge_golden_bams

        #expand(os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", ref_str + "_golden.bam.bai"), cov=covs, div=divs),
        # Expected output from index_golden_bams

        #expand(os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "regions", ref_str + "-{region}_golden.vcf.gz.tbi"), cov=covs, div=divs, region=regions),
        # Expected output from index_golden_vcfs_regions

        #expand(os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", ref_str + "_golden.vcf.gz"), cov=covs, div=divs),
        # Expected output from merge_golden_vcfs

        #expand(os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", ref_str + "_golden.vcf.gz.tbi"), cov=covs, div=divs),
        # Expected output from index_golden_vcfs_merged

        #expand(os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", ref_str + "_read1.fastp.fq.gz"), cov=covs, div=divs),
        #expand(os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", ref_str + "_read2.fastp.fq.gz"), cov=covs, div=divs),
        # Expected outputs from trim_simulated_reads

        #expand(os.path.join(outdir, "mapped-reads", "{cov}X", "{div}", ref_str + "-{cov}X-{div}.bam"), cov=covs, div=divs),
        # Expected outputs from map_simulated_reads

        #expand(os.path.join(outdir, "mapped-reads", "{cov}X", "{div}", ref_str + "-{cov}X-{div}.bam.bai"), cov=covs, div=divs),
        # Expected outputs from index_bams

        #expand(os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "regions", ref_str + "-{region}-{cov}X-{div}.gvcf.gz"), cov=covs, div=divs, region=regions),
        # Expected outputs from gatk_haplotypecaller

        #expand(os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "regions", ref_str + "-{region}-{cov}X-{div}.vcf.gz"), cov=covs, div=divs, region=regions),
        # Expected outputs from gatk_genotypegvcfs

        #expand(os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "regions", ref_str + "-{region}-{cov}X-{div}.vcf.gz.tbi"), cov=covs, div=divs, region=regions),
        # Expected output from index_vcfs_regions

        #expand(os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", ref_str + ".vcf.gz"), cov=covs, div=divs),
        # Expected output from merge_vcfs

        #expand(os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", ref_str + ".vcf.gz.tbi"), cov=covs, div=divs),
        # Expected output from index_vcfs_merged

## The final expected outputs should be listed in this rule. Only necessary to list final output from final rule, but I found it useful to list them 
## for all rules for debugging (can comment out outputs for rules you don't want to run), though there's also probably a better way to do this

#############################################################################
# Pipeline rules

rule extract_region_sim:
    input:
        ref = os.path.join(outdir, "simulated-genomes", "{cov}X", "{div}", ref_str + "-consensus.fa")
    output:
        os.path.join(outdir, "simulated-genomes", "{cov}X", "{div}", "regions", ref_str + "-{region}.fa")
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

rule index_region_sim:
    input:
        os.path.join(outdir, "simulated-genomes", "{cov}X", "{div}", "regions", ref_str + "-{region}.fa")
    output:
        #bwa_indices = multiext(os.path.join(outdir, "simulated-genomes", "{cov}X", "{div}", ref_str + "-consensus.fa"), ".amb", ".ann", ".bwt", ".pac", ".sa"),
        #samtools_index = os.path.join(outdir, "simulated-genomes", "{cov}X", "{div}", ref_str + "-consensus.fa.fai"),
        picard_dict = os.path.join(outdir, "simulated-genomes", "{cov}X", "{div}", "regions", ref_str + "-{region}.dict")
    log:
        os.path.join(outdir, "simulated-genomes", "{cov}X", "{div}", "regions", "logs", ref_str + "{region}-dict.log")
    shell:
        """
        rm -f {output.picard_dict}
        picard CreateSequenceDictionary R={input} O={output.picard_dict}
        """
# Index the simulated genome

######################################################################################################

def get_runtime_neat(wildcards, attempt):
    times = ["8:00:00", "8:00:00", "16:00:00", "20:00:00", "24:00:00"];
    return times[attempt];
# Sets the run time based on the number of times the job has been restarted

rule simulate_heterozygous_reads:
    input:
        os.path.join(outdir, "simulated-genomes", "{cov}X", "{div}", "regions", ref_str + "-{region}.fa")
    output:
        read1 = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", "regions", ref_str + "-{region}_read1.fq.gz"),
        read2 = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", "regions", ref_str + "-{region}_read2.fq.gz"),
        bam = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", "regions", ref_str + "-{region}_golden.bam"),
        vcf = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", "regions", ref_str + "-{region}_golden.vcf.gz"),
    params:
        prefix = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", "regions", ref_str + "-{region}"),
        het = "{het}",
        cov = "{cov}"
    log:
        os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", "logs", ref_str + "-{region}.log")
    resources:
        mem = "12g",
        time = get_runtime_neat
    shell:
        """
        python /n/home07/gthomas/projects/Mapping-simulations/pkgs/NEAT/gen_reads.py -r {input} -o {params.prefix} --bam --vcf -R 150 --pe 300 30 -c {params.cov} -M {params.het} &> {log}
        """
# Simulate reads per chromosome with varying levels of divergence (and possibly varying coverage) with NEAT
# This generates read pairs, a golden VCF and a golden BAM for each chromosome
# https://github.com/ncsa/NEAT

#################

def get_runtime_merge_fq(wildcards, attempt):
    times = ["8:00:00", "8:00:00", "12:00:00", "16:00:00", "24:00:00"];
    return times[attempt];
# Sets the run time based on the number of times the job has been restarted

rule merge_simulated_heterozygous_fastqs:
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

rule index_golden_vcfs_region:
    input:
        os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", "regions", ref_str + "-{region}_golden.vcf.gz")
    output:
        os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", "regions", ref_str + "-{region}_golden.vcf.gz.tbi")
    log:
        os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", "logs", ref_str + "-{region}-tabix.log")
    resources:
        mem = "2g",
        time = "2:00:00"
    shell:
        """
        tabix {input} &> {log}
        """
# Index the golden VCFs, otherwise bcftools concat errors

# Merge simulated FASTQ files by chromosome and read pair
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

rule map_simulated_reads:
    input:
        read1 = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", ref_str + "_read1.fastp.fq.gz"),
        read2 = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", ref_str + "_read2.fastp.fq.gz"),
        ref = REF,
        ref_index = multiext(REF, ".amb", ".ann", ".bwt", ".pac", ".sa")
    output:
        os.path.join(outdir, "mapped-reads", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h.bam")
    params:
        read_group = "@RG\\tID:" + ref_str + "-{cov}X-{div}d-{het}h\\tLB:NEAT\\tPL:ILLUMINA" + "\\tSM:" + ref_str + "-{cov}X-{div}d-{het}h"
    log:
        os.path.join(outdir, "mapped-reads", "{cov}X", "{div}", "logs", ref_str + "-{cov}X-{div}d-{het}h.log")
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
        os.path.join(outdir, "mapped-reads", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h.bam")
    output:
        os.path.join(outdir, "mapped-reads", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h.bam.bai")
    log:
        os.path.join(outdir, "mapped-reads", "{cov}X", "{div}", "logs", ref_str + "-{cov}X-{div}d-{het}h-index.log")
    resources:
        time = "4:00:00",
        mem = "4g"
    shell:
        """
        samtools index {input} 2> {log}
        """
# Index the mapped reads

#################

def get_runtime_gatk(wildcards, attempt):
    times = ["48:00:00", "56:00:00", "72:00:00", "96:00:00", "120:00:00"];
    return times[attempt];
# Sets the run time based on the number of times the job has been restarted

rule gatk_haplotypecaller: 
    input:
        bam = os.path.join(outdir, "mapped-reads", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h.bam"),
        bam_index = os.path.join(outdir, "mapped-reads", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h.bam.bai"),
        ref = REF
    output:
        os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "regions", ref_str + "-{region}-{cov}X-{div}d-{het}h.gvcf.gz")
    params:
        region = "{region}"
    log:
        os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "logs", ref_str + "-haplotypecaller-{region}-{cov}X-{div}d-{het}h.log"),
    resources:
        cpus = 4,
        mem="4g",
        time = get_runtime_gatk
    shell:
        """
        gatk HaplotypeCaller -R {input.ref} -I {input.bam} -stand-call-conf 30 --native-pair-hmm-threads {resources.cpus} -L {params.region} -ERC GVCF -O {output} &> {log}
        """
# Call variants with GATK HaplotypeCaller by input region

#################

def get_runtime_gatk(wildcards, attempt):
    times = ["4:00:00", "6:00:00", "12:00:00", "24:00:00", "32:00:00"];
    return times[attempt];
# Sets the run time based on the number of times the job has been restarted

rule gatk_genotypegvcfs: 
    input:
        gvcf = os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "regions", ref_str + "-{region}-{cov}X-{div}d-{het}h.gvcf.gz"),
        ref = REF
    output:
        os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "regions", ref_str + "-{region}-{cov}X-{div}d-{het}h.vcf.gz")
    params:
        region = "{region}"
    log:
        os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "logs", ref_str + "-genotypegvcfs-{region}-{cov}X-{div}d-{het}h.log"),
    resources:
        mem="4g",
        time = get_runtime_gatk
    shell:
        """
        gatk GenotypeGVCFs -R {input.ref} -V {input.gvcf} -O {output} --include-non-variant-sites &> {log}
        """
# Genotype all sites with GATK GenotypeGVCF by input region

#################

rule index_vcfs_regions:
    input:
        os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "regions", ref_str + "-{region}-{cov}X-{div}d-{het}h.vcf.gz")
    output:
        os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "regions", ref_str + "-{region}-{cov}X-{div}d-{het}h.vcf.gz.tbi")
    log:
        os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "logs", ref_str + "-{region}-{cov}X-{div}d-{het}h-tabix.log")
    resources:
        mem = "2g",
        time = "2:00:00"
    shell:
        """
        tabix {input} &> {log}
        """
# Index the called VCFs by region
## NOT USED

#################

rule merge_vcfs:
    input:
        vcf = expand(os.path.join(outdir, "called-variants", "gatk", "{{cov}}X", "{{div}}", "regions", ref_str + "-{region}-{{cov}}X-{{div}}d-{{het}}h.vcf.gz"), region=regions),
        index = expand(os.path.join(outdir, "called-variants", "gatk", "{{cov}}X", "{{div}}", "regions", ref_str + "-{region}-{{cov}}X-{{div}}d-{{het}}h.vcf.gz.tbi"), region=regions)
    output:
        os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h.vcf.gz"),
    log:
        os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "logs", ref_str + "-{cov}X-{div}d-{het}h-vcf-merge.log"),
    resources:
        cpus = 8,
        time = "2:00:00",
        mem = "4g"
    shell:
        """
        bcftools concat --threads {resources.cpus} -Oz -o {output} {input.vcf} &> {log}
        """
# Merge called VCF files for all chromosomes
## NOT USED

#################

rule index_vcfs_merged:
    input:
        os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h.vcf.gz")
    output:
        os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h.vcf.gz.tbi")
    log:
        os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "logs", ref_str + "-{cov}X-{div}d-{het}h-tabix.log")
    resources:
        mem = "2g",
        time = "2:00:00"
    shell:
        """
        tabix {input} &> {log}
        """
# Index the merged VCFs
## NOT USED

######################################################################################################

rule liftover_vcfs:
    input:
        vcf = os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "regions", ref_str + "-{region}-{cov}X-{div}d-{het}h.vcf.gz"),
        chain = os.path.join(outdir, "simulated-genomes", "{cov}X", "{div}", ref_str + "-consensus.chain"),
        sim_ref = os.path.join(outdir, "simulated-genomes", "{cov}X", "{div}", "regions", ref_str + "-{region}.fa"),
        sim_dict = os.path.join(outdir, "simulated-genomes", "{cov}X", "{div}", "regions", ref_str + "-{region}.dict")
    output:
        liftover_vcf = os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "regions", ref_str + "-{region}-{cov}X-{div}d-{het}h-liftover.vcf.gz"),
        liftover_index = os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "regions", ref_str + "-{region}-{cov}X-{div}d-{het}h-liftover.vcf.gz.tbi"),
        ## Not actually sure how LiftoverVcf creates the index...
        reject_file = os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "regions", ref_str + "-{region}-{cov}X-{div}d-{het}h-liftover-rejected.vcf.gz")
    log:
        os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "logs", ref_str + "-{region}-{cov}X-{div}d-{het}h-liftover-vcf.log")
    resources:
        mem = "24g"        
    shell:
        """
        picard LiftoverVcf I={input.vcf} O={output.liftover_vcf} CHAIN={input.chain} REJECT={output.reject_file} R={input.sim_ref} &> {log}
        """
# Liftover the coordinates of the called variants from the mouse reference to the simulated divergent genome (same as golden vcf here)
# https://gatk.broadinstitute.org/hc/en-us/articles/360036831351-LiftoverVcf-Picard-
## NOT USED, only comparing called variants to original golden variants, both of which are in ref coords

#################

rule liftover_bams:
    input:
        bam = os.path.join(outdir, "mapped-reads", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h.bam"),
        chain = os.path.join(outdir, "simulated-genomes", "{cov}X", "{div}", ref_str + "-consensus.chain")
    output:
        bam = os.path.join(outdir, "mapped-reads", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h-crossmap.bam"),
        bam_sorted = os.path.join(outdir, "mapped-reads", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h-crossmap.sorted.bam"),
        bam_index = os.path.join(outdir, "mapped-reads", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h-crossmap.sorted.bam.bai")
    params:
        out_prefix = os.path.join(outdir, "mapped-reads", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h-crossmap"),
    log:
        os.path.join(outdir, "mapped-reads", "{cov}X", "{div}", "logs", ref_str + "-{cov}X-{div}d-{het}h-crossmap.log")
    resources:
        mem = "24g",
    shell:
        """
        CrossMap.py bam -a {input.chain} {input.bam} {params.out_prefix} &> {log}
        """
# Liftover the coordinates of the mapped reads from the mouse reference to the simulated divergent genome (same as golden bam here)
# https://crossmap.sourceforge.net/#convert-bam-cram-sam-format-files 

#################

# rule index_liftover_bams:
#     input:
#         os.path.join(outdir, "mapped-reads", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h-crossmap.bam")
#     output:
#         os.path.join(outdir, "mapped-reads", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h-crossmap.bam.bai")
#     log:
#         os.path.join(outdir, "mapped-reads", "{cov}X", "{div}", "logs", ref_str + "-{cov}X-{div}d-{het}h-crossmap-index.log")
#     resources:
#         time = "4:00:00",
#         mem = "4g"
#     shell:
#         """
#         samtools index {input} 2> {log}
#         """
# Index the lifted over mapped reads
## NOT USED, crossmap sortes, indexes bam automatically

######################################################################################################

rule compare_vcfs:
    input:
        golden = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "divergent", "regions", ref_str + "-{region}_golden.vcf.gz"),
        golden_index = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "divergent", "regions", ref_str + "-{region}_golden.vcf.gz.tbi"), 
        query = os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "regions", ref_str + "-{region}-{cov}X-{div}d-{het}h.vcf.gz"),
        query_index = os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", "regions", ref_str + "-{region}-{cov}X-{div}d-{het}h.vcf.gz.tbi"),
    params:
        outdir = os.path.join(outdir, "summary-files", "{cov}X", "{div}", "regions"),
        het = "{het}",
        div = "{div}",
        cov = "{cov}",
        iteration = "1"
    output:
        summary = os.path.join(outdir, "summary-files", "{cov}X", "{div}", "regions", ref_str + "-{region}-{cov}X-{div}d-{het}h-compare-vcf-summary.csv"),
        snps = os.path.join(outdir, "summary-files", "{cov}X", "{div}", "regions", ref_str + "-{region}-{cov}X-{div}d-{het}h-compare-vcf-snps.csv")
    resources:
        mem = "4g",
        time = "24:00:00"
    retries: 0
    shell:
        """
        python /n/home07/gthomas/projects/Mapping-simulations/scripts/compare_vcfs_2.py {params.cov} {params.div} {params.het} {params.iteration} {input.golden} {input.query} {params.outdir} {output.summary} {output.snps}
        """
# Run the compare_vcfs script to get number of variants compared to golden file
# Use to combine files:
# for f in *.csv; do (cat "${f}"; echo); done | grep -v "#" | sort -r | uniq

#################

rule compare_bams:
    input:
        golden = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", ref_str + "_golden.bam"),
        golden_index = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", ref_str + "_golden.bam.bai"),
        query = os.path.join(outdir, "mapped-reads", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h-crossmap.sorted.bam"),
        query_index = os.path.join(outdir, "mapped-reads", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h-crossmap.sorted.bam.bai"),        
        sim_index = os.path.join(outdir, "simulated-genomes", "{cov}X", "{div}", ref_str + "-consensus.fa.fai")
    params:
        het = "{het}",
        div = "{div}",
        cov = "{cov}",
        reg = ",".join(regions),
        iteration = "1"
    output:
        summary = os.path.join(outdir, "summary-files", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h-compare-bam-summary.csv"),
    resources:
        mem = "24g",
        time = "8:00:00"
    shell:
        """
        python /n/home07/gthomas/projects/Mapping-simulations/scripts/compare_bams.py {params.reg} {params.cov} {params.div} {params.het} {params.iteration} {input.sim_index} {input.golden} {input.query} {output.summary}
        """
# Run the compare_vcfs script to get number of variants compared to golden file
# Use to combine files:
# for f in *.csv; do (cat "${f}"; echo); done | grep -v "#" | sort -r | uniq

# find . -type f -name *-vcf-summary.csv -not -name mm39-20X-vcf-summary.csv -exec awk 1 {} \; | grep -v "#" | sort -r | uniq > ../mm39-20X-vcf-summary.csv
# find . -type f -name *-vcf-snps.csv -not -name mm39-20X-vcf-snps.csv -exec awk 1 {} \; > ../mm39-20X-vcf-snps.csv
# find . -type f -name *-bam-summary.csv -exec awk 1 {} \; | grep -v "#" | sort -r | uniq > ../mm39-20X-bam-summary.csv

#################

def get_runtime_qualimap(wildcards, attempt):
    times = ["8:00:00", "12:00:00", "16:00:00", "20:00:00", "32:00:00"];
    return times[attempt];
# Sets the run time based on the number of times the job has been restarted
rule qualimap_bams:
    input:
        golden_bam = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", ref_str + "_golden.bam"),
        mapped_bam = os.path.join(outdir, "mapped-reads", "{cov}X", "{div}", ref_str + "-{cov}X-{div}.bam")
    output:
        golden_qualimap = os.path.join(outdir, "qualimap", "{cov}X", "{div}", "golden", "qualimapReport.html"),
        mapped_qualimap = os.path.join(outdir, "qualimap", "{cov}X", "{div}", "mapped", "qualimapReport.html")
    params:
        golden_outdir = os.path.join(outdir, "qualimap", "{cov}X", "{div}", "golden"),
        mapped_outdir = os.path.join(outdir, "qualimap", "{cov}X", "{div}", "mapped")
    log:
        golden_log = os.path.join(outdir, "qualimap", "{cov}X", "{div}", "logs", ref_str + "_golden.log"),
        mapped_log = os.path.join(outdir, "qualimap", "{cov}X", "{div}", "logs", ref_str + "_mapped.log")
    resources:
        cpus = 32,
        time = get_runtime_qualimap
    shell:
        """
        qualimap bamqc -bam {input.golden_bam} -nt {resources.cpus} -outdir {params.golden_outdir} -outformat html --java-mem-size=12G &> {log.golden_log}
        qualimap bamqc -bam {input.mapped_bam} -nt {resources.cpus} -outdir {params.mapped_outdir} -outformat html --java-mem-size=12G &> {log.mapped_log}
        """        
# Run QualiMap on the BAM files
## NOT USED

#################
#############################################################################


# time -p samtools sort -n -o mm39_golden.name.sorted.bam mm39_golden.bam
# name sorting for wub?