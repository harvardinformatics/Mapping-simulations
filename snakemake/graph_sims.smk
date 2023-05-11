#############################################################################
# Pipeline for read mapping simulations with varying divergence
#############################################################################

import os
import re

#############################################################################
# Example cmd for mouse genome

# snakemake -p -s graph_sims.smk --configfile ../simulation-configs/mm39-graph.yaml --profile profiles/slurm_profile/ --dryrun

# To generate rulegraph image:
# snakemake -p -s graph_sims.smk --configfile ../simulation-configs/mm39.yaml --profile profiles/slurm_profile/ --dryrun --rulegraph | dot -Tpng > dag.png

#############################################################################
# Reference file and path info

ref_str = config["ref_str"];
# The subfolder and abbreviation for the current reference genome

REF = config["ref_file"];
# The full reference genome file

REF_INDEX = config["ref_index"];
# The index for the reference genome, from samtools faidx

outdir = config["sim_outdir"];
# The directory for all output for the current reference genome

##########
# Other dependencies

indir = config["sim_indir"];
# The path to the files from the simulated genomes (output from simulate_genomes.smk)

NEAT_PATH = config["neat_path"];
# The path to the NEAT read simulator

VG_PATH = config["vg_path"];
# The path to the vg genome graph builder/mapper

###########

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
        expand(os.path.join(outdir, "summary-files", "{cov}X", ref_str + "-{cov}X-{het}h-bam-summary.csv"), cov=covs, het=hets)
## The final expected outputs should be listed in this rule. Only necessary to list final output from final rule, but I found it useful to list them 
## for all rules for debugging (can comment out outputs for rules you don't want to run), though there's also probably a better way to do this

#############################################################################
# Pipeline functions

def getRuntime(wildcards, attempt, max_time, base_time, multiplier):
    runtime = base_time + ( float(wildcards.div) * attempt * multiplier );
    if runtime > max_time:
        runtime = max_time;
    #print(base_time, wildcards.div, attempt, multiplier, runtime);
    return str(round(runtime)) + ":00:00";
# Sets the run time based on the number of times the job has been restarted
# Necessary for some rules that may time out based on the divergence of the simulation

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

######################################################################################################

rule simulate_heterozygous_reads:
    input:
        os.path.join(indir, "simulated-genomes", "{cov}X", "{div}", "regions", ref_str + "-{region}.fa")
    output:
        read1 = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", "regions", ref_str + "-{region}_read1.fq.gz"),
        read2 = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", "regions", ref_str + "-{region}_read2.fq.gz"),
        bam = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", "regions", ref_str + "-{region}_golden.bam"),
        vcf = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", "regions", ref_str + "-{region}_golden.vcf.gz"),
    params:
        prefix = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", "regions", ref_str + "-{region}"),
        het = "{het}",
        cov = "{cov}",
        neat_path = NEAT_PATH
    log:
        os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", "logs", ref_str + "-{region}.log")
    resources:
        mem = "12g",
        time = lambda wc, attempt: getRuntime(wc, attempt, max_time=168, base_time=6, multiplier=100)
    retries: 3
    shell:
        """
        python {params.neat_path} -r {input} -o {params.prefix} --bam --vcf -R 150 --pe 300 30 -c {params.cov} -M {params.het} &> {log}
        """
# Simulate reads per chromosome with varying levels of divergence (and possibly varying coverage) with NEAT
# This generates read pairs, a golden VCF and a golden BAM for each chromosome
# https://github.com/ncsa/NEAT

#################

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
        time = lambda wc, attempt: getRuntime(wc, attempt, max_time=168, base_time=6, multiplier=100)
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
#         os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", "regions", ref_str + "-{region}_golden.vcf.gz")
#     output:
#         os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", "regions", ref_str + "-{region}_golden.vcf.gz.tbi")
#     log:
#         os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", "logs", ref_str + "-{region}-tabix.log")
#     resources:
#         mem = "2g",
#         time = "2:00:00"
#     shell:
#         """
#         tabix {input} &> {log}
#         """
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
        mem="12g",
        time = lambda wc, attempt: getRuntime(wc, attempt, max_time=168, base_time=160, multiplier=200)
    retries: 3
    shell:
        """
        gatk HaplotypeCaller -R {input.ref} -I {input.bam} -stand-call-conf 30 --native-pair-hmm-threads {resources.cpus} -L {params.region} -ERC GVCF -O {output} &> {log}
        """
# Call variants with GATK HaplotypeCaller by input region

#################

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
        time = lambda wc, attempt: getRuntime(wc, attempt, max_time=168, base_time=6, multiplier=100)
    retries: 3
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

#################

rule build_graph:
    input:
        ref = REF,
        vcf = os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h.vcf.gz"),
        vcf_index = os.path.join(outdir, "called-variants", "gatk", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h.vcf.gz.tbi")
    output:
        os.path.join(outdir, "graphs", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h.giraffe.gbz")
    params:
        prefix = ref_str + "-{cov}X-{div}d-{het}h",
        vg_path = VG_PATH
    log:
        os.path.join(outdir, "graphs", "{cov}X", "{div}", "logs", ref_str + "-{cov}X-{div}d-{het}h-autoindex.log")
    resources:
        cpus = 12,
        mem = "24g"
    shell:
        """
        {params.vg_path} autoindex --workflow giraffe -t {resources.cpus} -r {input.ref} -v {input.vcf} -p {params.prefix}
        """

#################

rule map_graph:
    input:
        graph = os.path.join(outdir, "graphs", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h.giraffe.gbz"),
        read1 = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", ref_str + "_read1.fq.gz"),
        read2 = os.path.join(outdir, "simulated-reads", "{cov}X", "{div}", "heterozygous", "{het}", ref_str + "_read2.fq.gz"),
    output:
        os.path.join(outdir, "graphs", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h.bam")
    params:
        prefix = ref_str + "-{cov}X-{div}d-{het}h",
        vg_path = VG_PATH
    resources:
        cpus = 12,
        mem = "24g"
    shell:
        """
        {params.vg_path} giraffe -t {resources.cpus} -Z {input.graph} -f {input.read1} -f {input.read2} -o bam -p > {output}
        """

#################

rule liftover_bams:
    input:
        bam = os.path.join(outdir, "graphs", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h.bam"),
        chain = os.path.join(indir, "simulated-genomes", "{cov}X", "{div}", ref_str + "-consensus.chain")
    output:
        bam = os.path.join(outdir, "graphs", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h-crossmap.bam"),
        bam_sorted = os.path.join(outdir, "graphs", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h-crossmap.sorted.bam"),
        bam_index = os.path.join(outdir, "graphs", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h-crossmap.sorted.bam.bai")
    params:
        out_prefix = os.path.join(outdir, "graphs", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h-crossmap"),
    log:
        os.path.join(outdir, "graphs", "{cov}X", "{div}", "logs", ref_str + "-{cov}X-{div}d-{het}h-crossmap.log")
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
        query = os.path.join(outdir, "graphs", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h-crossmap.sorted.bam"),
        query_index = os.path.join(outdir, "graphs", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h-crossmap.sorted.bam.bai"),        
        sim_index = os.path.join(indir, "simulated-genomes", "{cov}X", "{div}", ref_str + "-consensus.fa.fai")
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

######################################################################################################

rule combine_summaries:
    input:
        bam_summary = expand(os.path.join(outdir, "summary-files", "{cov}X", "{div}", ref_str + "-{cov}X-{div}d-{het}h-compare-bam-summary.csv"), cov=covs, div=divs, het=hets)
    output:
        bam_out = os.path.join(outdir, "summary-files", "{cov}X", ref_str + "-{cov}X-{het}h-bam-summary.csv")
    params:
        indir = os.path.join(outdir, "summary-files", "{cov}X")
    resources:
        mem = "4g",
        time = "1:00:00"
    shell:
        """
        find {params.indir} -type f -name *-compare-bam-summary.csv -exec awk '/^[^#]/' {{}} \; | sort -r | uniq > {output.bam_out}
        """
# Combines the output files from compare_vcfs and compare_bams
# Originally ran separately as:
# find . -type f -name *-vcf-summary.csv -not -name mm39-20X-vcf-summary.csv -exec awk 1 {} \; | grep -v "#" | sort -r | uniq > ../mm39-20X-vcf-summary.csv
# find . -type f -name *-vcf-snps.csv -not -name mm39-20X-vcf-snps.csv -exec awk 1 {} \; > ../mm39-20X-vcf-snps.csv
# find . -type f -name *-bam-summary.csv -exec awk 1 {} \; | grep -v "#" | sort -r | uniq > ../mm39-20X-bam-summary.csv


#############################################################################

