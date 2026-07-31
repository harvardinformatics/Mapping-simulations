#############################################################################
# v2 reference-bias mapping simulations (mason2-based), single iteration.
# SNPs + small indels, replicate seeds, chr19 -> chr19 mapping.
#
# Simulate reads from a controlled diverged/diploid version of a chromosome, map
# them back to the reference, call variants, and compare BOTH read placements and
# variant calls against the mason truth (already in reference coordinates -- no
# CrossMap / liftover, unlike v1).
#
# Example launch (from an ACTIVATED mapping-sims env, in the snakemake/ dir):
#   snakemake -p -s simulate_and_map.smk \
#       --configfile ../simulation-configs/mm39-mason.yaml \
#       --profile profiles/slurm_profile/ --dryrun
#############################################################################

import os

#############################################################################
# Config

REF        = config["ref_file"]
REF_INDEX  = config["ref_index"]
ref_str    = config["ref_str"]
SCRIPTS    = config["scripts_dir"]

regions    = config["regions"]
covs       = config["covs"]
divs       = config["divs"]
hets       = config["hets"]
seeds      = config["seeds"]
READLEN    = config["read_length"]
FRAGMEAN   = config["fragment_mean"]
FILTER     = config["filter"]

INDEL_RATIO = config["indel_ratio"]
MIN_INDEL   = config["min_indel_size"]
MAX_INDEL   = config["max_indel_size"]

outdir     = config["sim_outdir"]
summarydir = config["summary_dir"]

#############################################################################
# Sanity checks

if not os.path.isfile(REF):
    raise FileNotFoundError("\n * ERROR: cannot find reference genome: " + REF + "\n")
if not os.path.isfile(REF_INDEX):
    raise FileNotFoundError("\n * ERROR: cannot find reference .fai index: " + REF_INDEX + "\n")

# div/het are decimals, seed is digits -- constrain so wildcards parse unambiguously
wildcard_constraints:
    div = r"[0-9.]+",
    het = r"[0-9.]+",
    cov = r"[0-9]+",
    seed = r"[0-9]+",
    region = r"[A-Za-z0-9]+"

#############################################################################
# Path templates

region_fa   = os.path.join(outdir, "ref-regions", ref_str + "-{region}.fa")

truth_base  = os.path.join(outdir, "truth", "{region}", "{div}", "{het}", "{seed}",
                           ref_str + "-{region}-{div}d-{het}h")
truth_vcf   = truth_base + "-truth.vcf.gz"            # mason format (FORMAT=".")
truth_norm  = truth_base + "-truth.norm.vcf.gz"       # normalized, for comparison

reads_dir   = os.path.join(outdir, "simulated-reads", "{cov}X", "{region}", "{div}", "{het}", "{seed}")
read1       = os.path.join(reads_dir, ref_str + "_read1.fq.gz")
read2       = os.path.join(reads_dir, ref_str + "_read2.fq.gz")
truth_bam   = os.path.join(reads_dir, ref_str + "_truth.bam")

mapped_dir  = os.path.join(outdir, "mapped-reads", "{cov}X", "{region}", "{div}", "{het}", "{seed}")
mapped_bam  = os.path.join(mapped_dir, ref_str + "_mapped.bam")

call_dir    = os.path.join(outdir, "called-variants", "{cov}X", "{region}", "{div}", "{het}", "{seed}")
called_norm = os.path.join(call_dir, ref_str + "_called-norm.vcf.gz")

summ_dir    = os.path.join(outdir, "summary-files", "{cov}X", "{region}", "{div}", "{het}", "{seed}")
bam_summary = os.path.join(summ_dir, ref_str + "-bam-summary.csv")
vcf_summary = os.path.join(summ_dir, ref_str + "-vcf-summary.csv")

combined_bam = os.path.join(summarydir, ref_str + "-mason-chr19-bam-summary.csv")
combined_vcf = os.path.join(summarydir, ref_str + "-mason-chr19-vcf-summary.csv")

#############################################################################
# Final targets

localrules: all, combine_bam_summaries, combine_vcf_summaries

rule all:
    input:
        combined_bam,
        combined_vcf

#############################################################################
# 1. Extract the region from the reference and build all needed indexes.

rule extract_region_ref:
    input:
        ref = REF
    output:
        fa    = region_fa,
        fai   = region_fa + ".fai",
        bwa   = multiext(region_fa, ".amb", ".ann", ".bwt", ".pac", ".sa"),
        dict_ = region_fa.replace(".fa", ".dict")
    params:
        region = "{region}"
    resources:
        cpus = 1, mem = "8g", time = "2:00:00", partition = "shared"
    shell:
        """
        samtools faidx {input.ref} {params.region} > {output.fa}
        samtools faidx {output.fa}
        bwa index {output.fa}
        rm -f {output.dict_}
        picard CreateSequenceDictionary R={output.fa} O={output.dict_}
        """

#############################################################################
# 2. Controlled truth VCF (SNP + indel; homozygous divergence + heterozygous).
#    Indel rates scale with the SNP rates via indel_ratio.

rule make_truth_vcf:
    input:
        fa = region_fa
    output:
        vcf = truth_vcf,
        tbi = truth_vcf + ".tbi"
    params:
        script = os.path.join(SCRIPTS, "make_truth_vcf.py"),
        div = "{div}", het = "{het}", seed = "{seed}",
        idiv = lambda w: "%.6f" % (float(w.div) * INDEL_RATIO),
        ihet = lambda w: "%.6f" % (float(w.het) * INDEL_RATIO),
        minind = MIN_INDEL, maxind = MAX_INDEL,
        plain = lambda w, output: output.vcf[:-3]
    resources:
        cpus = 1, mem = "8g", time = "1:00:00", partition = "shared"
    shell:
        """
        python {params.script} --reference {input.fa} \
            --divergence {params.div} --heterozygosity {params.het} \
            --indel-divergence {params.idiv} --indel-heterozygosity {params.ihet} \
            --min-indel-size {params.minind} --max-indel-size {params.maxind} \
            --seed {params.seed} --out {params.plain}
        bgzip -f {params.plain}
        tabix -p vcf {output.vcf}
        """

#############################################################################
# 3. Normalized truth for comparison. mason uses FORMAT="." (unreadable by
#    bcftools), so rewrite FORMAT to "GT", then left-align/normalize + split.

rule normalize_truth:
    input:
        vcf = truth_vcf,
        fa  = region_fa
    output:
        vcf = truth_norm,
        tbi = truth_norm + ".tbi"
    resources:
        cpus = 1, mem = "8g", time = "1:00:00", partition = "shared"
    shell:
        """
        zcat {input.vcf} \
            | awk 'BEGIN{{OFS="\\t"}} /^#/{{print; next}} {{$9="GT"; print}}' \
            | bcftools norm -f {input.fa} -m- -Oz -o {output.vcf}
        tabix -p vcf {output.vcf}
        """

#############################################################################
# 4. Simulate reads + truth alignment BAM with mason_simulator (mason-format VCF).

def sim_time(wildcards, attempt):
    return ["12:00:00", "16:00:00", "24:00:00"][attempt - 1]

rule simulate_reads:
    input:
        fa  = region_fa,
        fai = region_fa + ".fai",
        vcf = truth_vcf,
        tbi = truth_vcf + ".tbi"
    output:
        read1 = read1,
        read2 = read2,
        bam   = truth_bam,
        bai   = truth_bam + ".bai"
    params:
        cov = "{cov}", readlen = READLEN, fragmean = FRAGMEAN, seed = "{seed}",
        unsorted = lambda w, output: output.bam.replace(".bam", ".unsorted.bam")
    log:
        os.path.join(reads_dir, "logs", ref_str + "-mason.log")
    resources:
        cpus = 8, mem = "24g", time = sim_time, partition = "shared"
    shell:
        """
        LEN=$(cut -f2 {input.fai} | head -1)
        NFRAG=$(( {params.cov} * LEN / (2 * {params.readlen}) ))
        echo "# region length $LEN ; {params.cov}X -> $NFRAG fragments" > {log}
        mason_simulator \
            -ir {input.fa} -iv {input.vcf} -n $NFRAG \
            --illumina-read-length {params.readlen} \
            --fragment-mean-size {params.fragmean} \
            --num-threads {resources.cpus} --seed {params.seed} \
            -o {output.read1} -or {output.read2} -oa {params.unsorted} \
            --version-check 0 &>> {log}
        samtools sort -@ {resources.cpus} -o {output.bam} {params.unsorted} 2>> {log}
        samtools index {output.bam} 2>> {log}
        rm -f {params.unsorted}
        """

#############################################################################
# 5. Map reads back to the (chr19) reference with bwa mem (+ read group for GATK).

rule bwa_map:
    input:
        fa    = region_fa,
        bwa   = multiext(region_fa, ".amb", ".ann", ".bwt", ".pac", ".sa"),
        read1 = read1,
        read2 = read2
    output:
        bam = mapped_bam,
        bai = mapped_bam + ".bai"
    params:
        rg = r"@RG\tID:{cov}X-{region}-{div}-{het}-{seed}\tSM:" + ref_str + r"\tLB:lib1\tPL:ILLUMINA"
    log:
        os.path.join(mapped_dir, "logs", ref_str + "-bwa.log")
    resources:
        cpus = 8, mem = "16g", time = "12:00:00", partition = "shared"
    shell:
        """
        bwa mem -t {resources.cpus} -R '{params.rg}' {input.fa} {input.read1} {input.read2} 2> {log} \
            | samtools sort -@ {resources.cpus} -o {output.bam} -
        samtools index {output.bam}
        """

#############################################################################
# 6. Call variants: HaplotypeCaller -> GenotypeGVCFs -> filter -> normalize.
#    Keeps SNPs AND indels (compare_vcfs handles both).

rule gatk_call:
    input:
        fa    = region_fa,
        fai   = region_fa + ".fai",
        dict_ = region_fa.replace(".fa", ".dict"),
        bam   = mapped_bam,
        bai   = mapped_bam + ".bai"
    output:
        vcf = called_norm,
        tbi = called_norm + ".tbi"
    params:
        gvcf = lambda w, output: os.path.join(os.path.dirname(output.vcf), ref_str + "_g.vcf.gz"),
        raw  = lambda w, output: os.path.join(os.path.dirname(output.vcf), ref_str + "_raw.vcf.gz"),
        filt = FILTER
    log:
        os.path.join(call_dir, "logs", ref_str + "-gatk.log")
    resources:
        cpus = 4, mem = "16g", time = "24:00:00", partition = "shared"
    shell:
        """
        gatk HaplotypeCaller -R {input.fa} -I {input.bam} -ERC GVCF -O {params.gvcf} &> {log}
        gatk GenotypeGVCFs -R {input.fa} -V {params.gvcf} -O {params.raw} &>> {log}
        # drop low-confidence sites, then left-align/normalize + split multiallelics
        bcftools view -e '{params.filt}' {params.raw} 2>> {log} \
            | bcftools norm -f {input.fa} -m- -Oz -o {output.vcf} 2>> {log}
        tabix -p vcf {output.vcf}
        """

#############################################################################
# 7. Compare mapped read positions to truth (no liftover).

rule compare_bams:
    input:
        truth  = truth_bam,
        mapped = mapped_bam,
        script = os.path.join(SCRIPTS, "compare_bams.py")
    output:
        csv  = bam_summary,
        dict_ = bam_summary.replace("-bam-summary.csv", "-bam-summary.dict.txt")
    params:
        prefix = lambda w, output: output.csv[:-len("-bam-summary.csv")],
        region = "{region}", cov = "{cov}", div = "{div}", het = "{het}", seed = "{seed}"
    resources:
        cpus = 1, mem = "8g", time = "2:00:00", partition = "shared"
    shell:
        """
        python {input.script} --truth-bam {input.truth} --mapped-bam {input.mapped} \
            --region {params.region} --coverage {params.cov} \
            --divergence {params.div} --heterozygosity {params.het} --seed {params.seed} \
            --out-prefix {params.prefix}
        """

#############################################################################
# 8. Compare called variants (SNPs + indels) to normalized truth.

rule compare_vcfs:
    input:
        truth  = truth_norm,
        called = called_norm,
        script = os.path.join(SCRIPTS, "compare_vcfs.py")
    output:
        csv  = vcf_summary,
        dict_ = vcf_summary.replace("-vcf-summary.csv", "-vcf-summary.dict.txt")
    params:
        prefix = lambda w, output: output.csv[:-len("-vcf-summary.csv")],
        region = "{region}", cov = "{cov}", div = "{div}", het = "{het}", seed = "{seed}"
    resources:
        cpus = 1, mem = "8g", time = "2:00:00", partition = "shared"
    shell:
        """
        python {input.script} --truth-vcf {input.truth} --called-vcf {input.called} \
            --region {params.region} --coverage {params.cov} \
            --divergence {params.div} --heterozygosity {params.het} --seed {params.seed} \
            --out-prefix {params.prefix}
        """

#############################################################################
# 9. Combine per-condition summaries (all cov/region/div/het/seed) into summary-data/.

rule combine_bam_summaries:
    input:
        expand(bam_summary, cov=covs, region=regions, div=divs, het=hets, seed=seeds)
    output:
        combined_bam
    params:
        script = os.path.join(SCRIPTS, "combine_summaries.py")
    shell:
        "python {params.script} --inputs {input} --out {output} --kind bam"

rule combine_vcf_summaries:
    input:
        expand(vcf_summary, cov=covs, region=regions, div=divs, het=hets, seed=seeds)
    output:
        combined_vcf
    params:
        script = os.path.join(SCRIPTS, "combine_summaries.py")
    shell:
        "python {params.script} --inputs {input} --out {output} --kind vcf"

#############################################################################
