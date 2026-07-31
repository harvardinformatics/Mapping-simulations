#############################################################################
# Variant-caller benchmark for the v2 mason pipeline.
#
# Holds MAPPING FIXED (the existing bwa BAMs) and varies only the CALLER, to
# isolate how caller choice affects SNP/indel recovery as divergence increases:
#   gatk      - GATK HaplotypeCaller (reference-anchored De Bruijn assembly)
#   bcftools  - bcftools mpileup/call (pileup, no assembly)
#   freebayes - freebayes (haplotype-from-alignment, no k-mer reassembly)
#
# Each caller's NATIVE calls (no post-hoc MQ/DP hard filter) are normalized
# identically -- `bcftools norm -f ref -m- --atomize` (split multiallelics,
# decompose MNP/complex into atomic SNPs+indels) -- so SNP vs indel is classified
# the same way for every caller. Precision then captures each caller's FP tradeoff.
#
# GATK is NOT re-run: we reuse the raw GenotypeGVCFs VCF from the main sweep.
#
# Launch (from an activated mapping-sims env, in snakemake/):
#   snakemake -p -s caller_benchmark.smk \
#       --configfile ../simulation-configs/mm39-mason.yaml \
#       --profile profiles/slurm_profile/ --dryrun
#############################################################################

import os

ref_str   = config["ref_str"]
SCRIPTS   = config["scripts_dir"]
regions   = config["regions"]
covs      = config["covs"]
divs      = config["divs"]
hets      = config["hets"]
seeds     = config["seeds"]
outdir    = config["sim_outdir"]
summarydir = config["summary_dir"]

CALLERS = ["gatk", "bcftools", "freebayes"]

wildcard_constraints:
    caller = "gatk|bcftools|freebayes",
    div = r"[0-9.]+", het = r"[0-9.]+", cov = r"[0-9]+", seed = r"[0-9]+",
    region = r"[A-Za-z0-9]+"

# --- reused inputs from the main sweep --------------------------------------
region_fa  = os.path.join(outdir, "ref-regions", ref_str + "-{region}.fa")
truth_norm = os.path.join(outdir, "truth", "{region}", "{div}", "{het}", "{seed}",
                          ref_str + "-{region}-{div}d-{het}h-truth.norm.vcf.gz")
mapped_bam = os.path.join(outdir, "mapped-reads", "{cov}X", "{region}", "{div}", "{het}",
                          "{seed}", ref_str + "_mapped.bam")
gatk_raw   = os.path.join(outdir, "called-variants", "{cov}X", "{region}", "{div}", "{het}",
                          "{seed}", ref_str + "_raw.vcf.gz")

# --- benchmark outputs ------------------------------------------------------
CB        = os.path.join(outdir, "caller-benchmark", "{caller}", "{cov}X", "{region}",
                         "{div}", "{het}", "{seed}")
calls_norm = os.path.join(CB, "calls.norm.vcf.gz")            # {caller} wildcard
summ       = os.path.join(CB, ref_str + "-vcf-summary.csv")   # {caller} wildcard
combined   = os.path.join(summarydir, ref_str + "-mason-chr19-caller-benchmark.csv")

def norm_out(caller):
    return calls_norm.replace("{caller}", caller)

#############################################################################

localrules: all, combine_caller_benchmark

rule all:
    input:
        combined

# --- per-caller calling (each writes the uniform calls.norm.vcf.gz) ---------

rule call_gatk:
    # reuse the main sweep's raw GenotypeGVCFs output; just normalize it
    input:
        raw = gatk_raw,
        fa  = region_fa
    output:
        vcf = norm_out("gatk"),
        tbi = norm_out("gatk") + ".tbi"
    resources:
        cpus = 1, mem = "8g", time = "1:00:00", partition = "shared"
    shell:
        """
        bcftools norm -f {input.fa} -m- --atomize -Oz -o {output.vcf} {input.raw}
        tabix -p vcf {output.vcf}
        """

rule call_bcftools:
    input:
        bam = mapped_bam,
        fa  = region_fa
    output:
        vcf = norm_out("bcftools"),
        tbi = norm_out("bcftools") + ".tbi"
    log:
        norm_out("bcftools").replace("calls.norm.vcf.gz", "bcftools.log")
    resources:
        cpus = 4, mem = "8g", time = "8:00:00", partition = "shared"
    shell:
        """
        bcftools mpileup -f {input.fa} --min-MQ 20 -a FORMAT/AD,FORMAT/DP \
            --threads {resources.cpus} {input.bam} 2> {log} \
          | bcftools call -mv --threads {resources.cpus} 2>> {log} \
          | bcftools norm -f {input.fa} -m- --atomize -Oz -o {output.vcf} 2>> {log}
        tabix -p vcf {output.vcf}
        """

rule call_freebayes:
    input:
        bam = mapped_bam,
        fa  = region_fa
    output:
        vcf = norm_out("freebayes"),
        tbi = norm_out("freebayes") + ".tbi"
    log:
        norm_out("freebayes").replace("calls.norm.vcf.gz", "freebayes.log")
    resources:
        cpus = 1, mem = "16g", time = "12:00:00", partition = "shared"
    shell:
        """
        freebayes -f {input.fa} --min-mapping-quality 20 {input.bam} 2> {log} \
          | bcftools norm -f {input.fa} -m- --atomize -Oz -o {output.vcf} 2>> {log}
        tabix -p vcf {output.vcf}
        """

# --- score each caller against the shared normalized truth ------------------

rule score_caller:
    input:
        truth  = truth_norm,
        called = calls_norm,
        script = os.path.join(SCRIPTS, "compare_vcfs.py")
    output:
        csv   = summ,
        dict_ = summ.replace("-vcf-summary.csv", "-vcf-summary.dict.txt")
    params:
        prefix = lambda w, output: output.csv[:-len("-vcf-summary.csv")],
        region = "{region}", cov = "{cov}", div = "{div}", het = "{het}",
        seed = "{seed}", caller = "{caller}"
    resources:
        cpus = 1, mem = "8g", time = "2:00:00", partition = "shared"
    shell:
        """
        python {input.script} --truth-vcf {input.truth} --called-vcf {input.called} \
            --region {params.region} --coverage {params.cov} \
            --divergence {params.div} --heterozygosity {params.het} --seed {params.seed} \
            --caller {params.caller} --out-prefix {params.prefix}
        """

rule combine_caller_benchmark:
    input:
        expand(summ, caller=CALLERS, cov=covs, region=regions, div=divs, het=hets, seed=seeds)
    output:
        combined
    params:
        script = os.path.join(SCRIPTS, "combine_summaries.py")
    shell:
        "python {params.script} --inputs {input} --out {output} --kind vcf"

#############################################################################
