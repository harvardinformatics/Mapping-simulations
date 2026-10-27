#############################################################################
# s01 read-classification-by-iteration (add-on to workflow.smk)
#
# Mapping-side mechanism for the s01 result: as the pseudo-it consensus converges on
# the sample, reads should map more exactly (fewer shifted / soft-clipped / unmapped),
# which is *why* the caller then recovers more true SNPs. This workflow measures that.
#
# s01's simulate_reads did NOT keep a mason `-oa` truth alignment, so we regenerate it:
# mason is deterministic, so re-running mason_simulator with the SAME reference, truth
# VCF, read count, and seed reproduces the identical reads plus their true positions.
# We then classify each EXISTING per-iteration mapped BAM (already on disk from the main
# sweep) against that truth with compare_bams.py. Because the SNP-only consensus is the
# same length as the reference (identity coordinates), iterations 2..N are directly
# comparable to the reference-coordinate truth -- no liftover.
#
# This does NOT touch any output of workflow.smk: it only READS the mapped BAMs and
# adds a truth BAM + per-condition classification alongside them.
#
# Launch (from an activated mapping-sims env, in this directory):
#   snakemake -p -s read_classification.smk --configfile config/config.yaml \
#       --profile ../../snakemake/profiles/slurm_profile/ --dryrun
#############################################################################

import os

## --- config (same file as the main workflow) -----------------------------
REF        = config["ref_file"]
ref_str    = config["ref_str"]
REGION     = config["region"]
LIB        = config["lib_dir"]
outdir     = config["data_dir"]
resultsdir = config["results_dir"]

covs   = config["covs"]
divs   = config["divs"]
HET    = config["het"]
SEEDS  = config["seeds"]
READLEN  = config["read_length"]
FRAGMEAN = config["fragment_mean"]
NITER  = config["n_iterations"]
CALLERS = config["callers"]
TAG    = config["run_tag"]

ITERS = [str(i) for i in range(1, NITER + 1)]

wildcard_constraints:
    caller = "gatk|bcftools|freebayes",
    div = r"[0-9.]+", het = r"[0-9.]+", cov = r"[0-9]+", seed = r"[0-9]+",
    region = r"[A-Za-z0-9]+", iter = r"[0-9]+"

## --- path templates (mirror workflow.smk) --------------------------------
region_fa  = os.path.join(outdir, "ref-regions", ref_str + "-{region}.fa")

sim_dir    = os.path.join(outdir, "sim", "{cov}X", "{region}", "{div}", "{het}", "{seed}")
truth_vcf  = os.path.join(sim_dir, ref_str + "-truth.vcf.gz")     # mason format, as fed to -iv
truth_oa   = os.path.join(sim_dir, ref_str + "-truth.oa.bam")     # NEW: regenerated -oa alignment

iter_dir   = os.path.join(outdir, "iter", "{caller}", "{cov}X", "{region}", "{div}", "{het}", "{seed}", "iter{iter}")
mapped_bam = os.path.join(iter_dir, ref_str + "_mapped.bam")      # existing (from workflow.smk)
rc_prefix  = os.path.join(iter_dir, ref_str + "-readclass")
rc_csv     = rc_prefix + "-bam-summary.csv"

combined   = os.path.join(resultsdir, "s01-" + TAG + "-readclass-by-iter.csv")

localrules: all, combine_readclass

rule all:
    input:
        combined

## --- 1. regenerate the mason -oa truth BAM (identical reads, same seed) -----
rule truth_bam:
    input:
        fa  = region_fa, fai = region_fa + ".fai",
        vcf = truth_vcf, tbi = truth_vcf + ".tbi"
    output:
        # ~2.2 GB each; pure intermediate (only classify_reads consumes it) -> temp so
        # snakemake deletes it once all its classifications are done (~40 GB otherwise)
        bam = temp(truth_oa)
    params:
        readlen = READLEN, fragmean = FRAGMEAN
    log:
        os.path.join(os.path.dirname(truth_oa), "logs", ref_str + "-mason-oa.log")
    resources:
        cpus = 8, mem = "24g", time = "8:00:00", partition = "shared"
    shell:
        """
        mkdir -p $(dirname {log})
        LEN=$(cut -f2 {input.fai} | head -1)
        NFRAG=$(( {wildcards.cov} * LEN / (2 * {params.readlen}) ))
        TMP=$(mktemp -d)
        # same invocation as simulate_reads + `-oa` for the truth alignment; the reads
        # themselves go to scratch (we only keep the alignment).
        mason_simulator -ir {input.fa} -iv {input.vcf} -n $NFRAG \
            --illumina-read-length {params.readlen} --fragment-mean-size {params.fragmean} \
            --num-threads {resources.cpus} --seed {wildcards.seed} \
            -o $TMP/r1.fq.gz -or $TMP/r2.fq.gz -oa $TMP/truth.sam \
            --version-check 0 &> {log}
        # compare_bams streams the truth alignment (no random access) -> unsorted BAM is fine
        samtools view -b $TMP/truth.sam -o {output.bam}
        rm -rf $TMP
        """

## --- 2. classify each existing per-iteration mapped BAM vs truth ------------
rule classify_reads:
    input:
        truth  = truth_oa,
        mapped = mapped_bam,
        mason_vcf = truth_vcf,     # rebuilds the haplotypes -> each read's true span from its oP tag
        fa     = region_fa
    output:
        csv = rc_csv
    params:
        cmp_script = os.path.join(LIB, "compare_bams.py")
    log:
        os.path.join(iter_dir, "logs", ref_str + "-readclass.log")
    resources:
        cpus = 1, mem = "12g", time = "2:00:00", partition = "shared"   # peak ~6.5 GB (truth spans + haplotypes)
    shell:
        """
        mkdir -p $(dirname {log})
        PREF="{output.csv}"; PREF="${{PREF%-bam-summary.csv}}"   # compare_bams appends -bam-summary.csv
        python {params.cmp_script} \
            --truth-bam {input.truth} --mapped-bam {input.mapped} \
            --mason-vcf {input.mason_vcf} --reference {input.fa} \
            --region {wildcards.region} --coverage {wildcards.cov} \
            --divergence {wildcards.div} --heterozygosity {wildcards.het} \
            --seed {wildcards.seed} --caller {wildcards.caller} --iteration {wildcards.iter} \
            --out-prefix $PREF &> {log}
        """

## --- 3. concatenate into one table -----------------------------------------
rule combine_readclass:
    input:
        expand(rc_csv, caller=CALLERS, cov=covs, region=[REGION], div=divs,
               het=[HET], seed=SEEDS, iter=ITERS)
    output:
        combined
    params:
        script = os.path.join(LIB, "combine_summaries.py")
    shell:
        "python {params.script} --inputs {input} --out {output} --kind readclass-bam"
