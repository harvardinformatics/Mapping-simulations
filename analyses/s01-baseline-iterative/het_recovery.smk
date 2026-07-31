#############################################################################
# s01 heterozygous-site recovery by iteration (add-on to workflow.smk)
#
# The main s01 scoring measures DIVERGENCE (homozygous) recovery via the
# consensus-vs-reference diff. Heterozygous sites are deliberately NOT baked into the
# (haploid) consensus, so that scoring reports het_tp = 0 -- it cannot see het recovery.
#
# Het sites are the classic allelic-dropout signature of reference bias: a true 0/1
# site whose alt-bearing reads map poorly looks like homozygous reference. To track
# how iteration rescues them, we score the RAW per-iteration caller VCF (calls.norm,
# which carries 0/1 genotypes) against the truth with compare_vcfs.py, and read the
# het_* columns (het_tp / n_truth_het = het recall). Positions are identity (SNP-only,
# and het sites are never baked into the consensus, so their consensus base == the
# reference base) -> the caller VCF at every iteration is directly comparable to the
# reference-coordinate truth, no liftover.
#
# Note: the DIVERGENCE columns from THIS comparison are only meaningful at iteration 1;
# at iter>1 the homozygous divergence is already in the consensus and no longer re-called,
# so use the main summary's div_tp for divergence and this file's het_tp for het.
#
# This only READS existing calls.norm VCFs; it does not touch workflow.smk's outputs and
# needs no simulation (cheap: one compare_vcfs per condition).
#
# Launch (activated mapping-sims env, in this directory):
#   snakemake -p -s het_recovery.smk --configfile config/config.yaml \
#       --profile ../../snakemake/profiles/slurm_profile/ --dryrun
#############################################################################

import os

## --- config (same file as the main workflow) -----------------------------
ref_str    = config["ref_str"]
REGION     = config["region"]
LIB        = config["lib_dir"]
outdir     = config["data_dir"]
resultsdir = config["results_dir"]

covs   = config["covs"]
divs   = config["divs"]
HET    = config["het"]
SEEDS  = config["seeds"]
NITER  = config["n_iterations"]
CALLERS = config["callers"]
TAG    = config["run_tag"]

ITERS = [str(i) for i in range(1, NITER + 1)]

wildcard_constraints:
    caller = "gatk|bcftools|freebayes",
    div = r"[0-9.]+", het = r"[0-9.]+", cov = r"[0-9]+", seed = r"[0-9]+",
    region = r"[A-Za-z0-9]+", iter = r"[0-9]+"

## --- path templates (mirror workflow.smk) --------------------------------
sim_dir    = os.path.join(outdir, "sim", "{cov}X", "{region}", "{div}", "{het}", "{seed}")
truth_norm = os.path.join(sim_dir, ref_str + "-truth.norm.vcf.gz")     # normalized truth (has 0/1 GTs)

iter_dir   = os.path.join(outdir, "iter", "{caller}", "{cov}X", "{region}", "{div}", "{het}", "{seed}", "iter{iter}")
calls_norm = os.path.join(iter_dir, ref_str + "_calls.norm.vcf.gz")    # existing raw caller VCF
het_prefix = os.path.join(iter_dir, ref_str + "-hetscore")
het_csv    = het_prefix + "-vcf-summary.csv"

combined   = os.path.join(resultsdir, "s01-" + TAG + "-het-recovery-by-iter.csv")

localrules: all, combine_het

rule all:
    input:
        combined

## --- 1. score the raw caller VCF vs truth (read het_* columns) --------------
rule het_score:
    input:
        called = calls_norm, ctbi = calls_norm + ".tbi",
        truth  = truth_norm, ttbi = truth_norm + ".tbi"
    output:
        csv = het_csv
    params:
        cmp_script = os.path.join(LIB, "compare_vcfs.py")
    log:
        os.path.join(iter_dir, "logs", ref_str + "-hetscore.log")
    resources:
        cpus = 1, mem = "8g", time = "1:00:00", partition = "shared"
    shell:
        """
        mkdir -p $(dirname {log})
        PREF="{output.csv}"; PREF="${{PREF%-vcf-summary.csv}}"   # compare_vcfs appends -vcf-summary.csv
        python {params.cmp_script} \
            --truth-vcf {input.truth} --called-vcf {input.called} \
            --region {wildcards.region} --coverage {wildcards.cov} \
            --divergence {wildcards.div} --heterozygosity {wildcards.het} \
            --seed {wildcards.seed} --caller {wildcards.caller} --iteration {wildcards.iter} \
            --out-prefix $PREF &> {log}
        """

## --- 2. concatenate into one table -----------------------------------------
rule combine_het:
    input:
        expand(het_csv, caller=CALLERS, cov=covs, region=[REGION], div=divs,
               het=[HET], seed=SEEDS, iter=ITERS)
    output:
        combined
    params:
        script = os.path.join(LIB, "combine_summaries.py")
    shell:
        "python {params.script} --inputs {input} --out {output} --kind vcf"
