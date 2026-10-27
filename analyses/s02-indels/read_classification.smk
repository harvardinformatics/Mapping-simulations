#############################################################################
# s02 read-classification-by-iteration (indels) -- add-on to workflow.smk.
#
# Mapping-side mechanism, with indels. Classifies each iteration's mapped BAM against the
# mason -oa truth alignment (reference coords) with compare_bams.py. Iteration 1 maps to
# the reference (identity coords -> plain compare_bams). Iteration k>1 maps to the
# length-shifted consensus_{k-1}, so we point-lift each read's start/end back to reference
# coords via `--chain cons_{k-1}->ref` (the chainSwap'd chain from the previous build).
# The point-lift keeps indel-spanning reads (CrossMap bam would drop them); the end-check
# flags reads that span the wrong indel footprint (span_off).
#
# Consumes outputs of workflow.smk (truth.oa.bam, per-iteration BAMs, cons2ref chains);
# run AFTER the recovery workflow. Launch:  sbatch run_read_classification.sh
#############################################################################

import os

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

## --- paths (mirror workflow.smk) -----------------------------------------
sim_dir    = os.path.join(outdir, "sim", "{cov}X", "{region}", "{div}", "{het}", "{seed}")
truth_oa   = os.path.join(sim_dir, ref_str + "-truth.oa.bam")
truth_vcf  = os.path.join(sim_dir, ref_str + "-truth.vcf.gz")     # mason format, as fed to -iv

iter_dir   = os.path.join(outdir, "iter", "{caller}", "{cov}X", "{region}", "{div}", "{het}", "{seed}", "iter{iter}")
mapped_bam = os.path.join(iter_dir, ref_str + "_mapped.bam")
rc_prefix  = os.path.join(iter_dir, ref_str + "-readclass")
rc_csv     = rc_prefix + "-bam-summary.csv"

combined   = os.path.join(resultsdir, "s02-" + TAG + "-readclass-by-iter.csv")

def _base_ref(region):
    return os.path.join(outdir, "ref-regions", ref_str + "-" + region + ".fa")

# chain to lift THIS iteration's BAM (consensus_{k-1} coords) back to reference: iter1 has
# none (base-ref coords) -> dummy base_ref, ignored in-shell; else cons_{k-1}->ref chain.
def prev_chain(w):
    if int(w.iter) == 1:
        return _base_ref(w.region)
    prev = os.path.join(outdir, "iter", w.caller, w.cov + "X", w.region, w.div, w.het,
                        w.seed, "iter" + str(int(w.iter) - 1))
    return os.path.join(prev, "cons2ref.chain")

localrules: all, combine_readclass

rule all:
    input:
        combined

rule classify_reads:
    input:
        truth  = truth_oa,
        mapped = mapped_bam,
        chain  = prev_chain,
        mason_vcf = truth_vcf,     # rebuilds the haplotypes -> each read's true span from its oP tag
        fa     = lambda w: _base_ref(w.region)
    output:
        csv = rc_csv
    params:
        cmp = os.path.join(LIB, "compare_bams.py")
    log:
        os.path.join(iter_dir, "logs", ref_str + "-readclass.log")
    resources:
        cpus = 1, mem = "12g", time = "2:00:00", partition = "shared"   # peak ~6.5 GB (truth spans + haplotypes)
    shell:
        # iter1: identity coords, no chain. iter>1: point-lift start+end via cons_{k-1}->ref.
        r"""
        mkdir -p $(dirname {log})
        PREF="{output.csv}"; PREF="${{PREF%-bam-summary.csv}}"
        CHAINARG=""
        [ {wildcards.iter} -ne 1 ] && CHAINARG="--chain {input.chain}"
        python {params.cmp} \
            --truth-bam {input.truth} --mapped-bam {input.mapped} $CHAINARG \
            --mason-vcf {input.mason_vcf} --reference {input.fa} \
            --region {wildcards.region} --coverage {wildcards.cov} \
            --divergence {wildcards.div} --heterozygosity {wildcards.het} \
            --seed {wildcards.seed} --caller {wildcards.caller} --iteration {wildcards.iter} \
            --out-prefix $PREF &> {log}
        """

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
