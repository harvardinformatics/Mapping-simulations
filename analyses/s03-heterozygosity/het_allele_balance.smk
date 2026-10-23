#############################################################################
# s03 — reference-allele bias at heterozygous sites.
#
# Measures the reference-allele fraction (RAF) at true het SNP sites straight from the
# pileup of each iteration's BAM (allele_balance.py) -- the canonical reference-bias metric,
# independent of the variant caller. Consumes EXISTING BAMs (s01 or s02); no new simulation.
#
# Parameterized by config so the SAME workflow runs on either source:
#   - s01 (SNP-only): identity coordinates, no liftover (has_liftover: false).
#   - s02 (with indels): iter>1 BAMs are in consensus coords, so het positions are lifted
#     ref->consensus with that iteration's ref2cons chain (has_liftover: true).
#
# Launch:  sbatch run_allele_balance.sh   (runs both configs)
#############################################################################

import os

ref_str    = config["ref_str"]
REGION     = config["region"]
LIB        = config["lib_dir"]
SOURCE     = config["source_data_dir"]     # data/s01-baseline-iterative or data/s02-indels
resultsdir = config["results_dir"]
HAS_LIFT   = config["has_liftover"]

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

sim_dir    = os.path.join(SOURCE, "sim", "{cov}X", "{region}", "{div}", "{het}", "{seed}")
truth_vcf  = os.path.join(sim_dir, ref_str + "-truth.vcf.gz")
iter_dir   = os.path.join(SOURCE, "iter", "{caller}", "{cov}X", "{region}", "{div}", "{het}", "{seed}", "iter{iter}")
mapped_bam = os.path.join(iter_dir, ref_str + "_mapped.bam")
ab_csv     = os.path.join(iter_dir, ref_str + "-ab-summary.csv")

combined   = os.path.join(resultsdir, "s03-" + TAG + "-allele-balance-summary.csv")

# ref->consensus chain to lift het positions into the BAM's coords: only for s02 (has_liftover)
# and iter>1, where the BAM maps to consensus_{k-1} (== iter_{k-1}/ref2cons.chain).
def prev_chain(w):
    if HAS_LIFT and int(w.iter) > 1:
        return os.path.join(SOURCE, "iter", w.caller, w.cov + "X", w.region, w.div, w.het,
                            w.seed, "iter" + str(int(w.iter) - 1), "ref2cons.chain")
    return ""

localrules: all, combine_ab

rule all:
    input:
        combined

rule ab_pileup:
    input:
        truth = truth_vcf,
        bam   = mapped_bam,
        bai   = mapped_bam + ".bai"
    output:
        csv = ab_csv
    params:
        script = os.path.join(LIB, "allele_balance.py"),
        chain  = prev_chain,
        prefix = lambda w, output: output.csv[:-len("-ab-summary.csv")]
    log:
        os.path.join(iter_dir, "logs", ref_str + "-ab.log")
    resources:
        cpus = 1, mem = "10g", time = "2:00:00", partition = "shared"
    shell:
        r"""
        mkdir -p $(dirname {log})
        CH=""; [ -n "{params.chain}" ] && [ -f "{params.chain}" ] && CH="--chain {params.chain}"
        python {params.script} --truth-vcf {input.truth} --bam {input.bam} $CH \
            --region {wildcards.region} --coverage {wildcards.cov} \
            --divergence {wildcards.div} --heterozygosity {wildcards.het} --seed {wildcards.seed} \
            --caller {wildcards.caller} --iteration {wildcards.iter} \
            --out-prefix {params.prefix} &> {log}
        """

rule combine_ab:
    input:
        expand(ab_csv, caller=CALLERS, cov=covs, region=[REGION], div=divs,
               het=[HET], seed=SEEDS, iter=ITERS)
    output:
        combined
    params:
        script = os.path.join(LIB, "combine_summaries.py")
    shell:
        "python {params.script} --inputs {input} --out {output} --kind ab"
