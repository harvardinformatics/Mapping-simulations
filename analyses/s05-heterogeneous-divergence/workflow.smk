#############################################################################
# s05 — heterogeneous divergence (continuous autocorrelated field), all callers.
#
# Identical pipeline to s01 (SNP-only iterative mapping: iteration 1 is the single-shot
# baseline; iteration i>1 maps to the previous consensus; recovered divergence = the
# consensus-vs-reference SNP diff, scored vs truth). The ONLY change is how the truth
# VCF places fixed-divergence sites: instead of uniform placement, divergent (1|1) SNPs
# are drawn with probability proportional to an autocorrelated field (see
# analyses/lib/divergence_field.py), so divergence clusters/thins along the contig while
# the genome-wide mean stays pinned to div_target. Het sites remain uniform.
#
# The varying axis is the field tag {field} = "L{corr_len}-S{sigma}"; div_target is a
# fixed constant (not a wildcard). Each truth VCF also gets a per-site MEASURED local
# divergence table (annotate_local_divergence.py) that downstream analyses bin against.
#
# Because it is SNP-only, every consensus is the same length as the reference (identity
# coordinates, no liftover) -- same reasoning as s01.
#
# Launch (from an activated mapping-sims env):
#   snakemake -p -s workflow.smk --configfile config/config.yaml \
#       --profile ../../snakemake/profiles/slurm_profile/ --dryrun
#############################################################################

import os

## --- config --------------------------------------------------------------
REF        = config["ref_file"]
ref_str    = config["ref_str"]
REGION     = config["region"]
SCRIPTS    = config["scripts_dir"]     # s05-specific
LIB        = config["lib_dir"]         # shared analysis scripts
outdir     = config["data_dir"]
resultsdir = config["results_dir"]

covs   = config["covs"]
DIV    = str(config["div_target"])         # FIXED mean divergence (literal, not a wildcard);
                                           # str() so a --config div_target=0.01 override (parsed
                                           # as a float) still works in path templates
HET    = config["het"]
SEEDS  = config["seeds"]
READLEN  = config["read_length"]
FRAGMEAN = config["fragment_mean"]
IDIV   = config["indel_divergence"]        # 0 (SNP-only)
IHET   = config["indel_heterozygosity"]
NITER  = config["n_iterations"]
CALLERS = config["callers"]
MINMQ  = config["min_mapping_quality"]
GATK_SCATTER = config["gatk_scatter"]
TAG    = config["run_tag"]

FIELDS     = config["field_configs"]       # list of "L{corr}-S{sigma}[-C{cap}]" tags
FIELD_SEED = config["field_seed"]
LOCAL_WINDOWS = config["local_div_windows"]
HET_ON_FIELD = config.get("het_on_field", False)   # place het on the divergence field?

ITERS = [str(i) for i in range(1, NITER + 1)]


def parse_field(tag):
    """'L20000-S1.0-C0.20' -> (corr_len='20000', sigma='1.0', cap='0.20').
    Cap is optional: 'L20000-S1.0' -> (..., cap=None). sigma=0 = homogeneous."""
    parts = tag.split("-")
    corr = parts[0][1:]
    sigma = parts[1][1:]
    cap = parts[2][1:] if len(parts) > 2 and parts[2].startswith("C") else None
    return corr, sigma, cap


def field_opts(w):
    """make_truth_vcf field options for a config tag (corr/sigma/seed, optional cap,
    optional het-on-field). For sigma=0 make_truth_vcf ignores them (uniform)."""
    corr, sigma, cap = parse_field(w.field)
    o = ["--divergence-corr-len %s" % corr,
         "--divergence-sigma %s" % sigma,
         "--divergence-field-seed %d" % FIELD_SEED]
    if cap:
        o.append("--divergence-cap %s" % cap)
    if HET_ON_FIELD:
        o.append("--het-on-field")
    return " ".join(o)


wildcard_constraints:
    caller = "gatk|bcftools|freebayes",
    het = r"[0-9.]+", cov = r"[0-9]+", seed = r"[0-9]+",
    region = r"[A-Za-z0-9]+", iter = r"[0-9]+",
    field = r"L[0-9]+-S[0-9.]+(-C[0-9.]+)?"

## --- path templates ------------------------------------------------------
# div_target is inserted as a literal so the tree stays parallel to s01's, with the
# field tag added as the per-configuration subfolder (per GUIDE: group multi-knob runs).
region_fa  = os.path.join(outdir, "ref-regions", ref_str + "-{region}.fa")

sim_dir    = os.path.join(outdir, "sim", "{cov}X", "{region}", DIV, "{het}", "{field}", "{seed}")
read1      = os.path.join(sim_dir, ref_str + "_read1.fq.gz")
read2      = os.path.join(sim_dir, ref_str + "_read2.fq.gz")
truth_vcf  = os.path.join(sim_dir, ref_str + "-truth.vcf.gz")        # mason format
truth_norm = os.path.join(sim_dir, ref_str + "-truth.norm.vcf.gz")   # for scoring
truth_locald = os.path.join(sim_dir, ref_str + "-truth-locald.tsv")  # per-site local D

iter_dir   = os.path.join(outdir, "iter", "{caller}", "{cov}X", "{region}", DIV, "{het}", "{field}", "{seed}", "iter{iter}")
mapped_bam = os.path.join(iter_dir, ref_str + "_mapped.bam")
calls_norm = os.path.join(iter_dir, ref_str + "_calls.norm.vcf.gz")
cons_fa    = os.path.join(iter_dir, ref_str + "-consensus.fa")
recovered  = os.path.join(iter_dir, ref_str + "_recovered.vcf.gz")
summ       = os.path.join(iter_dir, ref_str + "-vcf-summary.csv")

combined   = os.path.join(resultsdir, "s05-" + TAG + "-iterative-summary.csv")

## --- recursive reference: what iteration i maps to / applies calls onto ----
def _ref_fa(w):
    i = int(w.iter)
    if i == 1:                                   # base reference (shared across callers)
        return os.path.join(outdir, "ref-regions", ref_str + "-" + w.region + ".fa")
    prev = os.path.join(outdir, "iter", w.caller, w.cov + "X", w.region, DIV, w.het,
                        w.field, w.seed, "iter" + str(i - 1))   # previous consensus
    return os.path.join(prev, ref_str + "-consensus.fa")

def ref_fa(w):    return _ref_fa(w)
def ref_bwt(w):   return _ref_fa(w) + ".bwt"
def ref_fai(w):   return _ref_fa(w) + ".fai"
def ref_dict(w):  return _ref_fa(w).replace(".fa", ".dict")

# the ORIGINAL reference (for the consensus-vs-reference diff), concrete per region
def base_ref(w):  return os.path.join(outdir, "ref-regions", ref_str + "-" + w.region + ".fa")

#############################################################################

localrules: all, combine_s05

rule all:
    input:
        combined,
        # always build the per-site local-divergence tables alongside the summary
        expand(truth_locald, cov=covs, region=[REGION], het=[HET], field=FIELDS, seed=SEEDS)

## --- 1. extract + index the reference region (base ref for iteration 1) ----
rule extract_region_ref:
    input:
        ref = REF
    output:
        fa    = region_fa,
        fai   = region_fa + ".fai",
        bwa   = multiext(region_fa, ".amb", ".ann", ".bwt", ".pac", ".sa"),
        dict_ = region_fa.replace(".fa", ".dict")
    resources:
        cpus = 1, mem = "8g", time = "2:00:00", partition = "shared"
    shell:
        """
        samtools faidx {input.ref} {wildcards.region} > {output.fa}
        samtools faidx {output.fa}
        bwa index {output.fa}
        rm -f {output.dict_}
        picard CreateSequenceDictionary R={output.fa} O={output.dict_}
        """

## --- 2. controlled truth VCF (SNP-only) with the divergence FIELD ----------
rule make_truth_vcf:
    input:
        fa = region_fa
    output:
        vcf = truth_vcf,
        tbi = truth_vcf + ".tbi"
    params:
        script = os.path.join(LIB, "make_truth_vcf.py"),
        div = DIV, idiv = IDIV, ihet = IHET,
        fieldopts = field_opts,
        plain = lambda w, output: output.vcf[:-3]
    resources:
        # The autocorrelated field is built by an FFT over the contig length
        # (61,420,004 = 2^2 * 53 * 289717). The large prime factor forces numpy's
        # FFT onto Bluestein's algorithm, which peaks ~13 GB; 24 GB gives headroom.
        # (A future optimization could pad the FFT to a 5-smooth length, but that
        # would change the field realization for a given seed, so it is deferred.)
        cpus = 1, mem = "24g", time = "1:00:00", partition = "shared"
    shell:
        """
        python {params.script} --reference {input.fa} \
            --divergence {params.div} --heterozygosity {wildcards.het} \
            --indel-divergence {params.idiv} --indel-heterozygosity {params.ihet} \
            {params.fieldopts} \
            --seed {wildcards.seed} --out {params.plain}
        bgzip -f {params.plain}
        tabix -p vcf {output.vcf}
        """

rule normalize_truth:
    input:
        vcf = truth_vcf, fa = region_fa
    output:
        vcf = truth_norm, tbi = truth_norm + ".tbi"
    resources:
        cpus = 1, mem = "8g", time = "1:00:00", partition = "shared"
    shell:
        """
        zcat {input.vcf} \
            | awk 'BEGIN{{OFS="\\t"}} /^#/{{print; next}} {{$9="GT"; print}}' \
            | bcftools norm -f {input.fa} -m- -Oz -o {output.vcf}
        tabix -p vcf {output.vcf}
        """

## --- 2b. per-site MEASURED local divergence (the downstream binning axis) ---
rule annotate_local_divergence:
    input:
        vcf = truth_vcf, tbi = truth_vcf + ".tbi"
    output:
        tsv = truth_locald
    params:
        script = os.path.join(LIB, "annotate_local_divergence.py"),
        windows = LOCAL_WINDOWS
    resources:
        cpus = 1, mem = "8g", time = "1:00:00", partition = "shared"
    shell:
        """
        python {params.script} --truth-vcf {input.vcf} \
            --windows {params.windows} --out {output.tsv}
        """

## --- 3. simulate reads (mason) --------------------------------------------
rule simulate_reads:
    input:
        fa = region_fa, fai = region_fa + ".fai",
        vcf = truth_vcf, tbi = truth_vcf + ".tbi"
    output:
        read1 = read1, read2 = read2
    params:
        readlen = READLEN, fragmean = FRAGMEAN
    log:
        os.path.join(sim_dir, "logs", ref_str + "-mason.log")
    resources:
        cpus = 8, mem = "24g", time = "16:00:00", partition = "shared"
    shell:
        """
        LEN=$(cut -f2 {input.fai} | head -1)
        NFRAG=$(( {wildcards.cov} * LEN / (2 * {params.readlen}) ))
        mason_simulator -ir {input.fa} -iv {input.vcf} -n $NFRAG \
            --illumina-read-length {params.readlen} --fragment-mean-size {params.fragmean} \
            --num-threads {resources.cpus} --seed {wildcards.seed} \
            -o {output.read1} -or {output.read2} --version-check 0 &> {log}
        """

## --- 4. map reads to the iteration reference ------------------------------
rule map_reads:
    input:
        fa = ref_fa, bwt = ref_bwt, read1 = read1, read2 = read2
    output:
        bam = mapped_bam, bai = mapped_bam + ".bai"
    params:
        rg = r"@RG\tID:{caller}-{cov}X-{region}-{field}-{seed}-i{iter}\tSM:" + ref_str + r"\tLB:lib1\tPL:ILLUMINA"
    resources:
        cpus = 8, mem = "16g", time = "12:00:00", partition = "shared"
    shell:
        """
        bwa mem -t {resources.cpus} -R '{params.rg}' {input.fa} {input.read1} {input.read2} 2>/dev/null \
            | samtools sort -@ {resources.cpus} -o {output.bam} -
        samtools index {output.bam}
        """

## --- 5. call variants (native output, normalized) -------------------------
rule call_variants:
    input:
        fa = ref_fa, fai = ref_fai, dict_ = ref_dict,   # dict only needed by gatk
        bam = mapped_bam, bai = mapped_bam + ".bai"
    output:
        vcf = calls_norm, tbi = calls_norm + ".tbi"
    params:
        minmq = MINMQ, scatter = GATK_SCATTER
    log:
        os.path.join(iter_dir, "logs", ref_str + "-call.log")
    resources:
        cpus = 8, mem = "48g", time = "24:00:00", partition = "shared"
    shell:
        r"""
        set -eo pipefail
        d=$(dirname {output.vcf})
        case {wildcards.caller} in
          gatk)
            mkdir -p $d/scatter
            LEN=$(cut -f2 {input.fai} | head -1)
            N={params.scatter}; CHUNK=$(( (LEN + N - 1) / N )); CONTIG={wildcards.region}
            pids=""
            for i in $(seq 1 $N); do
                s=$(( (i-1)*CHUNK + 1 )); e=$(( i*CHUNK )); [ $e -gt $LEN ] && e=$LEN
                gatk HaplotypeCaller -R {input.fa} -I {input.bam} -ERC GVCF \
                    -L $CONTIG:$s-$e -O $d/scatter/g_$i.vcf.gz --java-options -Xmx4g \
                    >> {log} 2>&1 &
                pids="$pids $!"
            done
            fail=0; for p in $pids; do wait $p || fail=1; done
            [ $fail -eq 0 ] || {{ echo "a scatter HaplotypeCaller failed" >> {log}; exit 1; }}
            gargs=""; for i in $(seq 1 $N); do gargs="$gargs -I $d/scatter/g_$i.vcf.gz"; done
            gatk GatherVcfs $gargs -O $d/g.vcf.gz >> {log} 2>&1
            tabix -p vcf $d/g.vcf.gz
            gatk GenotypeGVCFs -R {input.fa} -V $d/g.vcf.gz -O $d/raw.vcf.gz >> {log} 2>&1
            bcftools norm -f {input.fa} -m- --atomize -Oz -o {output.vcf} $d/raw.vcf.gz 2>> {log}
            rm -rf $d/scatter ;;
          bcftools)
            bcftools mpileup -f {input.fa} --min-MQ {params.minmq} -a FORMAT/AD,FORMAT/DP \
                --threads {resources.cpus} {input.bam} 2> {log} \
              | bcftools call -mv --threads {resources.cpus} 2>> {log} \
              | bcftools norm -f {input.fa} -m- --atomize -Oz -o {output.vcf} 2>> {log} ;;
          freebayes)
            freebayes -f {input.fa} --min-mapping-quality {params.minmq} {input.bam} 2> {log} \
              | bcftools norm -f {input.fa} -m- --atomize -Oz -o {output.vcf} 2>> {log} ;;
        esac
        tabix -p vcf {output.vcf}
        """

## --- 6. build the iteration consensus (reference for the next iteration) ----
rule build_consensus:
    input:
        fa = ref_fa, calls = calls_norm, tbi = calls_norm + ".tbi"
    output:
        fa    = cons_fa,
        fai   = cons_fa + ".fai",
        bwa   = multiext(cons_fa, ".amb", ".ann", ".bwt", ".pac", ".sa"),
        dict_ = cons_fa.replace(".fa", ".dict")
    resources:
        cpus = 1, mem = "8g", time = "2:00:00", partition = "shared"
    shell:
        """
        # Apply only HOMOZYGOUS SNP calls to the consensus (standard pseudo-it), same as
        # s01: SNP-only keeps the consensus the same length as the reference (identity
        # coordinates), and hom-only avoids IUPAC ambiguity codes that would break the
        # next iteration's mapping/calling.
        d=$(dirname {output.fa})
        bcftools view -v snps -g hom -Oz -o $d/calls.homsnps.vcf.gz {input.calls}
        tabix -f -p vcf $d/calls.homsnps.vcf.gz
        bcftools consensus -f {input.fa} -o {output.fa} $d/calls.homsnps.vcf.gz
        samtools faidx {output.fa}
        bwa index {output.fa}
        rm -f {output.dict_}
        picard CreateSequenceDictionary R={output.fa} O={output.dict_}
        """

## --- 7. score: consensus-vs-reference SNP diff, compared to truth ----------
rule score_iteration:
    input:
        base = base_ref,
        cons = cons_fa,
        truth = truth_norm,
        diff_script = os.path.join(SCRIPTS, "consensus_diff_snps.py"),
        cmp_script = os.path.join(LIB, "compare_vcfs.py")
    output:
        csv = summ,
        dict_ = summ.replace("-vcf-summary.csv", "-vcf-summary.dict.txt")
    params:
        recov_plain = lambda w, output: os.path.join(os.path.dirname(output.csv), ref_str + "_recovered.vcf"),
        prefix = lambda w, output: output.csv[:-len("-vcf-summary.csv")],
        div = DIV
    resources:
        cpus = 1, mem = "8g", time = "2:00:00", partition = "shared"
    shell:
        """
        python {input.diff_script} --reference {input.base} --consensus {input.cons} \
            --out {params.recov_plain}
        bgzip -f {params.recov_plain}
        tabix -p vcf {params.recov_plain}.gz
        python {input.cmp_script} --truth-vcf {input.truth} --called-vcf {params.recov_plain}.gz \
            --region {wildcards.region} --coverage {wildcards.cov} \
            --divergence {params.div} --heterozygosity {wildcards.het} --seed {wildcards.seed} \
            --caller {wildcards.caller} --iteration {wildcards.iter} \
            --extra field={wildcards.field} --out-prefix {params.prefix}
        """

## --- 8. combine all per-condition summaries -------------------------------
rule combine_s05:
    input:
        expand(summ, caller=CALLERS, cov=covs, region=[REGION], het=[HET],
               field=FIELDS, seed=SEEDS, iter=ITERS)
    output:
        combined
    params:
        script = os.path.join(LIB, "combine_summaries.py")
    shell:
        "python {params.script} --inputs {input} --out {output} --kind vcf"

#############################################################################
