#############################################################################
# s01 — baseline vs iterative mapping (SNP-only), all three callers.
#
# Iterative (pseudo-it) mapping: iteration 1 maps reads to the reference (this IS
# the single-shot baseline); iteration i>1 maps to the consensus built from the
# previous iteration's calls. Because the baseline is SNP-only, every consensus is
# the SAME LENGTH as the reference and shares its coordinates (verified: the
# bcftools consensus chain is a single gapless block) -- so NO liftover is needed.
#
# At each iteration the pipeline's recovered divergence = where the consensus differs
# from the reference (an exact position-wise SNP diff), scored vs truth with
# compare_vcfs. Each caller iterates independently (the consensus depends on its calls).
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
SCRIPTS    = config["scripts_dir"]     # s01-specific
LIB        = config["lib_dir"]         # shared analysis scripts
outdir     = config["data_dir"]
resultsdir = config["results_dir"]

covs   = config["covs"]
divs   = config["divs"]
HET    = config["het"]
SEEDS  = config["seeds"]
READLEN  = config["read_length"]
FRAGMEAN = config["fragment_mean"]
IDIV   = config["indel_divergence"]        # 0 for the SNP-only baseline
IHET   = config["indel_heterozygosity"]
NITER  = config["n_iterations"]
CALLERS = config["callers"]
MINMQ  = config["min_mapping_quality"]
GATK_SCATTER = config["gatk_scatter"]
TAG    = config["run_tag"]

ITERS = [str(i) for i in range(1, NITER + 1)]

wildcard_constraints:
    caller = "gatk|bcftools|freebayes",
    div = r"[0-9.]+", het = r"[0-9.]+", cov = r"[0-9]+", seed = r"[0-9]+",
    region = r"[A-Za-z0-9]+", iter = r"[0-9]+"

## --- path templates ------------------------------------------------------
region_fa  = os.path.join(outdir, "ref-regions", ref_str + "-{region}.fa")

sim_dir    = os.path.join(outdir, "sim", "{cov}X", "{region}", "{div}", "{het}", "{seed}")
read1      = os.path.join(sim_dir, ref_str + "_read1.fq.gz")
read2      = os.path.join(sim_dir, ref_str + "_read2.fq.gz")
truth_vcf  = os.path.join(sim_dir, ref_str + "-truth.vcf.gz")        # mason format
truth_norm = os.path.join(sim_dir, ref_str + "-truth.norm.vcf.gz")   # for scoring

iter_dir   = os.path.join(outdir, "iter", "{caller}", "{cov}X", "{region}", "{div}", "{het}", "{seed}", "iter{iter}")
mapped_bam = os.path.join(iter_dir, ref_str + "_mapped.bam")
calls_norm = os.path.join(iter_dir, ref_str + "_calls.norm.vcf.gz")
cons_fa    = os.path.join(iter_dir, ref_str + "-consensus.fa")
recovered  = os.path.join(iter_dir, ref_str + "_recovered.vcf.gz")
summ       = os.path.join(iter_dir, ref_str + "-vcf-summary.csv")

combined   = os.path.join(resultsdir, "s01-" + TAG + "-iterative-summary.csv")

## --- recursive reference: what iteration i maps to / applies calls onto ----
def _ref_fa(w):
    i = int(w.iter)
    if i == 1:                                   # base reference (shared across callers)
        return os.path.join(outdir, "ref-regions", ref_str + "-" + w.region + ".fa")
    prev = os.path.join(outdir, "iter", w.caller, w.cov + "X", w.region, w.div, w.het,
                        w.seed, "iter" + str(i - 1))   # previous iteration's consensus
    return os.path.join(prev, ref_str + "-consensus.fa")

def ref_fa(w):    return _ref_fa(w)
def ref_bwt(w):   return _ref_fa(w) + ".bwt"
def ref_fai(w):   return _ref_fa(w) + ".fai"
def ref_dict(w):  return _ref_fa(w).replace(".fa", ".dict")

# the ORIGINAL reference (for the consensus-vs-reference diff), concrete per region
def base_ref(w):  return os.path.join(outdir, "ref-regions", ref_str + "-" + w.region + ".fa")

#############################################################################

localrules: all, combine_s01

rule all:
    input:
        combined

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

## --- 2. controlled truth VCF (SNP-only) + normalized copy ------------------
rule make_truth_vcf:
    input:
        fa = region_fa
    output:
        vcf = truth_vcf,
        tbi = truth_vcf + ".tbi"
    params:
        script = os.path.join(LIB, "make_truth_vcf.py"),
        idiv = IDIV, ihet = IHET, plain = lambda w, output: output.vcf[:-3]
    resources:
        cpus = 1, mem = "8g", time = "1:00:00", partition = "shared"
    shell:
        """
        python {params.script} --reference {input.fa} \
            --divergence {wildcards.div} --heterozygosity {wildcards.het} \
            --indel-divergence {params.idiv} --indel-heterozygosity {params.ihet} \
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
        rg = r"@RG\tID:{caller}-{cov}X-{region}-{div}-{het}-{seed}-i{iter}\tSM:" + ref_str + r"\tLB:lib1\tPL:ILLUMINA"
    resources:
        cpus = 8, mem = "16g", time = "12:00:00", partition = "shared"
    shell:
        """
        bwa mem -t {resources.cpus} -R '{params.rg}' {input.fa} {input.read1} {input.read2} 2>/dev/null \
            | samtools sort -@ {resources.cpus} -o {output.bam} -
        samtools index {output.bam}
        """

## --- 5. call variants (native output, normalized). One rule, dispatched on the
##        {caller} wildcard so the recursive per-caller reference resolves. --------
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
            # scatter HaplotypeCaller across N equal intervals of the contig, run in
            # parallel on this node, gather the GVCFs, then a single GenotypeGVCFs.
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
        # Apply only HOMOZYGOUS SNP calls to the consensus (standard pseudo-it):
        #  - SNP-only: the baseline truth has no indels, so caller indels are false
        #    positives; dropping them keeps the consensus the SAME LENGTH as the
        #    reference (identity coordinates, no liftover). (s02 will apply real indels
        #    and use liftover.)
        #  - homozygous only: het sites are polymorphisms, not fixed differences, and
        #    bcftools consensus renders them as IUPAC ambiguity codes (e.g. M), which are
        #    non-ACGT and break mapping/calling in the next iteration. We update the
        #    reference only with fixed (1/1) differences, keeping a clean ACGT consensus.
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
        prefix = lambda w, output: output.csv[:-len("-vcf-summary.csv")]
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
            --divergence {wildcards.div} --heterozygosity {wildcards.het} --seed {wildcards.seed} \
            --caller {wildcards.caller} --iteration {wildcards.iter} --out-prefix {params.prefix}
        """

## --- 8. combine all per-condition summaries -------------------------------
rule combine_s01:
    input:
        expand(summ, caller=CALLERS, cov=covs, region=[REGION], div=divs, het=[HET],
               seed=SEEDS, iter=ITERS)
    output:
        combined
    params:
        script = os.path.join(LIB, "combine_summaries.py")
    shell:
        "python {params.script} --inputs {input} --out {output} --kind vcf"

#############################################################################
