#############################################################################
# s02 — indels + coordinate liftover (recovery scoring). All three callers.
#
# Iterative (pseudo-it) mapping WITH indels. Unlike s01 (SNP-only, identity coords), the
# consensus is length-shifted by indels, so iteration>1 lives in consensus coordinates.
#
# LIFTOVER ARCHITECTURE -- rebuild each consensus from the ORIGINAL reference using the
# ACCUMULATED reference-coordinate recovered calls (so every iteration has ONE clean
# ref<->consensus chain; no chain composition):
#
#   iter 1: map reads to the reference -> call -> the homozygous SNP+indel calls are
#           already reference coords -> that IS the accumulated recovered set acc_1.
#   iter k: map reads to consensus_{k-1} -> call (consensus_{k-1} coords) -> CrossMap the
#           calls back to reference via chain_rev_{k-1} (cons_{k-1}->ref) -> merge the new
#           homozygous calls into acc_{k-1} -> acc_k.
#   consensus_k = reference + acc_k  (bcftools consensus, emits ref->cons chain; chainSwap
#           -> cons->ref chain_rev_k used by iteration k+1).
#
# The accumulated reference-coordinate homozygous set acc_k IS the recovered divergence at
# iteration k -> scored directly vs truth with compare_vcfs (SNP and indel separately).
# Het recovery is scored from the lifted RAW caller VCF (het is never baked into the
# haploid consensus).
#
# Launch:  sbatch run_sweep.sh   (from this directory)
#############################################################################

import os

## --- config --------------------------------------------------------------
REF        = config["ref_file"]
ref_str    = config["ref_str"]
REGION     = config["region"]
SCRIPTS    = config["scripts_dir"]
LIB        = config["lib_dir"]
outdir     = config["data_dir"]
resultsdir = config["results_dir"]

covs   = config["covs"]
divs   = config["divs"]
HET    = config["het"]
SEEDS  = config["seeds"]
READLEN  = config["read_length"]
FRAGMEAN = config["fragment_mean"]
IRATIO = config["indel_snp_ratio"]
MINISZ = config["min_indel_size"]
MAXISZ = config["max_indel_size"]
IHET   = config["indel_heterozygosity"]
NITER  = config["n_iterations"]
CALLERS = config["callers"]
MINMQ  = config["min_mapping_quality"]
GATK_SCATTER = config["gatk_scatter"]
CHAINSWAP_ENV = config["chainswap_env"]
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
truth_vcf  = os.path.join(sim_dir, ref_str + "-truth.vcf.gz")       # mason format
truth_norm = os.path.join(sim_dir, ref_str + "-truth.norm.vcf.gz")  # GT-tagged, normalized
truth_oa   = os.path.join(sim_dir, ref_str + "-truth.oa.bam")       # mason -oa (ref coords)

iter_dir   = os.path.join(outdir, "iter", "{caller}", "{cov}X", "{region}", "{div}", "{het}", "{seed}", "iter{iter}")
mapped_bam = os.path.join(iter_dir, ref_str + "_mapped.bam")
calls_norm = os.path.join(iter_dir, ref_str + "_calls.norm.vcf.gz")       # consensus_{k-1} coords
lifted     = os.path.join(iter_dir, ref_str + "_calls.lifted.vcf.gz")     # reference coords
acc        = os.path.join(iter_dir, ref_str + "_recovered.acc.vcf.gz")    # accumulated recovered (ref)
cons_fa    = os.path.join(iter_dir, ref_str + "-consensus.fa")
chain_fwd  = os.path.join(iter_dir, "ref2cons.chain")
chain_rev  = os.path.join(iter_dir, "cons2ref.chain")
recov_summ = os.path.join(iter_dir, ref_str + "-recovery-vcf-summary.csv")
het_summ   = os.path.join(iter_dir, ref_str + "-het-vcf-summary.csv")

combined_recov = os.path.join(resultsdir, "s02-" + TAG + "-recovery-summary.csv")
combined_het   = os.path.join(resultsdir, "s02-" + TAG + "-het-summary.csv")

## --- recursion helpers ---------------------------------------------------
def _base_ref(region):
    return os.path.join(outdir, "ref-regions", ref_str + "-" + region + ".fa")

def _prev_dir(w):
    return os.path.join(outdir, "iter", w.caller, w.cov + "X", w.region, w.div, w.het,
                        w.seed, "iter" + str(int(w.iter) - 1))

# reference that iteration k maps to / applies calls onto: iter1 = base ref, else cons_{k-1}
def map_ref(w):
    return _base_ref(w.region) if int(w.iter) == 1 else os.path.join(_prev_dir(w), ref_str + "-consensus.fa")
def map_ref_fai(w):  return map_ref(w) + ".fai"
def map_ref_bwt(w):  return map_ref(w) + ".bwt"
def map_ref_dict(w): return map_ref(w).replace(".fa", ".dict")

# chain to lift THIS iteration's calls/reads back to reference: iter1 none (dummy = base
# ref, ignored in-shell), else cons_{k-1}->ref chain from the previous build.
def prev_chain(w):
    return _base_ref(w.region) if int(w.iter) == 1 else os.path.join(_prev_dir(w), "cons2ref.chain")

# accumulated recovered set to extend: iter1 none (dummy = base ref, ignored), else acc_{k-1}
def prev_acc(w):
    return _base_ref(w.region) if int(w.iter) == 1 else os.path.join(_prev_dir(w), ref_str + "_recovered.acc.vcf.gz")

def base_ref(w):      return _base_ref(w.region)
def base_ref_fai(w):  return _base_ref(w.region) + ".fai"
def base_ref_dict(w): return _base_ref(w.region).replace(".fa", ".dict")   # Picard target dict

localrules: all, combine_recovery, combine_het

rule all:
    input:
        combined_recov, combined_het

## --- 1. extract + index the reference region ------------------------------
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

## --- 2. controlled truth (SNPs + indels), normalized with a real GT --------
rule make_truth_vcf:
    input:
        fa = region_fa, fai = region_fa + ".fai"
    output:
        vcf = truth_vcf, tbi = truth_vcf + ".tbi"
    params:
        script = os.path.join(LIB, "make_truth_vcf.py"),
        het = HET, minisz = MINISZ, maxisz = MAXISZ, ihet = IHET,
        idiv = lambda w: "%.6f" % (IRATIO * float(w.div)),   # Di = ratio * D
        plain = lambda w, output: output.vcf[:-3]            # make_truth_vcf writes plain VCF
    resources:
        cpus = 1, mem = "8g", time = "1:00:00", partition = "shared"
    shell:
        """
        python {params.script} --reference {input.fa} \
            --divergence {wildcards.div} --heterozygosity {params.het} \
            --indel-divergence {params.idiv} --indel-heterozygosity {params.ihet} \
            --min-indel-size {params.minisz} --max-indel-size {params.maxisz} \
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
        # mason-format truth has FORMAT='.'; tag a real GT (col 9) so bcftools/compare_vcfs
        # can read the genotypes, then normalize (split multiallelics, atomize).
        r"""
        zcat {input.vcf} \
          | awk 'BEGIN{{OFS="\t"}} /^#/{{print; next}} {{$9="GT"; print}}' \
          | bcftools norm -f {input.fa} -m- --atomize -Oz -o {output.vcf}
        tabix -p vcf {output.vcf}
        """

## --- 3. simulate reads + mason -oa truth alignment (reference coords) -------
rule simulate_reads:
    input:
        fa = region_fa, fai = region_fa + ".fai",
        vcf = truth_vcf, tbi = truth_vcf + ".tbi"
    output:
        read1 = read1, read2 = read2, oa = truth_oa
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
        TMP=$(mktemp -d)
        mason_simulator -ir {input.fa} -iv {input.vcf} -n $NFRAG \
            --illumina-read-length {params.readlen} --fragment-mean-size {params.fragmean} \
            --num-threads {resources.cpus} --seed {wildcards.seed} \
            -o {output.read1} -or {output.read2} -oa $TMP/truth.sam --version-check 0 &> {log}
        samtools view -b $TMP/truth.sam -o {output.oa}
        rm -rf $TMP
        """

## --- 4. map reads to the iteration reference (cons_{k-1}; iter1 = base ref) --
rule map_reads:
    input:
        fa = map_ref, bwt = map_ref_bwt, read1 = read1, read2 = read2
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

## --- 5. call variants (native output, normalized; SNPs AND indels kept) -----
rule call_variants:
    input:
        fa = map_ref, fai = map_ref_fai, dict_ = map_ref_dict,
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

## --- 6. lift this iteration's calls consensus_{k-1} -> reference ------------
rule lift_calls:
    input:
        calls = calls_norm, tbi = calls_norm + ".tbi",
        chain = prev_chain, base = base_ref, base_fai = base_ref_fai, base_dict = base_ref_dict
    output:
        vcf = lifted, tbi = lifted + ".tbi"
    log:
        os.path.join(iter_dir, "logs", ref_str + "-lift.log")
    resources:
        cpus = 1, mem = "20g", time = "4:00:00", partition = "shared"
    shell:
        # iter1: calls are already reference coords -> copy. iter>1: Picard LiftoverVcf
        # through the cons_{k-1}->ref chain -- it re-derives indel alleles from the
        # reference (CrossMap mangled ~40% of indels into invalid records). ~1% of calls
        # are genuinely un-liftable (indels straddling a baked-indel boundary) and go to a
        # reject file. Re-normalize the lifted set against the reference.
        r"""
        if [ {wildcards.iter} -eq 1 ]; then
            cp {input.calls} {output.vcf}; tabix -f -p vcf {output.vcf}
        else
            d=$(dirname {output.vcf})
            rm -f $d/lifted.picard.vcf.gz $d/lifted.picard.vcf.gz.tbi $d/lift.reject.vcf.gz  # clear stale
            # freebayes emits up to ~2M records at high divergence; give Picard a big heap
            # and spill the sort to disk (low MAX_RECORDS_IN_RAM) so LiftoverVcf doesn't OOM.
            picard -Xmx16g LiftoverVcf I={input.calls} O=$d/lifted.picard.vcf.gz \
                CHAIN={input.chain} REJECT=$d/lift.reject.vcf.gz R={input.base} \
                MAX_RECORDS_IN_RAM=100000 WARN_ON_MISSING_CONTIG=true &> {log}
            bcftools norm -f {input.base} -m- --atomize -Oz -o {output.vcf} $d/lifted.picard.vcf.gz 2>> {log}
            tabix -f -p vcf {output.vcf}
            rm -f $d/lifted.picard.vcf.gz $d/lifted.picard.vcf.gz.tbi
        fi
        """

## --- 7. accumulate recovered (homozygous) calls in reference coords --------
rule accumulate:
    input:
        lifted = lifted, tbi = lifted + ".tbi",
        prev = prev_acc, base = base_ref, base_fai = base_ref_fai
    output:
        vcf = acc, tbi = acc + ".tbi"
    resources:
        cpus = 1, mem = "8g", time = "1:00:00", partition = "shared"
    shell:
        # keep only homozygous SNP+indel calls (fixed differences); iter1 = this round's
        # homs, iter>1 = previous accumulation UNION this round's homs (concat, sort, dedup).
        r"""
        d=$(dirname {output.vcf})
        bcftools view -g hom -Oz -o $d/this.hom.vcf.gz {input.lifted}
        tabix -f -p vcf $d/this.hom.vcf.gz
        if [ {wildcards.iter} -eq 1 ]; then
            bcftools norm -f {input.base} -m- -Oz -o {output.vcf} $d/this.hom.vcf.gz
        else
            bcftools concat -a {input.prev} $d/this.hom.vcf.gz 2>/dev/null \
              | bcftools sort -Oz -o $d/acc.sorted.vcf.gz
            # dedup exact duplicates (both inputs are already split by lift_calls' norm, so
            # -m- is unnecessary and bcftools 1.17 forbids combining -d with -m-)
            bcftools norm -f {input.base} -d exact -Oz -o {output.vcf} $d/acc.sorted.vcf.gz
            rm -f $d/acc.sorted.vcf.gz
        fi
        tabix -f -p vcf {output.vcf}
        rm -f $d/this.hom.vcf.gz $d/this.hom.vcf.gz.tbi
        """

## --- 8. rebuild the consensus from the reference + accumulated recovered ----
rule build_consensus:
    input:
        base = base_ref, acc = acc, tbi = acc + ".tbi"
    output:
        fa    = cons_fa,
        fai   = cons_fa + ".fai",
        bwa   = multiext(cons_fa, ".amb", ".ann", ".bwt", ".pac", ".sa"),
        dict_ = cons_fa.replace(".fa", ".dict"),
        fwd   = chain_fwd,
        rev   = chain_rev
    resources:
        cpus = 1, mem = "8g", time = "2:00:00", partition = "shared"
    shell:
        # consensus = reference + all accumulated homozygous recovered variants; bcftools
        # emits the ref->cons chain, chainSwap flips it to cons->ref for the next iteration.
        # (chainSwap is available directly in the mapping-sims env.)
        """
        bcftools consensus -f {input.base} -c {output.fwd} -o {output.fa} {input.acc}
        samtools faidx {output.fa}
        bwa index {output.fa}
        rm -f {output.dict_}
        picard CreateSequenceDictionary R={output.fa} O={output.dict_}
        chainSwap {output.fwd} {output.rev}
        """

## --- 9. score: accumulated recovered vs truth (divergence, SNP + indel) -----
rule score_recovery:
    input:
        acc = acc, tbi = acc + ".tbi", truth = truth_norm,
        cmp = os.path.join(LIB, "compare_vcfs.py")
    output:
        csv = recov_summ,
        dict_ = recov_summ.replace("-vcf-summary.csv", "-vcf-summary.dict.txt")
    params:
        prefix = lambda w, output: output.csv[:-len("-vcf-summary.csv")]
    resources:
        cpus = 1, mem = "8g", time = "1:00:00", partition = "shared"
    shell:
        """
        python {input.cmp} --truth-vcf {input.truth} --called-vcf {input.acc} \
            --region {wildcards.region} --coverage {wildcards.cov} \
            --divergence {wildcards.div} --heterozygosity {wildcards.het} --seed {wildcards.seed} \
            --caller {wildcards.caller} --iteration {wildcards.iter} --out-prefix {params.prefix}
        """

## --- 10. score: lifted raw caller VCF vs truth (het recovery) ---------------
rule score_het:
    input:
        lifted = lifted, tbi = lifted + ".tbi", truth = truth_norm,
        cmp = os.path.join(LIB, "compare_vcfs.py")
    output:
        csv = het_summ,
        dict_ = het_summ.replace("-vcf-summary.csv", "-vcf-summary.dict.txt")
    params:
        prefix = lambda w, output: output.csv[:-len("-vcf-summary.csv")]
    resources:
        cpus = 1, mem = "8g", time = "1:00:00", partition = "shared"
    shell:
        """
        python {input.cmp} --truth-vcf {input.truth} --called-vcf {input.lifted} \
            --region {wildcards.region} --coverage {wildcards.cov} \
            --divergence {wildcards.div} --heterozygosity {wildcards.het} --seed {wildcards.seed} \
            --caller {wildcards.caller} --iteration {wildcards.iter} --out-prefix {params.prefix}
        """

## --- 11. combine ----------------------------------------------------------
rule combine_recovery:
    input:
        expand(recov_summ, caller=CALLERS, cov=covs, region=[REGION], div=divs,
               het=[HET], seed=SEEDS, iter=ITERS)
    output:
        combined_recov
    params:
        script = os.path.join(LIB, "combine_summaries.py")
    shell:
        "python {params.script} --inputs {input} --out {output} --kind vcf"

rule combine_het:
    input:
        expand(het_summ, caller=CALLERS, cov=covs, region=[REGION], div=divs,
               het=[HET], seed=SEEDS, iter=ITERS)
    output:
        combined_het
    params:
        script = os.path.join(LIB, "combine_summaries.py")
    shell:
        "python {params.script} --inputs {input} --out {output} --kind vcf"
