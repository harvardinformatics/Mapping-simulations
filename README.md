# Mapping-simulations

Simulations of **reference bias**: how divergence between a sample and the reference genome causes
true variants to be lost during read mapping and variant calling, and whether **iterative mapping**
(pseudo-it style: call, bake the calls into a consensus, remap) recovers them. Reads are simulated with
`mason2` from controlled truth variants on mouse **chr19** (mm39), mapped with `bwa mem`, and called
with three callers (GATK HaplotypeCaller, bcftools, freebayes). The project charter is
[`GUIDE.md`](GUIDE.md).

## Repository layout

| path | what it holds |
|---|---|
| `analyses/<section>/` | The current analyses, one section per question (table below). Each has its workflow and run scripts, `config/`, `scripts/`, `results/` (tables with data dictionaries), a notebook (`.rmd` + rendered `.html`) and a `MANIFEST.yaml`. |
| `analyses/lib/` | Shared tools: truth-VCF generation, divergence fields, scoring of calls and read placements, the per-variant read-fate ledger, provenance helpers, the manifest checker and the re-run planner. |
| `data/` | All simulation work files (references, truth VCFs, reads, BAMs, calls, consensus sequences). **Not in git.** |
| `docs/` | HTML pages from earlier phases. `docs/mason_results.html` is the first, exploratory mason2 pass. |
| `scripts/`, `snakemake/`, `simulation-configs/`, `summary-data/` | The exploratory mason2 pipeline behind `docs/mason_results.html` (`snakemake/simulate_and_map.smk`). The analyses in `analyses/` do not use it. |
| `mapping-sims.yml`, `renv.yml` | Conda environments for the pipeline and the notebooks (see [Software](#software)). |
| `*/v1/` | The archived first version, based on the NEAT read simulator, with its environment file in `snakemake/v1/`. Its run instructions are in this README's git history. |
| `evolution/` | 2023 Evolution meeting talk and notes. |
| `pseudo-it/` | Separate git repository (pseudo-it 4 development), ignored by this one. |
| `AGENTS.md`, `CLAUDE.md` | Rules for LLM agents working in this repository. |

## Analysis sections

| section | question | notebook |
|---|---|---|
| `s01-baseline-iterative` | Does iterative mapping recover more true divergence than single-shot mapping, for each caller? SNPs only, uniform divergence. | `s01-results.rmd`, `s01-variant-fate-panel.rmd` |
| `s02-indels` | The same with indels added, using consensus liftover back to reference coordinates. | `s02-results.rmd` |
| `s03-heterozygosity` | Does divergence bias heterozygous sites toward the reference allele, by what read-level mechanism, and does iteration restore the balance? | `s03-results.rmd` |
| `s04-softclip` | Does disabling soft-clipping (`--dont-use-soft-clipped-bases`) change GATK's recall? | `s04-results.rmd` |
| `s05-heterogeneous-divergence` | At a fixed genome-wide mean, how does patchy (spatially autocorrelated) divergence change recovery, and is recall predicted by local divergence? | `s05-results.rmd` |
| `s06-filtering` | Can hard filters (QUAL, DP, MQ) recover the lost divergence, or only buy precision? | `s06-results.rmd` |
| `s07-svs` | How do structural variants, which no caller here calls, depress small-variant recovery near their breakpoints? | `s07-results.rmd` |
| `s08-mappability` | How much does low reference mappability (repeats) reduce recovery, independent of divergence? | `s08-results.rmd` |
| `synthesis` | All sections together, plus a short figure set for coauthors. | `synthesis.rmd`, `brief.rmd` |

## How the simulations are made

Each section's Snakemake workflow (or run scripts, for sections that reuse earlier data) does:

1. **Reference:** extract chr19 from mm39 (GRCm39).
2. **Truth variants** (`analyses/lib/make_truth_vcf.py`): fixed differences at divergence rate D
   (genotype 1|1, random alternative alleles) plus heterozygous SNPs at 0.5% (0|1). Depending on
   the section, also indels (1-10 bp, half insertions, indel:SNP about 0.12), structural variants
   (s07), or a spatially varying divergence field with its mean fixed at D (s05).
3. **Reads:** `mason_simulator` applies the truth VCF and simulates 150 bp paired-end reads (fragment
   mean 300 bp, 30x, seed 42; s01 also seeds 43 and 44). Its `-oa` output records each read's true
   origin, which the read-level analyses use as ground truth.
4. **Mapping and iteration:** `bwa mem` to the reference, then 3 iterations. Each calls variants,
   writes them into the reference with `bcftools consensus`, and remaps to that consensus.
   Iteration 1 is ordinary single-shot mapping. With indels (s02 on), results are lifted back to
   reference coordinates for scoring (Picard for VCFs, a point-lift for BAMs).
5. **Calling:** GATK HaplotypeCaller (GVCF mode, then GenotypeGVCFs), bcftools mpileup/call, and
   freebayes. Reads below mapping quality 20 are ignored by all three: bcftools and freebayes get
   `min_mapping_quality: 20` from the config, and 20 is GATK's default read filter.
6. **Scoring:** `compare_vcfs.py` (recall of truth variants), `compare_bams.py` (read placement vs
   truth) and `variant_read_fate.py` (for each variant, where the reads that truly cover it went).

Run scripts are submitted to Slurm from inside their section directory, because their log paths are
relative, for example:

```bash
cd analyses/s01-baseline-iterative && sbatch run_sweep.sh
```

Use `--dryrun` on any Snakemake command before running it for real. Every table's exact producing
command is in its section's `MANIFEST.yaml` (see below).

## Software

Two conda environments:

- **`mapping-sims`** for the pipeline, from `mapping-sims.yml`: mason 2.0.9, bwa 0.7.17, samtools
  and bcftools 1.17, GATK 4.3.0.0, freebayes 1.3.6, Picard 2.27.5, GenMap 1.3.0, snakemake 7.20.0,
  Python 3.9 with pysam 0.21.0 and numpy 1.24. The file is an export of the working environment, so
  every package is pinned.
- **`renv`** for the notebooks, from `renv.yml`: R 4.5.3 with rmarkdown, knitr, ggplot2, dplyr,
  tidyr, readr, patchwork, ggalluvial, ggrepel and kableExtra, and pandoc 3.8.2. It lists only what
  the notebooks load, at the versions in use; creating an environment from it has not been tested yet.

```bash
conda env create -f mapping-sims.yml
conda env create -f renv.yml
```

The v1 pipeline's 2021 environment file is kept at `snakemake/v1/mapping-sims.yml`.

Render a notebook from the project root:

```bash
conda run -n renv Rscript -e 'rmarkdown::render("analyses/s04-softclip/s04-results.rmd", knit_root_dir = getwd(), quiet = TRUE)'
```

## Provenance and data dictionaries

Three things record where every result in `analyses/` came from:

- **Data dictionaries.** Every table has a `<table>.dict.txt` next to it. It defines the columns and
  starts with a provenance block: the producing script, the exact command, and the git commit the
  table was written at. Tables made before 2026-10-26 have dictionaries marked "backfill": they were
  documented afterwards, with the commands recovered from the run scripts and session logs, and
  `analyses/lib/backfill/` records how.
- **Section manifests.** `analyses/<section>/MANIFEST.yaml` is the section's top-level data
  dictionary. For every data file, figure and notebook it lists what it is, the script and exact
  command that produce it, its inputs, the notebooks that read it, whether git tracks it, and its
  status: `current`, `earlier-version`, `artifact` or `unused`, with a note explaining any status
  other than `current`. Check all of them with:

  ```bash
  conda run -n mapping-sims python analyses/lib/check_manifest.py
  ```

  It fails if a file is missing from its manifest, a table has no dictionary, or a notebook reads a
  file that no manifest lists. Run it after adding or regenerating outputs.
- **Notebooks.** Each rendered notebook starts with a box giving the commit it was rendered from and
  any uncommitted edits at render time.

Large per-site intermediate tables are not in git (`git: ignored` in the manifests). Before
re-rendering a notebook that reads one, regenerate it from `data/` with the command in its
manifest entry.

## Re-running after a change

The manifests together form a dependency graph, and `analyses/lib/rerun.py` walks it to list what a
change requires, in run order. It only prints a plan; it runs nothing.

```bash
# everything downstream of a changed file (code, config or input)
conda run -n mapping-sims python analyses/lib/rerun.py plan --changed analyses/lib/truth_haplotypes.py

# outputs whose code or inputs changed since the commit their dictionary (or provenance box) records
conda run -n mapping-sims python analyses/lib/rerun.py plan --stale
```

For each step the plan gives the reason, the exact command, whether it runs on Slurm or locally, and
warnings: steps that delete files (for example Snakemake `temp()` outputs), and steps where Snakemake
will not see the change and needs `--forcerun`. Dependencies are tracked per file, so a change
anywhere in a file lists everything that uses it; read the plan and strike steps the change cannot
affect. Tables whose dictionaries predate provenance tracking can only be judged through their inputs
by `--stale`.

To keep provenance clean when acting on a plan:

1. Commit the code change first, so every regenerated dictionary records a real commit.
2. Run the listed commands stage by stage, waiting for each stage's Slurm jobs to finish.
3. Check the manifests (`check_manifest.py`) and commit the regenerated outputs.
4. Render the affected notebooks from that clean commit, the brief last (rendering it rewrites its
   exported figures, which would otherwise show as uncommitted edits in the other notebooks'
   provenance boxes), and commit the HTML.
