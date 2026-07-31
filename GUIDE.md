# Reference bias and iterative mapping simulations

Author: Gregg Thomas
Date: July 29, 2026

This document exists to outline the goals and rules for the analysis of the mapping simulation project hosted in the project repo. It is intended both for human readers to understand how the analyses were organized and performed and to provide guides and guardrails for the LLMs performing the analyses.

Many of these analyses have already been performed but in an exploratory and disorganized way. For reproducibility, both now and later, these will be redone here in a more deliberate fashion.

## Overarching guidelines

- The project's git root is `/n/holylfs05/LABS/informatics/Lab/projects/gwct/Mapping-simulations/`, which is referred to in this document as **project root**.
- LLMs reading this must abide by the guides set in the project root's `AGENTS.md` file.
- An LLM must never edit this file, `GUIDE.md`, or the aforementioned `AGENTS.md` file.
- LLMs tend to create many one-off or temporary files. This should be kept to a minimum in this folder because all analyses and scripts should be deliberate and focused. Temporary files used by the LLM must go in `.llm-tmp/`, as per AGENTS.md. Any other files that must be retained but are not required for the analyses should go in `etc/`.

### Metadata

For every saved data table or data structure, a data dictionary must be created. A data dictionary is a document that outlines the structure, content, and variable definitions for a dataset or collection of data.

See [https://datamanagement.hms.harvard.edu/collect-analyze/documentation-metadata/data-dictionary](here) for more guidance and links to examples.

Each table or data file should have its own accompanying dictionary, and each section should have a top-level dictionary that lists all tables or data in that section.

Likewise, a README should be created and maintained in the project folder to outline the structure of the folder and the data processing and analyses.

### Scripts

- All scripts must be written to be completely understandable by a human reader for reproducibility. This means an organized structure and spacing and, importantly, comments that explain each step of the analysis. Not necessarily one comment per line, but there should be no shortage of explanatory comment blocks and explanatory line comments.

### Notebooks

For each sub-section of the analysis (see below), the LLM must create an accompanying notebook with plots and text explanations. These notebooks are not meant to be products for publication, but rather for context for human reviewers of the code and analyses.

- The text in these notebooks should be concise and understandable.
- The text should describe the methods used with enough detail to be entirely reproducible without the use of an LLM by a human reading the notebook, however they are not intended to be full Methods sections. Bullet points will suffice.
- The notebook layout should follow this notebook: `docs/scripts/markdown/mason_results.rmd`. The content should be wide, with a sticky toc on the left. Notably, LLMs have struggled with this because it is accomplished via two style blocks, one in the header and one in the markdown text. Both must be included.

### Assumptions and thresholds

During the course of any analysis, assumptions must be made. These must be clearly documented in any notebooks, metadata, and scripts (as comments).

Likewise, certain decisions regarding thresholds for filtering, statistical cutoffs, figure display, etc. must be clearly documented and, whenever possible, available as options to tweak easily in the scripts. These could live in a single block at the top of the script.

Resulting filenames and/or metadata should be labeled with the thresholds that vary. Re-running with a new threshold should never overwrite another threshold's result. The complete parameter set (including thresholds not in the name, plus seeds) is always recorded in that output's accompanying metadata/dictionary. When many thresholds vary at once, group the run under a short per-configuration subfolder rather than a long filename.

## Analysis outline

The overall goal of the project is to determine how divergence affects read mapping and SNP calling, and to test whether iterative mapping recovers more true divergence.

The completed analyses thus far include pilot simulations of chr19 of the mouse genome and subsequent analysis, including analysis of different variant callers.

### Next priorities

These steps should have self contained results notebooks, but there must be a main notebook that compares results between them and between the completed analysis.

The main priority is to determine if iterative mapping improves results. I am not sure whether this would be easiest to do first, or to address the things below first.

### Simulation effects

In addition to divergence, we must compare the effects of various other parameters:

1. Indels
2. Heterozygosity
3. Structural variation
4. Heterogenous divergence

I am not sure the best way to do this: we could have one baseline simulation without indels or SVs, homogenous divergence, and one value of heterozygosity, then do comparisons off of that while adding the other parameters, but that doesn't really get at how the parameters interact.

### Mapping

Does disabling soft clipping affect results?

### Variant filtering

How does a sweep of different filtering criteria (e.g. mapping quality, read depth, etc.) affect results?

### Low priority

1. Whole-genome mapping: I think just using chr19 for speed is best until we nail down the pipeline and simulation parameters
2. Comparison to v1: I'm curious, but ultimately it doesn't matter if we've found something that works now. 

