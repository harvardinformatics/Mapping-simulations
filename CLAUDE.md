# Claude Code project instructions

**`AGENTS.md` in this directory is the binding rulebook for all work in this project. Follow it at
all times.** It is imported below so its rules are always in context — not only when a session
happens to read the file.

@AGENTS.md

The project charter — file/code organization priorities and data-analysis conventions — is
`GUIDE.md`, imported below.

@GUIDE.md

## Non-negotiable hard stops (inlined from AGENTS.md as a safeguard)

Ask first and wait for explicit confirmation before ANY of these — never do them silently:

- **Deleting any file** — no `rm`, `file.remove()`, `find -delete`, or scripted equivalents
  (including deletions embedded inside scripts). *"Never ever delete a file without asking permission."*
- `chmod`, package/dependency installation, network/API access, `git push`/`git pull`, or
  git-destructive commands (`git reset --hard`, `git clean -fd`, `git checkout --`).
- Modifying `AGENTS.md`, `GUIDE.md`, tests, or any non-script/non-config file.

Other standing rules: never insert silent fallbacks or hide errors (scripts should fail loudly);
never fabricate input data or outputs; put LLM-only temp/diagnostic files in `.llm-tmp/` (git-ignored),
not scattered in the project; accompany every plot/table/analysis with plain-text interpretation and a
data dictionary, and be explicit about binning/filtering in both the notebook and the plot. Use the
`mapping-sims` conda env (`renv` for R); if it won't activate, run via `conda run -n <env> ...`.
