# provenance.r -- R twin of analyses/lib/provenance.py.
#
# 1. prov_write_dictionary(): for R scripts that write tables. Writes <table>.dict.txt starting with
#    the same provenance block as the Python tools (Produced by / Command / Git commit / Written),
#    followed by the caller's field definitions.
# 2. render_provenance(): for notebooks. Prints which git commit a notebook was rendered from and
#    whether the notebook or any analysis code had uncommitted edits at render time.
#
# Run everything from the project root (as the notebooks and scripts already require).
# Git failures stop with an error (no silent fallback): provenance must never be guessed.

# --- git helpers ---------------------------------------------------------------------------------

prov_git <- function(...) {
  # run a git command in the project repo; stop() if git fails. --no-optional-locks: never take the
  # index lock (other jobs may be running git at the same time)
  out <- suppressWarnings(system2("git", c("--no-optional-locks", ...), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) {
    stop("git ", paste(c(...), collapse = " "), " failed: ", paste(out, collapse = "\n"))
  }
  out
}

prov_git_commit <- function() prov_git("rev-parse", "--short=12", "HEAD")

prov_script_state <- function(script) {
  # how a file relates to HEAD: committed, modified, or untracked
  status <- prov_git("status", "--porcelain", "--", script)
  if (length(status) == 0) return("script as committed")
  if (startsWith(status[1], "??")) return("script NOT YET COMMITTED (untracked)")
  "script has UNCOMMITTED EDITS relative to this commit"
}

prov_now <- function() format(Sys.time(), "%Y-%m-%d %H:%M")

prov_rscript_command <- function() {
  # the running Rscript's own command line, as typed from the project root
  args <- commandArgs(trailingOnly = FALSE)
  file <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(file) != 1) stop("prov_rscript_command(): not running under Rscript <file>")
  paste(c("Rscript", file, commandArgs(trailingOnly = TRUE)), collapse = " ")
}

# --- data dictionaries ---------------------------------------------------------------------------

prov_header <- function(table, script, command, backfill = FALSE, notes = character()) {
  # the 'Data dictionary for ...' title plus the provenance block (same format as provenance.py)
  lines <- c(sprintf("Data dictionary for %s", table), "",
             sprintf("Produced by : %s", script),
             sprintf("Command     : %s", command))
  if (backfill) {
    if (!file.exists(table)) stop("backfill: table ", table, " does not exist")
    mtime <- format(file.mtime(table), "%Y-%m-%d %H:%M")
    lines <- c(lines,
               sprintf("Git commit  : not recorded (table written %s, before provenance tracking)", mtime),
               sprintf("Dictionary  : added %s at commit %s (backfill; the table was not regenerated)",
                       prov_now(), prov_git_commit()))
  } else {
    lines <- c(lines,
               sprintf("Git commit  : %s (%s)", prov_git_commit(), prov_script_state(script)),
               sprintf("Written     : %s", prov_now()))
  }
  if (length(notes) > 0) lines <- c(lines, sprintf("Note        : %s", notes))
  paste0(paste(lines, collapse = "\n"), "\n\n")
}

prov_write_dictionary <- function(table, body, script, command = prov_rscript_command(),
                                  backfill = FALSE, notes = character()) {
  # write <table>.dict.txt = provenance block + body (character vector of lines, or one string)
  path <- paste0(table, ".dict.txt")
  writeLines(c(sub("\n$", "", prov_header(table, script, command, backfill, notes)), body), path)
  invisible(path)
}

# --- notebooks -----------------------------------------------------------------------------------

render_provenance <- function(notebook) {
  # Markdown lines for a notebook's provenance box: the commit it was rendered from, and whether the
  # notebook itself or any tracked analysis code had uncommitted edits at render time (if so, the
  # rendered output may not match that commit exactly).
  commit <- prov_git_commit()
  branch <- prov_git("rev-parse", "--abbrev-ref", "HEAD")
  nb_state <- prov_script_state(notebook)
  dirty <- prov_git("status", "--porcelain", "--untracked-files=no", "--", "analyses")
  dirty_files <- sub("^...", "", dirty)
  lines <- c(sprintf("- **Rendered from commit** `%s` (branch `%s`) on %s", commit, branch, prov_now()),
             sprintf("- **This notebook:** %s", sub("^script", "notebook", nb_state)))
  if (length(dirty_files) == 0) {
    lines <- c(lines, "- **Analysis code:** no uncommitted edits under `analyses/`")
  } else {
    shown <- head(dirty_files, 8)                     # keep the box short; git status has the rest
    more <- if (length(dirty_files) > 8) sprintf(" and %d more", length(dirty_files) - 8) else ""
    lines <- c(lines, sprintf("- **Uncommitted edits under `analyses/` at render time (%d files):** %s%s",
                              length(dirty_files), paste0("`", shown, "`", collapse = ", "), more))
  }
  c(lines, "- Each data file's producing script and command: the section's `MANIFEST.yaml`")
}
