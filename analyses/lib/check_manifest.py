#!/usr/bin/env python3
"""
check_manifest.py -- check that every data file, figure and notebook in each analysis section is
documented in that section's MANIFEST.yaml, and that the manifests agree with the files on disk,
with git, and with what the notebooks actually read.

WHY: each MANIFEST.yaml is the section's top-level data dictionary (GUIDE.md) AND its provenance
record: for every output it names the script and the exact command that produce it, the inputs, and
the notebooks that use it. This checker keeps that record honest as files are added or changed.

Usage (from the project root):
    conda run -n mapping-sims python analyses/lib/check_manifest.py              # all sections
    conda run -n mapping-sims python analyses/lib/check_manifest.py s04-softclip # some sections
Exit status 1 if any ERROR is found; NOTICEs alone exit 0.

---------------------------------------------------------------------------------------------------
MANIFEST.yaml SCHEMA (all paths relative to the project root)

section:  s04-softclip                      # must equal the directory name
question: one line -- what this section asks
upstream:                                    # optional: data this section reads from elsewhere
  - path: data/s01-baseline-iterative/       #   (typically the git-ignored data/ tree)
    from: how it is made (e.g. "s01 workflow: bash analyses/s01-baseline-iterative/run_sweep.sh")
notebooks:
  - path:    analyses/s04-softclip/s04-results.rmd
    html:    analyses/s04-softclip/s04-results.html
    render:  exact render command
outputs:                                     # one entry per file, or per family of files
  - path:    analyses/s04-softclip/results/softclip-toggle-{D}.csv   # {name} placeholders ...
    values:  {D: ["0.06", "0.10"]}                                  # ... expanded from values
    # glob: analyses/x/results/fig-*.png    (instead of path, for irregular families: then a
    #                                        missing member cannot be detected)
    what:    one line: what the file holds
    kind:    table | figure | caption | other   (default from the extension: csv/tsv/vcf/bed/json
                                                 -> table, png/pdf/svg -> figure, else other)
    dictionary: template of its data dictionary (default: <path>.dict.txt; tables only)
    script:  the script that writes it (an existing file in the repo), or `none` when the producing
             code no longer exists (only with a non-current status and a note saying so), or
             `hand-written` for a file no script writes, e.g. a caption (with a note)
    command: the exact command that produced it ({name} placeholders allowed)
    inputs:  [list of input paths or templates]                     # optional
    used_by: [notebooks/scripts that read it; [] if none]
    git:     tracked | ignored        # ignored = too large / regenerable; NOT committed
    status:  current | earlier-version | artifact | unused        # default current
    note:    free text; REQUIRED when status is not current
ignore:                                      # optional: files that are not outputs, with a reason
  - glob:   analyses/s02-indels/results/*.status
    reason: job-status marker written by an ad hoc background loop

status meanings:
  current          made by the listed script/command as the code stands now
  earlier-version  made by an earlier version of the script (note: what differs / was verified)
  artifact         known to be wrong; kept for the record, do not use (note: why)
  unused           nothing reads it (candidate for removal -- ask before deleting)

---------------------------------------------------------------------------------------------------
CHECKS
ERRORS (exit 1):
  UNLISTED        a file in the section matches no outputs/notebooks/ignore entry
  MISSING         an expanded `path` with git: tracked does not exist
  GIT             git state disagrees with the entry (says tracked but is ignored, or vice versa)
  NO-DICTIONARY   a table exists but its dictionary does not
  BAD-ENTRY       missing required field, unknown status/kind/git value, script/used_by not found
  UNDOCUMENTED-REF  a notebook reads (by literal path) a file no manifest lists
  USED_BY         a notebook reads a listed file, but the entry's used_by does not name it
NOTICES (exit 0):
  ABSENT          a git: ignored file is not present here (regenerate with its command)
  UNCOMMITTED     a file is neither committed nor git-ignored yet
  DYNAMIC-REF     a notebook path built with sprintf/glue placeholders (matched as a pattern)

Scope of "files in a section": everything under analyses/<section>/ except code (*.py *.R *.r
*.sh *.smk *.css), data dictionaries (*.dict.txt, checked as dictionaries), index files (*.tbi
*.csi *.fai), the MANIFEST.yaml itself, and the config/, slurm-logs/, .snakemake/ and __pycache__/
directories.
"""
import fnmatch
import glob
import itertools
import os
import re
import subprocess
import sys

import yaml

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
ANALYSES = "analyses"
NOT_SECTIONS = {"lib", "main-comparison"}           # analyses/ subdirectories without outputs

SKIP_DIRS = {"config", "slurm-logs", ".snakemake", "__pycache__"}
CODE_EXT = (".py", ".R", ".r", ".sh", ".smk", ".css")
INDEX_EXT = (".tbi", ".csi", ".fai")
TABLE_EXT = (".csv", ".tsv", ".vcf", ".vcf.gz", ".bed", ".bed.gz", ".json")
FIGURE_EXT = (".png", ".pdf", ".svg")

STATUSES = {"current", "earlier-version", "artifact", "unused"}
KINDS = {"table", "figure", "caption", "other"}
GITS = {"tracked", "ignored"}
REQUIRED = ("what", "script", "command", "used_by", "git")

# literal project paths in notebooks: "analyses/..." or 'data/...'
PATH_LITERAL = re.compile(r"""["']((?:analyses|data)/[^"'\s]+?)["']""")
# sprintf-style (%s, %d, %.2f, %05.1f ...) and glue-style ({x}) placeholders
PLACEHOLDER = re.compile(r"%[-+ 0#]*\d*(?:\.\d+)?[sdfgiex]|\{[^}]*\}")


# --- helpers -------------------------------------------------------------------------------------

def git_lines(*args):
    """Lines of a git command's output (run in the project root); raises if git fails."""
    out = subprocess.run(["git", "--no-optional-locks", "-C", ROOT] + list(args),
                         check=True, capture_output=True, text=True).stdout
    return set(l for l in out.splitlines() if l)


def default_kind(path):
    if path.endswith(TABLE_EXT):
        return "table"
    if path.endswith(FIGURE_EXT):
        return "figure"
    return "other"


def expand(template, values):
    """All paths from a {name} template and its values dict (cartesian product)."""
    names = sorted(set(re.findall(r"\{(\w+)\}", template)))
    if not names:
        return [template]
    missing = [n for n in names if n not in (values or {})]
    if missing:
        raise ValueError("no values for placeholder(s) %s" % ", ".join(missing))
    out = []
    for combo in itertools.product(*[[str(v) for v in values[n]] for n in names]):
        sub = dict(zip(names, combo))
        out.append(re.sub(r"\{(\w+)\}", lambda m: sub[m.group(1)], template))
    return out


def fill(template, path_template, path):
    """Fill a companion template (e.g. a dictionary path) with the placeholder values of `path`."""
    names = re.findall(r"\{(\w+)\}", path_template)
    if not names:
        return template
    rx = "^" + re.sub(r"\\\{(\w+)\\\}", r"(?P<\1>[^/]+?)", re.escape(path_template)) + "$"
    m = re.match(rx, path)
    return re.sub(r"\{(\w+)\}", lambda x: m.group(x.group(1)), template)


def dynamic_hits(ref, exact, globs):
    """Manifest entries matching a notebook path built with sprintf/glue placeholders.

    A placeholder first matches only hyphen-free values (callers, divergences like 0.10), so
    'calling-metrics-%s.csv' matches calling-metrics-0.10.csv but NOT calling-metrics-by-fate.csv.
    Only if that matches nothing is it widened to any value without '/' (needed for placeholders
    whose values contain hyphens, e.g. divergence-field tags like L20000-S1.0)."""
    parts = PLACEHOLDER.split(ref)
    for value_rx in (r"[^/-]+", r"[^/]+"):
        rx = re.compile("^" + value_rx.join(re.escape(p) for p in parts) + "$")
        hits = [e for f, e in exact.items() if rx.match(f)]
        # glob entries: compare the notebook pattern with the glob as text (placeholder ~ '*')
        hits += [e for e in globs if fnmatch.fnmatch(e["glob"], PLACEHOLDER.sub("*", ref))]
        if hits:
            return hits
    return []


def section_files(section):
    """Files in scope for a section (see the module docstring)."""
    files = []
    base = os.path.join(ANALYSES, section)
    for dirpath, dirnames, filenames in os.walk(base):
        dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS]
        for f in filenames:
            p = os.path.join(dirpath, f)
            if (f == "MANIFEST.yaml" or f.endswith(".dict.txt") or f.endswith(CODE_EXT)
                    or f.endswith(INDEX_EXT)):
                continue
            files.append(p)
    return sorted(files)


# --- main check ----------------------------------------------------------------------------------

def main():
    os.chdir(ROOT)
    all_sections = sorted(d for d in os.listdir(ANALYSES)
                          if os.path.isdir(os.path.join(ANALYSES, d)) and d not in NOT_SECTIONS)
    sections = sys.argv[1:] or all_sections
    errors, notices = [], []
    err = lambda code, msg: errors.append("%-16s %s" % (code, msg))
    note = lambda code, msg: notices.append("%-16s %s" % (code, msg))

    tracked = git_lines("ls-files", ANALYSES)
    ignored = git_lines("ls-files", "--others", "--ignored", "--exclude-standard", ANALYSES)

    # 1. load every manifest (all sections: notebook references can cross sections)
    manifests = {}
    for s in all_sections:
        mp = os.path.join(ANALYSES, s, "MANIFEST.yaml")
        if not os.path.exists(mp):
            if s in sections:
                err("BAD-ENTRY", "%s: no MANIFEST.yaml" % s)
            continue
        with open(mp) as fh:
            manifests[s] = yaml.safe_load(fh)
        if manifests[s].get("section") != s:
            err("BAD-ENTRY", "%s: section: is %r" % (mp, manifests[s].get("section")))

    # 2. index outputs: expanded file -> entry, plus glob entries
    exact, globs, notebooks = {}, [], {}
    for s, m in manifests.items():
        for nb in m.get("notebooks") or []:
            for k in ("path", "html", "render"):
                if not nb.get(k):
                    err("BAD-ENTRY", "%s notebook entry lacks %s" % (s, k))
            notebooks[nb.get("path")] = (s, nb)
        for e in m.get("outputs") or []:
            e["_section"] = s
            label = e.get("path") or e.get("glob")
            if not label:
                err("BAD-ENTRY", "%s: output entry with neither path nor glob" % s)
                continue
            if e.get("path"):
                try:
                    e["_files"] = expand(e["path"], e.get("values"))
                except ValueError as x:
                    err("BAD-ENTRY", "%s: %s" % (label, x))
                    continue
                for f in e["_files"]:
                    exact[f] = e
            else:
                globs.append(e)

    def entry_for(path):
        if path in exact:
            return exact[path]
        for e in globs:
            if fnmatch.fnmatch(path, e["glob"]):
                return e
        return None

    # 3. validate entries of the sections being checked
    for s in sections:
        m = manifests.get(s)
        if not m:
            continue
        for nb in m.get("notebooks") or []:
            for k in ("path", "html"):
                if nb.get(k) and not os.path.exists(nb[k]):
                    err("MISSING", "notebook %s" % nb[k])
        for e in m.get("outputs") or []:
            label = e.get("path") or e.get("glob")
            if not label or "_files" not in e and not e.get("glob"):
                continue
            for k in REQUIRED:
                if k not in e:
                    err("BAD-ENTRY", "%s: no %s" % (label, k))
            status = e.get("status", "current")
            kind = e.get("kind", default_kind(label))
            if status not in STATUSES:
                err("BAD-ENTRY", "%s: status %r" % (label, status))
            if status != "current" and not e.get("note"):
                err("BAD-ENTRY", "%s: status %s needs a note" % (label, status))
            if kind not in KINDS:
                err("BAD-ENTRY", "%s: kind %r" % (label, kind))
            if e.get("git") not in GITS:
                err("BAD-ENTRY", "%s: git %r" % (label, e.get("git")))
            if e.get("script") == "hand-written":
                if not e.get("note"):          # e.g. a caption: no script writes it
                    err("BAD-ENTRY", "%s: script hand-written needs a note (who/when/checked)" % label)
            elif e.get("script") == "none":
                # the producing code no longer exists: allowed only when the entry says so
                if status == "current" or not e.get("note"):
                    err("BAD-ENTRY", "%s: script none needs a non-current status and a note" % label)
            elif e.get("script") and not os.path.exists(e["script"]):
                err("BAD-ENTRY", "%s: script %s not found" % (label, e["script"]))
            for u in e.get("used_by") or []:
                if not os.path.exists(u):
                    err("BAD-ENTRY", "%s: used_by %s not found" % (label, u))
            # per-file checks
            files = e.get("_files") or sorted(glob.glob(e["glob"]))
            if e.get("glob") and not files and e.get("git") == "tracked":
                err("MISSING", "%s: no file matches" % label)
            for f in files:
                if not os.path.exists(f):
                    if e.get("git") == "tracked":
                        err("MISSING", f)
                    else:
                        note("ABSENT", "%s  (regenerate: %s)" % (f, fill(e.get("command", "?"), label, f)))
                    continue
                if f in tracked and e.get("git") == "ignored":
                    err("GIT", "%s is committed but the manifest says ignored" % f)
                elif f in ignored and e.get("git") == "tracked":
                    err("GIT", "%s is git-ignored but the manifest says tracked" % f)
                elif f not in tracked and f not in ignored:
                    note("UNCOMMITTED", f)
                if kind == "table":
                    d = fill(e["dictionary"], label, f) if e.get("dictionary") else f + ".dict.txt"
                    if not os.path.exists(d):
                        err("NO-DICTIONARY", "%s (expected %s)" % (f, d))

        # 4. every file on disk must be listed (or ignored with a reason)
        nb_files = set()
        for nb in m.get("notebooks") or []:
            nb_files.update([nb.get("path"), nb.get("html")])
        ign = [(i.get("glob"), i.get("reason")) for i in m.get("ignore") or []]
        for g, r in ign:
            if not g or not r:
                err("BAD-ENTRY", "%s: ignore entries need glob and reason" % s)
        for f in section_files(s):
            if f in nb_files or entry_for(f) or any(g and fnmatch.fnmatch(f, g) for g, _ in ign):
                continue
            err("UNLISTED", f)

    # 5. notebook references: every literal path a notebook reads must be listed, with used_by
    for nbpath, (s, nb) in sorted(notebooks.items()):
        if s not in sections or not nbpath or not os.path.exists(nbpath):
            continue
        with open(nbpath) as fh:
            text = fh.read()
        for ref in sorted(set(PATH_LITERAL.findall(text))):
            if (not ref.startswith(ANALYSES + "/") or ref.endswith(CODE_EXT)
                    or ref.endswith((".rmd", ".html", "/")) or os.path.isdir(ref)
                    or "." not in os.path.basename(ref)):
                continue                       # data/ inputs, sourced code, notebooks, directories
            if PLACEHOLDER.search(ref):
                hits = dynamic_hits(ref, exact, globs)
                note("DYNAMIC-REF", "%s reads %s" % (nbpath, ref))
                if not hits:
                    err("UNDOCUMENTED-REF", "%s reads %s (no listed file matches)" % (nbpath, ref))
                for e in {id(e): e for e in hits}.values():
                    if nbpath not in (e.get("used_by") or []):
                        err("USED_BY", "%s reads %s but %s's used_by lacks it"
                            % (nbpath, ref, e.get("path") or e.get("glob")))
                continue
            e = entry_for(ref)
            if e is None:
                if any(fnmatch.fnmatch(ref, g) for mm in manifests.values()
                       for g in [i.get("glob") for i in mm.get("ignore") or []] if g):
                    continue
                err("UNDOCUMENTED-REF", "%s reads %s" % (nbpath, ref))
            elif nbpath not in (e.get("used_by") or []):
                err("USED_BY", "%s reads %s but its entry's used_by lacks it" % (nbpath, ref))

    # 6. report
    for line in notices:
        print("NOTICE  " + line)
    for line in errors:
        print("ERROR   " + line)
    print("\nchecked %d section(s): %d error(s), %d notice(s)" % (len(sections), len(errors), len(notices)))
    sys.exit(1 if errors else 0)


if __name__ == "__main__":
    main()
