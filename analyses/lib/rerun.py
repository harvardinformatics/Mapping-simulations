#!/usr/bin/env python3
"""
rerun.py -- plan the re-runs that a change requires, from the section MANIFEST.yaml files.

WHY: every analysis output is listed in its section's MANIFEST.yaml with the command that makes it,
its inputs and the notebooks that read it (check_manifest.py keeps that record honest). Together the
manifests are a dependency graph. This tool walks it, so that "what must be regenerated after I
changed X, and in what order?" is answered by the record instead of by memory.

This is the PLANNER only: it prints a plan and runs nothing.

Usage (from the project root):
  python analyses/lib/rerun.py plan --changed analyses/lib/truth_haplotypes.py [FILE ...]
      everything downstream of the changed file(s), in run order
  python analyses/lib/rerun.py plan --stale
      outputs whose producing code or inputs changed since the commit recorded in their data
      dictionary (tables) or in their provenance box (notebooks), plus everything downstream
  add --json OUT to also write the plan as JSON (for the executor, not built yet)

HOW THE GRAPH IS BUILT
  nodes  : one per manifest output entry (a file or a {placeholder} family of files) and one per
           notebook.
  edges  : producer -> consumer when a consumer entry's `inputs` names the producer's files (same
           template text, a listed file, or a placeholder match as in check_manifest.py); when a
           producer's `used_by` names a notebook or the consumer's script; and notebook -> output
           for files a notebook writes (an entry whose script is the .rmd).
  code   : each node's code dependencies are found from its script and command, then followed
           recursively: Python imports of modules in the same directory or analyses/lib, R/Rmd
           source() calls, and file names (*.py, *.R, *.sh, *.smk, *.yaml) in shell scripts and
           Snakefiles (the files they run). Dependencies are files: a change anywhere in a file
           marks everything that uses it, whether or not the changed function is called.
  stale  : (--stale) a node is stale when, since the commit recorded for it, its code or a tracked
           input changed (git diff against that commit, so uncommitted edits count), or a git-ignored
           input file is newer than it. The recorded commit is the dictionary's "Git commit", or, when
           the dictionary was later rewritten without regenerating the table (provenance.py
           table_commit), the rewrite commit: the rewrite records that the code changes up to then
           left the table unchanged. A table made with uncommitted edits has no pinned code and is
           reported as such.

ASSUMPTIONS AND LIMITS (also printed with the plan where they apply)
  - IGNORE_FOR_STALENESS: provenance.py only writes dictionary headers, so a change there does not
    make a table stale.
  - SHALLOW_IMPORTS: modules imported only for their column names / dictionary text (e.g.
    combine_summaries.py importing compare_bams.py) count as dependencies, but what they import
    does not.
  - Comments in shell scripts and Snakefiles are ignored when looking for files they run.
  - Entries with status unused / artifact are reported but not planned.
  - A step's command is the manifest's `command` up to an explanatory "(...)" note; commands with
    {placeholders} are expanded with the entry's values.
  - Snakemake re-runs a rule when its own code or declared inputs change, not when a module that
    its script imports changes: such steps are flagged as needing --forcerun.
  - Steps whose code uses Snakemake temp() outputs or `rm` (other than of a scratch $TMP directory)
    are flagged as deleting files.
  - A notebook that writes tracked files (the brief's figures) is rendered after every other
    notebook, so their provenance boxes see a clean tree.
"""
import argparse
import glob
import json
import os
import re
import subprocess
import sys
from collections import defaultdict

import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_manifest import (ANALYSES, NOT_SECTIONS, PLACEHOLDER, ROOT, expand,   # shared manifest helpers
                            fill, dynamic_hits)

# --- settings ------------------------------------------------------------------------------------
IGNORE_FOR_STALENESS = {"analyses/lib/provenance.py"}     # changes there touch dictionary headers only
# Imports used only for column names / dictionary text: the imported file is a dependency (its text
# goes into the dictionary), but its own imports are not followed (they cannot change the table).
SHALLOW_IMPORTS = {"analyses/lib/combine_summaries.py": {"compare_bams.py", "compare_vcfs.py"}}
SKIP_STATUSES = {"unused", "artifact"}                      # reported, not planned
CODE_SEARCH_DIRS = ("analyses/lib",)                        # where a bare file name is looked up
NOT_CODE = {"MANIFEST.yaml"}                                # never a code dependency
CODE_TOKEN = re.compile(r"[\w./${}-]*?([\w.-]+\.(?:py|R|r|sh|smk|yaml))\b")
PY_IMPORT = re.compile(r"^\s*(?:from\s+([\w.]+)\s+import|import\s+([\w., ]+))", re.M)
R_SOURCE = re.compile(r"""source\(\s*["']([^"']+)["']""")
DICT_COMMIT = re.compile(r"^Git commit\s*:\s*([0-9a-f]{7,40})\b(.*)$", re.M)
DICT_REWRITE = re.compile(r"^Dictionary\s*:\s*rewritten .* at commit ([0-9a-f]{7,40})\b", re.M)
HTML_COMMIT = re.compile(r"Rendered from commit\W*(?:<[^>]+>\W*)*([0-9a-f]{7,40})")


def git(*args):
    return subprocess.run(["git", "--no-optional-locks", "-C", ROOT] + list(args), check=True,
                          capture_output=True, text=True).stdout


# --- 1. manifests -> nodes ----------------------------------------------------------------------

class Node:
    def __init__(self, kind, label, section, entry):
        self.kind, self.label, self.section, self.entry = kind, label, section, entry
        self.files, self.code = [], set()
        self.parents, self.children = set(), set()
        self.changed = []                                    # changed files that make it stale

    def __repr__(self):
        return self.label


def load_nodes():
    nodes, by_file = {}, {}
    for mf in sorted(glob.glob(os.path.join(ANALYSES, "*", "MANIFEST.yaml"))):
        section = mf.split(os.sep)[1]
        if section in NOT_SECTIONS:
            continue
        m = yaml.safe_load(open(mf)) or {}
        for nb in m.get("notebooks") or []:
            n = Node("notebook", nb["path"], section, nb)
            n.files = [nb["html"]]
            nodes[n.label] = n
        for e in m.get("outputs") or []:
            label = e.get("path") or e.get("glob")
            n = Node("output", label, section, e)
            n.files = (expand(e["path"], e.get("values")) if e.get("path")
                       else sorted(glob.glob(e["glob"])))
            nodes[label] = n
            for f in n.files:
                by_file[f] = n
    return nodes, by_file


# --- 2. code dependencies --------------------------------------------------------------------------

def resolve(token, here, section):
    """A code file name mentioned in `here` -> an existing repo path, or None."""
    token = token.strip("'\"")
    if os.path.basename(token) in NOT_CODE:
        return None
    cands = [token, os.path.join(os.path.dirname(here), token)]
    if section:
        cands.append(os.path.join(ANALYSES, section, token))
    cands += [os.path.join(d, os.path.basename(token)) for d in CODE_SEARCH_DIRS]
    if section:
        cands += [os.path.join(ANALYSES, section, sub, os.path.basename(token)) for sub in ("scripts", "config")]
    for c in cands:
        c = os.path.normpath(c)
        if os.path.isfile(c) and c.startswith(ANALYSES + os.sep):
            return c
    return None


def code_text(f):
    """A file's text; for shell scripts and Snakefiles without comments (a '#' at a line start or
    after whitespace, to the end of the line)."""
    text = open(f, errors="replace").read()
    if f.endswith((".sh", ".smk")):
        text = "\n".join(re.sub(r"(^|\s)#.*$", "", line) for line in text.splitlines())
    return text


def code_closure(start, section):
    """All code files reachable from `start` (a file path), including itself. A file reached only
    through a SHALLOW_IMPORTS edge is included but not followed; reached directly anywhere, it is."""
    seen, followed, todo = set(), set(), [(start, True)]
    while todo:
        f, follow = todo.pop()
        if not f or not os.path.isfile(f):
            continue
        seen.add(f)
        if not follow or f in followed:                    # shallow import, or already followed
            continue
        followed.add(f)
        text = code_text(f)
        refs = set()
        if f.endswith(".py"):
            for a, b in PY_IMPORT.findall(text):
                for mod in ([a] if a else [x.strip().split(" as ")[0] for x in b.split(",")]):
                    refs.add(mod.split(".")[0] + ".py")
        if f.endswith((".R", ".r", ".rmd", ".Rmd")):
            refs.update(R_SOURCE.findall(text))
        if f.endswith((".sh", ".smk")):                    # the files a run script / Snakefile runs
            refs.update(m.group(1) if "/" not in m.group(0) else m.group(0).lstrip("$")
                        for m in CODE_TOKEN.finditer(text))
        for r in refs:
            p = resolve(r, f, section)
            if p and p not in followed:
                todo.append((p, os.path.basename(p) not in SHALLOW_IMPORTS.get(f, ())))
    return seen


def command_core(cmd):
    """The runnable part of a manifest command and its explanatory note: a trailing balanced
    '( ... )' group preceded by a space is a note, not part of the command."""
    cmd = " ".join((cmd or "").split())
    if cmd.endswith(")"):
        depth = 0
        for i in range(len(cmd) - 1, -1, -1):
            depth += {")": 1, "(": -1}.get(cmd[i], 0)
            if depth == 0:
                if i > 0 and cmd[i - 1] == " " and " " in cmd[i:]:
                    return cmd[:i].rstrip(), cmd[i:]
                break
    return cmd, ""


def command_files(cmd, section):
    """Code files named in a command (relative to the section dir after `cd analyses/<sec> &&`)."""
    out = set()
    for m in CODE_TOKEN.finditer(cmd):
        p = resolve(m.group(0), os.path.join(ANALYSES, section, "x"), section)
        if p:
            out.add(p)
    return out


def attach_code(nodes):
    for n in nodes.values():
        if n.kind == "notebook":
            n.code = code_closure(n.label, n.section)
            continue
        script = n.entry.get("script")
        starts = set()
        if script and os.path.isfile(script):
            starts.add(script)
        starts |= command_files(command_core(n.entry.get("command"))[0], n.section)
        for s in starts:
            n.code |= code_closure(s, n.section)


# --- 3. edges ---------------------------------------------------------------------------------------

def input_producers(inp, nodes, by_file):
    """Output nodes whose files an `inputs` item names."""
    if inp.startswith("("):                                   # explanatory note, not a path
        return []
    if inp in nodes and nodes[inp].kind == "output":         # same template text
        return [nodes[inp]]
    if inp in by_file:
        return [by_file[inp]]
    if PLACEHOLDER.search(inp):
        exact = {f: n.entry for f, n in by_file.items()}
        hits = dynamic_hits(inp, exact, [])
        return list({id(e): by_file_entry(e, nodes) for e in hits}.values())
    return []


def by_file_entry(e, nodes):
    return nodes[e.get("path") or e.get("glob")]


def link(nodes, by_file):
    def edge(a, b):
        if a is not b:
            a.children.add(b); b.parents.add(a)
    script_users = defaultdict(list)                          # script path -> entries it produces
    for n in nodes.values():
        if n.kind == "output" and n.entry.get("script"):
            script_users[n.entry["script"]].append(n)
    for n in nodes.values():
        if n.kind != "output":
            continue
        for inp in n.entry.get("inputs") or []:
            for p in input_producers(str(inp), nodes, by_file):
                edge(p, n)
        for u in n.entry.get("used_by") or []:
            if u in nodes and nodes[u].kind == "notebook":
                edge(n, nodes[u])
            for c in script_users.get(u, []):                 # used_by names another entry's script
                edge(n, c)
        if n.entry.get("script") in nodes and nodes[n.entry["script"]].kind == "notebook":
            edge(nodes[n.entry["script"]], n)                 # a file the notebook writes


# --- 4. what is stale -------------------------------------------------------------------------------

def downstream(seeds):
    out, todo = set(seeds), list(seeds)
    while todo:
        for c in todo.pop().children:
            if c not in out:
                out.add(c); todo.append(c)
    return out


def matches_input(path, n):
    """Does a changed file match one of n's inputs (template or literal)?"""
    for inp in n.entry.get("inputs") or []:
        inp = str(inp)
        if inp == path:
            return True
        if PLACEHOLDER.search(inp):
            rx = "^" + "[^/]+".join(re.escape(p) for p in PLACEHOLDER.split(inp)) + "$"
            if re.match(rx, path):
                return True
    return False


def seeds_changed(nodes, changed):
    """Nodes directly affected by changed files, with the reason."""
    why = {}
    for n in nodes.values():
        for f in changed:
            if f in n.code:
                why.setdefault(n, []).append("uses %s" % f); n.changed.append(f)
            elif n.kind == "output" and (f in n.files or matches_input(f, n)):
                why.setdefault(n, []).append("reads %s" % f)
    return why


def recorded_commit(n):
    """Commits to judge n's files against, from dictionaries / provenance boxes (see STALE above),
    and how many files have none. A commit is returned as (hash, made_with_uncommitted_edits)."""
    commits, missing = set(), 0
    if n.kind == "notebook":
        html = n.files[0]
        if os.path.exists(html):
            text = open(html, errors="replace").read()
            m = HTML_COMMIT.search(text)
            if m:
                return {(m.group(1), "UNCOMMITTED EDITS" in text.split("Rendered from commit", 1)[1][:600])}, 0
        return set(), 1
    for f in n.files:
        d = fill(n.entry["dictionary"], n.label, f) if n.entry.get("dictionary") else f + ".dict.txt"
        if not os.path.exists(d):
            missing += 1
            continue
        text = open(d, errors="replace").read()
        m, rw = DICT_COMMIT.search(text), DICT_REWRITE.search(text)
        if rw:
            commits.add((rw.group(1), False))
        elif m:
            commits.add((m.group(1), "UNCOMMITTED" in m.group(2)))
        else:
            missing += 1
    return commits, missing


def seeds_stale(nodes):
    """Nodes whose code or inputs changed since their recorded commit (see STALE in the docstring)."""
    why, unknown = {}, []
    for n in nodes.values():
        commits, missing = recorded_commit(n)
        if not commits:
            unknown.append(n)
            continue
        literal = [str(i) for i in (n.entry.get("inputs") or []) if not PLACEHOLDER.search(str(i))
                   and os.path.isfile(str(i))]
        tracked_in = [f for p in n.parents if p.entry.get("git") == "tracked" for f in p.files if os.path.exists(f)]
        ignored_in = [f for p in n.parents if p.entry.get("git") == "ignored" for f in p.files if os.path.exists(f)]
        watch = sorted((n.code - IGNORE_FOR_STALENESS) | set(literal) | set(tracked_in))
        for c, dirty in sorted(commits):
            if dirty:
                why.setdefault(n, []).append("made at %s with uncommitted edits: its code is not pinned" % c[:7])
            changed = [l for l in git("diff", "--name-only", c, "--", *watch).splitlines() if l] if watch else []
            if changed:
                why.setdefault(n, []).append("changed since %s: %s" % (c[:7], ", ".join(changed)))
                n.changed += [f for f in changed if f in n.code]
        mine = [os.path.getmtime(f) for f in n.files if os.path.exists(f)]
        newer = [f for f in ignored_in if mine and os.path.getmtime(f) > min(mine)]
        if newer:
            why.setdefault(n, []).append("newer git-ignored input: %s" % ", ".join(sorted(newer)[:3])
                                         + (" and %d more" % (len(newer) - 3) if len(newer) > 3 else ""))
    return why, unknown


# --- 5. plan ----------------------------------------------------------------------------------------

def warnings_for(n):
    w = []
    smk = sorted(f for f in n.code if f.endswith(".smk"))
    for f in smk + sorted(f for f in n.code if f.endswith(".sh")):
        text = code_text(f)
        if "temp(" in text:
            w.append("deletes files: %s marks outputs temp()" % f)
        rms = [line.strip() for line in text.splitlines()
               if re.search(r"(^|[;&|]|\bthen|\bdo)\s*rm\s", line)
               and not re.search(r"rm\s+(-\w+\s+)*\"?\$\{?(TMP|TMPDIR|tmp)\b", line)]
        if rms:
            w.append("deletes files: %s runs rm in %d place(s), e.g. %s" % (f, len(rms), rms[0][:60]))
    if smk and n.changed:
        hidden = [f for f in n.changed if not f.endswith((".smk", ".yaml"))]
        cmd = n.entry.get("command") or ""
        if hidden and "--forcerun" not in cmd:
            w.append("needs --forcerun: Snakemake does not see %s (imported, not a rule file)"
                     % ", ".join(hidden))
    return w


def stages(stale):
    """Topological levels of the stale subgraph; notebooks that write tracked files go last."""
    stale = set(stale)
    notebooks = {n for n in stale if n.kind == "notebook"}
    writers = {n for n in notebooks if any(c.kind == "output" for c in n.children)}
    deps = {n: (n.parents & stale) | ((notebooks - writers) if n in writers else set()) for n in stale}
    level, todo = {}, set(stale)
    while todo:
        ready = [n for n in todo if not (deps[n] & todo)]
        if not ready:
            raise SystemExit("dependency cycle among: %s" % ", ".join(sorted(map(str, todo))))
        k = len(set(level.values()))
        for n in ready:
            level[n] = k
        todo -= set(ready)
    out = defaultdict(list)
    for n, k in level.items():
        out[k].append(n)
    return [sorted(out[k], key=lambda n: (n.kind, n.section, n.label)) for k in sorted(out)]


def where_runs(cmds, n):
    """slurm if the command submits a job itself or through the run script it calls."""
    if any("sbatch" in c for c in cmds):
        return "slurm"
    for c in cmds:
        for f in command_files(c, n.section):
            if f.endswith(".sh") and re.search(r"^[^#]*\bsbatch\b", open(f, errors="replace").read(), re.M):
                return "slurm"
    return "local"


def steps(n):
    """Concrete commands for one node."""
    if n.kind == "notebook":
        return [" ".join((n.entry.get("render") or "").split())], ""
    core, note = command_core(n.entry.get("command"))
    names = set(re.findall(r"\{(\w+)\}", core))
    vals = {k: v for k, v in (n.entry.get("values") or {}).items() if k in names}
    if names and set(vals) == names:
        return sorted(set(expand(core, vals))), note
    return [core], note + (" [placeholders %s not in values]" % sorted(names - set(vals)) if names - set(vals) else "")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("plan", help="print what to re-run, in order (runs nothing)")
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--changed", nargs="+", metavar="FILE", help="files that changed (repo paths)")
    g.add_argument("--stale", action="store_true", help="find stale outputs from recorded commits")
    p.add_argument("--json", metavar="OUT", help="also write the plan as JSON")
    a = ap.parse_args()
    os.chdir(ROOT)

    nodes, by_file = load_nodes()
    attach_code(nodes)
    link(nodes, by_file)

    changed = [os.path.normpath(f) for f in (a.changed or [])]
    for f in changed:
        if not os.path.exists(f):
            raise SystemExit("--changed %s: no such file" % f)
    if a.changed:
        why, unknown = seeds_changed(nodes, changed), []
    else:
        why, unknown = seeds_stale(nodes)
    stale = downstream(why)
    skipped = sorted((n for n in stale if n.entry.get("status") in SKIP_STATUSES), key=str)
    plan_nodes = [n for n in stale if n.entry.get("status") not in SKIP_STATUSES]

    # --- report ---
    print("# rerun plan (%s)" % ("changed: " + ", ".join(changed) if a.changed else "stale since recorded commits"))
    print("# %d output/notebook entries to re-run; nothing has been run\n" % len(plan_nodes))
    out_json, seen_cmd = [], {}                                # command -> (first node, stage)
    for k, st in enumerate(stages(plan_nodes), 1):
        print("== stage %d" % k)
        for n in st:
            cmds, note = steps(n)
            reason = "; ".join(why[n]) if n in why else "input: " + ", ".join(
                sorted(str(x) for x in n.parents if x in stale))
            where = where_runs(cmds, n)
            key = tuple(cmds)
            if key in seen_cmd:                                # one command makes several entries
                print("   %-78s  (made by the same command as %s, stage %d)" % ((n.label,) + seen_cmd[key]))
                continue
            seen_cmd[key] = (n.label, k)
            print("   %s  [%s, %s]" % (n.label, n.section, where))
            print("      why : %s" % reason)
            for c in cmds:
                print("      run : %s" % c)
            if note:
                print("      note: %s" % note)
            warns = warnings_for(n)
            for w in warns:
                print("      WARN: %s" % w)
            out_json.append({"stage": k, "node": n.label, "kind": n.kind, "section": n.section,
                             "where": where, "commands": cmds, "why": reason, "warnings": warns})
        print()
    if skipped:
        print("== not planned (status unused/artifact): %s\n" % ", ".join(map(str, skipped)))
    tables = [n for n in unknown if n.kind == "output" and n.entry.get("kind", "table") == "table"
              and n.files and n.files[0].endswith((".csv", ".tsv"))]
    if tables and not a.changed:
        print("== tables with no recorded commit (judged only through their inputs): %d" % len(tables))
        for n in sorted(tables, key=str):
            print("   %s" % n)
    if a.json:
        json.dump(out_json, open(a.json, "w"), indent=1)


if __name__ == "__main__":
    main()
