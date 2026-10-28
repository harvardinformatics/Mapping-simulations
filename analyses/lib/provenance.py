#!/usr/bin/env python3
"""
provenance.py -- the provenance block at the top of every data dictionary in this project.

Every table a tool writes gets a <table>.dict.txt whose first lines say exactly what made it:

    Data dictionary for analyses/s04-softclip/results/softclip-toggle-0.06.csv

    Produced by : analyses/s04-softclip/run_softclip_toggle.sh
    Command     : sbatch analyses/s04-softclip/run_softclip_toggle.sh 0.06
    Git commit  : 38557cb1a2b3 (script as committed)
    Written     : 2026-10-26 10:30

Two ways to use it:
  * from a Python tool (the command line is taken from sys.argv):
        import provenance
        d.write(provenance.header(out_path))
  * from a shell run script, with the field definitions as a heredoc on stdin:
        python analyses/lib/provenance.py --table "$OUT" --script analyses/x/run_x.sh \\
            --command "sbatch analyses/x/run_x.sh $D" <<'EOF'
        ...field definitions...
        EOF

BACKFILL: a dictionary written AFTER its table already existed (older outputs that never had one)
uses backfill=True / --backfill. The block then says the commit that produced the table was not
recorded, gives the table's own modification time, and says when the dictionary was added. The
--command given for a backfill must be the command that originally produced the table.

REWRITE: a dictionary rewritten for a table whose producing commit IS known (from its previous
dictionary) uses table_commit=<commit> / --table-commit, with backfill. The block then gives that
commit and the table's modification time, and the commit the dictionary was rewritten at -- e.g.
after a fix to the dictionary text only, when re-running the table would reproduce it unchanged.

Git failures raise (no silent fallback): a dictionary must never claim a provenance it can't verify.
All paths are written relative to the project root.
"""
import argparse
import datetime
import os
import shlex
import subprocess
import sys

# project root = two levels above this file (analyses/lib/provenance.py)
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))


def relpath(path):
    """Path relative to the project root when it lies inside it, else the absolute path."""
    ap = os.path.abspath(path)
    if ap == PROJECT_ROOT or ap.startswith(PROJECT_ROOT + os.sep):
        return os.path.relpath(ap, PROJECT_ROOT)
    return ap


def _git(*args):
    """Run a git command in the project repo; raises CalledProcessError if git fails.
    --no-optional-locks: many pipeline jobs write dictionaries in parallel, and a plain
    `git status` may take the index lock to refresh it, colliding with the others."""
    return subprocess.run(["git", "--no-optional-locks", "-C", PROJECT_ROOT] + list(args),
                          check=True, capture_output=True, text=True).stdout


def git_commit():
    """Short hash of the project repo's current HEAD."""
    return _git("rev-parse", "--short=12", "HEAD").strip()


def script_state(script):
    """How the producing script relates to HEAD: committed, modified, or untracked."""
    status = _git("status", "--porcelain", "--", relpath(script)).strip()
    if status == "":
        return "script as committed"
    if status.startswith("??"):
        return "script NOT YET COMMITTED (untracked)"
    return "script has UNCOMMITTED EDITS relative to this commit"


def current_command():
    """The command line of the running Python tool, as it would be typed from the project root."""
    return "python " + relpath(sys.argv[0]) + (" " + shlex.join(sys.argv[1:]) if len(sys.argv) > 1 else "")


def _now():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M")


def header(table, script=None, command=None, backfill=False, notes=(), table_commit=None):
    """The 'Data dictionary for ...' title plus the provenance block, ending in a blank line.

    table    : path of the table being documented
    script   : producing script (default: the running Python tool, sys.argv[0])
    command  : exact producing command (default: the running tool's own command line)
    backfill : the table already existed and is NOT being regenerated -- see the module docstring
    notes    : extra 'Note' lines (e.g. that a table came from an earlier version of the script)
    table_commit : with backfill, the commit whose script produced the table (see REWRITE above)
    """
    script = script if script is not None else sys.argv[0]
    if backfill and command is None:
        raise ValueError("backfill dictionaries must be given the ORIGINAL producing command")
    if table_commit is not None and not backfill:
        raise ValueError("table_commit is only for a dictionary rewritten without regenerating the table")
    command = command if command is not None else current_command()
    lines = ["Data dictionary for %s" % relpath(table), "",
             "Produced by : %s" % relpath(script),
             "Command     : %s" % command]
    if backfill:
        if not os.path.exists(table):
            raise FileNotFoundError("backfill: table %s does not exist" % table)
        mtime = datetime.datetime.fromtimestamp(os.path.getmtime(table)).strftime("%Y-%m-%d %H:%M")
    if backfill and table_commit is not None:
        # the commit must exist in this repo (fails loudly otherwise)
        made_at = _git("rev-parse", "--short=12", "--verify", table_commit + "^{commit}").strip()
        lines += ["Git commit  : %s (table written %s by the script as committed there)" % (made_at, mtime),
                  "Dictionary  : rewritten %s at commit %s (%s; the table was not regenerated)"
                  % (_now(), git_commit(), script_state(script))]
    elif backfill:
        lines += ["Git commit  : not recorded (table written %s, before provenance tracking)" % mtime,
                  "Dictionary  : added %s at commit %s (backfill; the table was not regenerated)"
                  % (_now(), git_commit())]
    else:
        lines += ["Git commit  : %s (%s)" % (git_commit(), script_state(script)),
                  "Written     : %s" % _now()]
    for n in notes:
        lines.append("Note        : %s" % n)
    return "\n".join(lines) + "\n\n"


def write_dictionary(table, body, replace_existing=True, **kw):
    """Write <table>.dict.txt = header(table, **kw) + body. With replace_existing=False an existing
    dictionary is an error (protects hand-written dictionaries from a backfill)."""
    path = table + ".dict.txt"
    if not replace_existing and os.path.exists(path):
        raise FileExistsError("%s already exists (pass --replace to overwrite it)" % path)
    with open(path, "w") as d:
        d.write(header(table, **kw))
        d.write(body if body.endswith("\n") else body + "\n")
    return path


def main():
    ap = argparse.ArgumentParser(description="Write <table>.dict.txt: provenance block + body from stdin.")
    ap.add_argument("--table", required=True, help="table being documented")
    ap.add_argument("--script", required=True, help="script that produced the table")
    ap.add_argument("--command", required=True, help="exact command that produced the table")
    ap.add_argument("--backfill", action="store_true",
                    help="the table already existed and is not being regenerated")
    ap.add_argument("--note", action="append", default=[], help="extra Note line (repeatable)")
    ap.add_argument("--table-commit", metavar="COMMIT",
                    help="with --backfill: the commit that produced the table (a dictionary rewrite)")
    ap.add_argument("--replace", action="store_true",
                    help="with --backfill: allow overwriting an existing dictionary")
    ap.add_argument("--append-definitions", metavar="MODULE",
                    help="append the FIELD_DEFINITIONS of an analyses/lib module (e.g. compare_vcfs) "
                         "after the stdin body -- for tables assembled from that tool's rows")
    a = ap.parse_args()
    if not os.path.exists(a.table):
        sys.exit("ERROR: table %s does not exist" % a.table)
    body = sys.stdin.read()
    if not body.strip():
        sys.exit("ERROR: empty dictionary body on stdin")
    if a.append_definitions:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        module = __import__(a.append_definitions)
        body = (body.rstrip("\n") + "\n\nColumn definitions (from analyses/lib/%s.py):\n\n"
                % a.append_definitions + module.FIELD_DEFINITIONS)
    path = write_dictionary(a.table, body, replace_existing=(a.replace or not a.backfill),
                            script=a.script, command=a.command, backfill=a.backfill, notes=a.note,
                            table_commit=a.table_commit)
    print("wrote %s" % relpath(path))


if __name__ == "__main__":
    main()
