#!/usr/bin/env python3
"""Platform checks for the LOSAT release binary (workflow platform-checks.yml).

Standard library only; it runs with the stock Python of the Linux, Windows and macOS runners.

run: for each case of platform_check_cases.json, run the binary with the performance switches on
(the default) and with LOSAT_X_DEFAULT=0 (the reference paths), at 1 and 4 threads, and require
exit code 0 and one output (one SHA-256) over all four runs.  Writes DIR/rows.tsv with the columns
case, threads, mode, exit_code, sha256, seconds.  Exit code 1 on any difference or failure.

compare: reads the rows.tsv of several platforms and requires, for every (case, threads, mode),
the same SHA-256 on all of them.  No expected hash is stored: the frozen hashes are checked on
Linux by ci.yml (ci_fast_regressions.py).  Exit code 1 when a table is missing or empty, the
tables do not hold the same keys, or any hash differs; the differing rows are printed.

usage:
  platform_checks.py run --losat BIN --out DIR [--cases FILE] [--limit N] [--timeout SECONDS]
  platform_checks.py compare --table ID=PATH [--table ID=PATH ...]
"""
import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CASES = Path(__file__).resolve().parent / "platform_check_cases.json"
THREADS = (1, 4)
MODES = ("default", "reference")
HEADER = ["case", "threads", "mode", "exit_code", "sha256", "seconds"]


def run_command(args: argparse.Namespace) -> int:
    losat = Path(args.losat).resolve()
    if not losat.exists() and losat.with_name(losat.name + ".exe").exists():
        losat = losat.with_name(losat.name + ".exe")
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    cases = json.loads(Path(args.cases).read_text())
    if args.limit is not None:
        cases = cases[: args.limit]
    # Every LOSAT_* variable of the caller is dropped; the reference mode sets the only one used.
    base = {k: v for k, v in os.environ.items() if not k.startswith("LOSAT_")}
    rows, failures = [], []
    print(f"platform {platform.platform()} machine {platform.machine()} cases {len(cases)}", flush=True)
    for case in cases:
        argv = [str(ROOT / a) if a.startswith("LOSAT/tests/") else a for a in case["argv"]]
        results = {}
        for threads in THREADS:
            for mode in MODES:
                env = dict(base)
                if mode == "reference":
                    env["LOSAT_X_DEFAULT"] = "0"
                target = out / "run.out"
                if target.exists():
                    target.unlink()
                start = time.time()
                try:
                    done = subprocess.run([str(losat), *argv, "-num_threads", str(threads), "-out", str(target)],
                                          cwd=ROOT / "LOSAT", env=env, capture_output=True, timeout=args.timeout)
                    code = str(done.returncode)
                    if done.returncode != 0:
                        sys.stdout.write(done.stderr.decode(errors="replace")[-600:] + "\n")
                except subprocess.TimeoutExpired:
                    code = "timeout"
                seconds = time.time() - start
                sha = hashlib.sha256(target.read_bytes()).hexdigest() if target.exists() else "missing"
                results[(threads, mode)] = (code, sha)
                rows.append([case["case_id"], str(threads), mode, code, sha, f"{seconds:.1f}"])
        same = len(set(results.values())) == 1 and next(iter(results.values()))[0] == "0"
        if not same:
            failures.append(case["case_id"])
        detail = " ".join(f"n{t}:{m}:{h[:10]}:rc{c}" for (t, m), (c, h) in sorted(results.items()))
        print(f"{'SAME' if same else 'DIFFER'}\t{case['case_id']}\t{detail}", flush=True)
    (out / "rows.tsv").write_text("\n".join("\t".join(r) for r in [HEADER, *rows]) + "\n")
    print(f"default = LOSAT_X_DEFAULT=0 at 1 and 4 threads, exit code 0: {len(cases) - len(failures)}/{len(cases)}; "
          f"failures: {failures}")
    return 1 if failures else 0


def read_table(path: Path) -> dict:
    lines = path.read_text().splitlines()
    if not lines or lines[0].split("\t") != HEADER:
        raise ValueError(f"{path}: not a platform-checks table")
    table = {}
    for line in lines[1:]:
        case, threads, mode, code, sha, _seconds = line.split("\t")
        table[(case, threads, mode)] = (code, sha)
    if not table:
        raise ValueError(f"{path}: no rows")
    return table


def compare_command(args: argparse.Namespace) -> int:
    tables, problems = {}, []
    for spec in args.table:
        name, _, path = spec.partition("=")
        try:
            tables[name] = read_table(Path(path))
        except (OSError, ValueError) as error:
            problems.append(f"{name}: table unusable ({error})")
    if len(tables) + len(problems) < 2:
        problems.append("fewer than two platforms given")
    names = sorted(tables)
    keys = sorted(set().union(*(t.keys() for t in tables.values()))) if tables else []
    differing = 0
    lines = []
    for key in keys:
        entries = {n: tables[n].get(key) for n in names}
        if any(e is None for e in entries.values()):
            problems.append(f"{'/'.join(key)}: missing on " + ", ".join(n for n, e in entries.items() if e is None))
        elif len({e for e in entries.values()}) != 1 or any(e[0] != "0" for e in entries.values()):
            differing += 1
            lines.append(f"DIFFER {' '.join(key)}: " + "; ".join(f"{n} rc{e[0]} {e[1][:12]}" for n, e in entries.items()))
    for line in lines + [f"PROBLEM {p}" for p in problems]:
        print(line)
    summary = (f"platforms compared: {', '.join(names) or 'none'}; keys: {len(keys)}; "
               f"differing: {differing}; problems: {len(problems)}")
    print(summary)
    step_summary = os.environ.get("GITHUB_STEP_SUMMARY")
    if step_summary:
        with open(step_summary, "a") as handle:
            handle.write("### Cross-platform hash comparison\n\n" + summary + "\n\n")
            if lines or problems:
                handle.write("```\n" + "\n".join(lines + problems) + "\n```\n")
    return 1 if differing or problems else 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run", help="default against LOSAT_X_DEFAULT=0 on this platform")
    run.add_argument("--losat", required=True)
    run.add_argument("--out", required=True)
    run.add_argument("--cases", default=str(CASES))
    run.add_argument("--limit", type=int, default=None, help="run only the first N cases (smoke runs)")
    run.add_argument("--timeout", type=float, default=1800, help="seconds per LOSAT run")
    run.set_defaults(func=run_command)
    compare = sub.add_parser("compare", help="compare the tables of several platforms")
    compare.add_argument("--table", action="append", required=True, metavar="ID=PATH")
    compare.set_defaults(func=compare_command)
    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
