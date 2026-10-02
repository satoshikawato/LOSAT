#!/usr/bin/env python3
"""Fast output regressions for pull requests (no NCBI executable is run).

Runs one LOSAT release binary over the regression cases of
docs/evidence/losat_web_e1a/capture_outputs.py and checks, per case:
- output, normalized stderr, exit code and command against the S02 baseline
  (docs/evidence/losat_web_e1a/baseline/hashes.tsv);
- the frozen expected hash (Gate A, TLOSAN Stage G), with the known mismatches of
  frozen_mismatch_allowlist.json (a stale entry fails);
the outfmt 0 fixtures (docs/evidence/losat_web_e2a/check_losat.py) of the
selected programs at 1 and 4 threads, for BLASTN the frozen NCBI outputs of
blastn_regression_fixtures.py (preliminary hit lists, batches, split, ambiguity,
lowercase masking) and for TBLASTX those of tblastx_regression_fixtures.py (query
batches, warnings, hit lists, ambiguity, empty inputs, failed writes).

Programs are selected from the changed paths (--changed-from): a program's own
directory selects it and every program that imports it; any other engine or
fixture path selects all programs; documentation alone selects none. In the
default fast mode, TBLASTX runs only its short manifest cases (FAST_TBLASTX); the
nightly workflow passes --all-cases.

Usage:
  ci_fast_regressions.py --losat BIN --out DIR [--programs all|LIST | --changed-from REV] [--all-cases] [--jobs N]
  ci_fast_regressions.py --select-only --changed-from REV
"""
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
TESTS = REPO / "LOSAT" / "tests"
CAPTURE_DIR = REPO / "docs" / "evidence" / "losat_web_e1a"
FIXTURE_DIR = REPO / "docs" / "evidence" / "losat_web_e2a"
BASELINE = CAPTURE_DIR / "baseline" / "hashes.tsv"
sys.path.insert(0, str(CAPTURE_DIR))
sys.path.insert(0, str(TESTS))

import capture_outputs  # noqa: E402
import certify_platform_native_v010 as certify  # noqa: E402
from frozen_allowlist import Allowlist  # noqa: E402

PROGRAMS = capture_outputs.PROGRAMS
FIXTURE_PROGRAMS = ("blastn", "blastp", "tblastn", "tblastx")
REGRESSION_FIXTURES = {"blastn": "BLASTN", "tblastx": "TBLASTX"}
# TBLASTX manifest cases that finish in under 30 s (2026-10-01, native release binary);
# the other 13 (46 s to over 10 min) run in the nightly workflow.
FAST_TBLASTX = frozenset({
    "p02_mje_mela", "p03_mela_pemojnva", "p04_pemojnva_pesemjnv", "p07_lvmjnv_trcumjnv",
    "p08_trcumjnv_mellatmjnv", "p12_lc738874_lc738875_default", "p13_mela_mje_reverse",
})
PROGRAM_DIR = re.compile(r"^LOSAT/src/algorithm/(blastn|blastp|blastx|tblastn|tblastx)/")
STAGE_G_PREFIX = "/mnt/c/Users/genom/GitHub/LOSAT/"


def program_dependents() -> dict[str, set[str]]:
    """Programs whose code imports each program directory (transitively)."""
    imports = {}
    for program in PROGRAMS:
        text = "".join(path.read_text(errors="replace")
                       for path in (REPO / "LOSAT/src/algorithm" / program).rglob("*.rs"))
        imports[program] = set(re.findall(r"crate::algorithm::(blastn|blastp|blastx|tblastn|tblastx)\b", text)) - {program}
    dependents = {program: {program} for program in PROGRAMS}
    changed = True
    while changed:
        changed = False
        for user, used in imports.items():
            for program in PROGRAMS:
                if used & dependents[program] and user not in dependents[program]:
                    dependents[program].add(user)
                    changed = True
    return dependents


def select_programs(paths: list[str]) -> set[str]:
    dependents = program_dependents()
    selected: set[str] = set()
    for path in paths:
        match = PROGRAM_DIR.match(path)
        if match:
            selected |= dependents[match.group(1)]
        elif path.startswith(("LOSAT/", "docs/evidence/losat_web_e1a/", "docs/evidence/losat_web_e2a/",
                              "docs/evidence/tlosan_stage_g/", "docs/evidence/tlosan_stage_d/")) \
                or path == ".github/workflows/ci.yml":
            return set(PROGRAMS)
    return selected


def changed_paths(revision: str) -> list[str]:
    output = subprocess.check_output(["git", "diff", "--name-only", f"{revision}...HEAD"], cwd=REPO, text=True)
    return [line for line in output.splitlines() if line]


def stage_lexical_fixtures() -> None:
    """Copies LOSAT/tests/fasta to the Gate A lexical root that frozen outfmt 0/7 bytes name."""
    target = Path(certify.HISTORICAL_LEXICAL_ROOT) / "tests" / "fasta"
    for source in sorted((TESTS / "fasta").iterdir()):
        if not source.is_file():
            continue
        destination = target / source.name
        if destination.exists():
            if destination.read_bytes() != source.read_bytes():
                raise SystemExit(f"{destination} differs from {source}")
            continue
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)


def read_tsv(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    with open(path, newline="") as handle:
        return {(row["program"], row["case_id"]): row for row in csv.DictReader(handle, delimiter="\t")}


def run_capture(cases, losat: Path, out: Path, jobs: int) -> list[dict[str, str]]:
    timing = open(out / "timing.tsv", "w")
    timing.write("program\tcase_id\tseconds\n")

    def one(case):
        start = time.perf_counter()
        row = capture_outputs.run_case(case, losat, out)
        timing.write(f"{case.program}\t{case.case_id}\t{time.perf_counter() - start:.2f}\n")
        timing.flush()
        return row

    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as pool:
        rows = list(pool.map(one, cases))
    timing.close()
    with open(out / "hashes.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=capture_outputs.FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return rows


def check_rows(rows, baseline, allowlist: Allowlist, programs: set[str], all_cases: bool) -> tuple[list[str], list[str]]:
    failures, allowed = [], []
    for row in rows:
        key = (row["program"], row["case_id"])
        expected = baseline.get(key)
        if expected is None:
            failures.append(f"not in the S02 baseline: {key}")
        else:
            for field in ("command", "exit", "output_sha256", "stderr_sha256"):
                if expected[field] != row[field]:
                    failures.append(f"{field} differs from the S02 baseline: {key}")
        if row["expected_sha256"]:
            verdict = allowlist.classify(row["program"], row["case_id"], row["expected_sha256"], row["output_sha256"])
            if verdict == "allowed":
                allowed.append(f"{row['program']}/{row['case_id']}")
            elif verdict == "stale":
                failures.append(f"allow-list entry now matches its frozen hash (remove it): {key}")
            elif verdict == "mismatch":
                failures.append(f"frozen hash mismatch ({row['expected_source']}): {key}")
    if all_cases:
        executed = {(row["program"], row["case_id"]) for row in rows}
        for key in allowlist.unexecuted(executed, programs):
            failures.append(f"allow-list entry was not executed (remove it): {key}")
        for key in baseline:
            if key[0] in programs and key not in executed:
                failures.append(f"baseline case was not executed: {key}")
    return failures, allowed


def run_fixtures(losat: Path, programs: set[str], out: Path) -> list[str]:
    failures = []
    selected = [program for program in FIXTURE_PROGRAMS if program in programs]
    if not selected:
        return failures
    for threads in (1, 4):
        log = out / f"outfmt0-fixtures-n{threads}.tsv"
        result = subprocess.run([sys.executable, str(FIXTURE_DIR / "check_losat.py"), "--losat", str(losat),
                                 "--programs", ",".join(selected), "--threads", str(threads)],
                                cwd=FIXTURE_DIR, capture_output=True, text=True)
        log.write_text(result.stdout + result.stderr)
        if result.returncode != 0:
            failures.append(f"outfmt 0 fixtures differ at {threads} thread(s): see {log.name}")
    for program, name in REGRESSION_FIXTURES.items():
        if program not in programs:
            continue
        log = out / f"{program}-regression-fixtures.tsv"
        result = subprocess.run([sys.executable, str(TESTS / f"{program}_regression_fixtures.py"), "check",
                                 "--losat", str(losat), "--out", str(log)], capture_output=True, text=True)
        if result.returncode != 0:
            failures += [f"{name} regression fixture {line}" for line in result.stdout.splitlines()
                         if "\t" in line and not line.endswith("\tsame")]
            failures.append(f"{name} regression fixtures differ: see {log.name}")
    return failures


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--losat", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--programs", default="")
    parser.add_argument("--changed-from", default="")
    parser.add_argument("--all-cases", action="store_true")
    parser.add_argument("--select-only", action="store_true")
    parser.add_argument("--jobs", type=int, default=os.cpu_count() or 4)
    args = parser.parse_args()

    if args.changed_from:
        paths = changed_paths(args.changed_from)
        programs = select_programs(paths)
        reason = f"{len(paths)} changed paths since {args.changed_from}"
    elif args.programs in ("", "all"):
        programs, reason = set(PROGRAMS), "all programs"
    else:
        programs, reason = set(args.programs.split(",")), "given programs"
        if programs - set(PROGRAMS):
            parser.error(f"unknown programs: {sorted(programs - set(PROGRAMS))}")
    ordered = [program for program in PROGRAMS if program in programs]
    print(f"selected programs: {','.join(ordered) or 'none'} ({reason})", flush=True)
    if args.select_only:
        return 0
    if args.losat is None or args.out is None:
        parser.error("--losat and --out are required")
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    summary = {"programs": ordered, "selection": reason, "all_cases": args.all_cases, "failures": [], "allowed": []}
    if ordered:
        losat = args.losat.resolve()
        stage_lexical_fixtures()
        if "tblastn" in programs and not Path(STAGE_G_PREFIX).is_dir():
            raise SystemExit(f"TBLASTN Stage G inputs are named by {STAGE_G_PREFIX}; link the checkout there first")
        cases = capture_outputs.all_cases(programs)
        if not args.all_cases:
            cases = [case for case in cases if case.program != "tblastx" or case.case_id in FAST_TBLASTX]
        started = time.perf_counter()
        rows = run_capture(cases, losat, out, args.jobs)
        summary["capture_seconds"] = round(time.perf_counter() - started, 1)
        summary["cases"] = len(rows)
        failures, allowed = check_rows(rows, read_tsv(BASELINE), Allowlist(), programs, args.all_cases)
        started = time.perf_counter()
        failures += run_fixtures(losat, programs, out)
        summary["fixture_seconds"] = round(time.perf_counter() - started, 1)
        summary["failures"], summary["allowed"] = failures, allowed
        summary["losat_sha256"] = capture_outputs.sha256_bytes(losat.read_bytes())
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    for line in summary["allowed"]:
        print(f"allowed known mismatch: {line}")
    for line in summary["failures"]:
        print(f"FAILED: {line}")
    print(f"{summary.get('cases', 0)} cases, {len(summary['failures'])} failures, "
          f"{len(summary['allowed'])} allowed known mismatches -> {out / 'summary.json'}")
    return 1 if summary["failures"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
