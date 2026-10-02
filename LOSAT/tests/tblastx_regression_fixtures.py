#!/usr/bin/env python3
"""Frozen NCBI BLAST+ 2.17.0 TBLASTX outputs for fast pull-request regressions.

The cases cover TBLASTX paths that the Gate A manifest (outfmt 6 of whole genomes) and
the outfmt 0/7 fixtures (LOSAT/tests/outfmt0_manifest.tsv) do not reach: query batches
set by BATCH_SIZE, batches without a valid query, the order of the warnings and the
report (`.merged` cases), hit lists that overflow with ties, ambiguous subjects at other
thread counts, CTOOLKIT_COMPATIBLE, empty inputs and failed writes.

The format is that of blastn_regression_fixtures.py (whose helpers it uses): a case may
set environment variables (`env`) for both programs, every other variable that changes
NCBI's batches or report is unset, `oracle_env` is NCBI's configuration for an approved
exception, and a case id ending in `.merged` runs with standard error merged into
standard output (`2>&1`).

- generate: writes the inputs to LOSAT/tests/fixtures/tblastx_regression/inputs/ (the
  files are committed); the other inputs are those of the outfmt 0/7 fixtures.
- freeze --oracle TBLASTX: runs NCBI BLAST+ (comparison oracle only) from LOSAT/ and
  writes <case>.out, <case>.err and the hash columns of manifest.tsv.
- check --losat LOSAT: runs `LOSAT tblastx <argv> <losat_extra>` from LOSAT/ and compares
  stdout, stderr and the exit status with the frozen files.

Usage:
  tblastx_regression_fixtures.py generate
  tblastx_regression_fixtures.py freeze --oracle /path/to/tblastx
  tblastx_regression_fixtures.py check --losat BIN [--jobs N] [--out TSV]
"""
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from blastn_regression_fixtures import FIELDS, REPORT_ENV, run_case, sha256  # noqa: E402

ENGINE = Path(__file__).resolve().parents[1]
FIXTURES = ENGINE / "tests/fixtures/tblastx_regression"
INPUTS = FIXTURES / "inputs"
MANIFEST = FIXTURES / "manifest.tsv"

I = "tests/fixtures/tblastx_regression/inputs"
O = "tests/fasta/outfmt0"
SUBJECTS = f"-subject {O}/tblastx_multi_subject.fasta"
M = f"-query {O}/tblastx_multi_query.fasta {SUBJECTS}"
B = f"-query {O}/tblastx_batch_query.fasta {SUBJECTS}"
V = f"-query {O}/tblastx_invalid_query.fasta {SUBJECTS}"
U = f"-query {O}/tblastx_unsearched_query.fasta {SUBJECTS}"
A = f"-query {O}/tblastx_ambig_query.fasta -subject {O}/tblastx_ambig_subject.fasta"
C = f"-query {O}/tblastx_code4_query.fasta -subject {O}/tblastx_code4_subject.fasta -query_gencode 4"
Y = f"-query {O}/tblastx_many_query.fasta -subject {O}/tblastx_many_subject.fasta"
# (case_id, NCBI argv after `tblastx`, extra LOSAT-only arguments[, environment[, NCBI environment]])
_CASES = [
    # Query batches (blast_input_aux.cpp GetQueryBatchSize; the linking cutoffs use each
    # batch's average query length and smallest Lambda).
    # Queries of 1500, 20000, 800 and 6000 nt: the default batches are (1500 20000) and
    # (800 6000); 700 searches each query alone, 22000 the first three together, 100000 all.
    ("env.batch700", f"{B} -outfmt 6", "", "BATCH_SIZE=700"),
    ("env.batch22000_fmt0", f"{B} -outfmt 0", "", "BATCH_SIZE=22000"),
    ("env.batch100000_fmt7", f"{B} -outfmt 7", "", "BATCH_SIZE=100000"),
    ("env.batch700_threads4", f"{B} -outfmt 6", "-num_threads 4", "BATCH_SIZE=700"),
    ("env.batch_negative", f"{B} -outfmt 6", "", "BATCH_SIZE=-1"),
    # NCBI checks for an empty query before it reads BATCH_SIZE (tblastx_app.cpp:132-137);
    # LOSAT rejects a BATCH_SIZE that is not an integer (NCBI's CStringException), after it.
    ("env.batch_text_empty_query", f"-query {I}/empty.fa {SUBJECTS} -outfmt 6", "", "BATCH_SIZE=abc"),
    # A batch size of 0: the empty first batch fails after the outfmt 0 prolog (exit 3).
    ("env.batch0_fmt0", f"{B} -outfmt 0", "", "BATCH_SIZE=0"),
    ("env.batch0_fmt7", f"{B} -outfmt 7", "", "BATCH_SIZE=0"),
    ("env.batch0_empty_query_fmt0", f"-query {I}/empty.fa {SUBJECTS} -outfmt 0", "", "BATCH_SIZE=0"),
    # A batch of only the all-N query between searched batches (no search, warnings), then a
    # searched batch with the 2-nt query (no warning).
    ("env.batch300_invalid_fmt0.merged", f"{V} -outfmt 0", "", "BATCH_SIZE=300"),
    ("env.batch300_invalid_fmt7.merged", f"{V} -outfmt 7", "", "BATCH_SIZE=300"),
    ("env.batch300_invalid_fmt6.merged", f"{V} -outfmt 6", "", "BATCH_SIZE=300"),
    # Warnings and the report in one stream: the unsearched batch's warnings before its
    # reports, the -max_target_seqs warning before the prolog.
    ("warnings.unsearched_fmt0.merged", f"{U} -outfmt 0", ""),
    ("warnings.unsearched_fmt7.merged", f"{U} -outfmt 7", ""),
    ("warnings.unsearched_fmt6.merged", f"{U} -outfmt 6", ""),
    ("warnings.mts2_fmt0.merged", f"{C} -max_target_seqs 2 -outfmt 0", ""),
    ("warnings.mts4_fmt7.merged", f"{C} -max_target_seqs 4 -outfmt 7", ""),
    # Empty inputs: "Query is Empty!" (tblastx_app.cpp) and an empty subject (exit 3).
    ("empty.query_fmt0", f"-query {I}/empty.fa {SUBJECTS} -outfmt 0", ""),
    ("empty.query_fmt7", f"-query {I}/empty.fa {SUBJECTS} -outfmt 7", ""),
    ("empty.blank_query_fmt6", f"-query {I}/blank.fa {SUBJECTS} -outfmt 6", ""),
    ("empty.subject_fmt0", f"-query {O}/tblastx_code4_query.fasta -subject {I}/empty.fa -outfmt 0", ""),
    ("empty.subject_fmt6", f"-query {O}/tblastx_code4_query.fasta -subject {I}/empty.fa -outfmt 6", ""),
    # Hit lists of 260 subjects: overflow (e-value order, ties by subject order) and threads.
    ("hitlist.max1", f"{Y} -max_target_seqs 1 -outfmt 6", ""),
    ("hitlist.max2_fmt7", f"{Y} -max_target_seqs 2 -outfmt 7", ""),
    ("hitlist.max5", f"{Y} -max_target_seqs 5 -outfmt 6", ""),
    ("hitlist.max250", f"{Y} -max_target_seqs 250 -outfmt 6", ""),
    ("hitlist.default_threads4", f"{Y} -outfmt 6", "-num_threads 4"),
    ("hitlist.max5_fmt0_threads2", f"{Y} -max_target_seqs 5 -outfmt 0", "-num_threads 2"),
    # Ambiguous subjects (random ncbi2na bases in the preliminary search) at other thread counts.
    ("ambig.fmt6", f"{A} -outfmt 6", ""),
    ("ambig.fmt6_threads4", f"{A} -outfmt 6", "-num_threads 4"),
    ("ambig.fmt0_threads2", f"{A} -outfmt 0", "-num_threads 2"),
    # Word and window options in the other formats.
    ("options.thrwin_fmt6", f"{M} -threshold 12 -window_size 10 -outfmt 6", ""),
    ("options.thrwin_fmt7", f"{M} -threshold 16 -window_size 60 -outfmt 7", ""),
    ("options.segno_fmt7", f"{C} -seg no -outfmt 7", ""),
    # showdefline.cpp kBits is "(bits)" when CTOOLKIT_COMPATIBLE is set (also empty); the
    # tabular formats do not read it.
    ("ctoolkit.fmt0", f"{C} -outfmt 0", "", "CTOOLKIT_COMPATIBLE=1"),
    ("ctoolkit.empty_fmt0", f"{C} -max_target_seqs 3 -outfmt 0", "", "CTOOLKIT_COMPATIBLE="),
    ("ctoolkit.fmt7", f"{C} -outfmt 7", "", "CTOOLKIT_COMPATIBLE=1"),
    # A .ncbirc (found through $HOME) with entries that change no output.
    ("ncbirc.harmless_fmt0", f"{C} -outfmt 0", "",
     "HOME=tests/fixtures/blastn_regression/ncbirc_home BLAST_USAGE_REPORT=0"),
    # A failed outfmt 0 write ("BLAST failed to write output", exit 6; Linux /dev/full),
    # also with the warnings of an unsearched batch (the stream fails before the query is read).
    ("write.devfull_fmt0", f"{C} -outfmt 0 -out /dev/full", ""),
    ("write.devfull_unsearched_fmt0", f"{U} -outfmt 0 -out /dev/full", ""),
]
CASES = [(*case, *[""] * (5 - len(case))) for case in _CASES]


def command_generate(_args) -> int:
    INPUTS.mkdir(parents=True, exist_ok=True)
    (INPUTS / "empty.fa").write_bytes(b"")
    (INPUTS / "blank.fa").write_bytes(b"\n  \n\t\n")
    return 0


def read_manifest() -> list[dict[str, str]]:
    lines = [line for line in MANIFEST.read_text().splitlines() if not line.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def command_freeze(args) -> int:
    if any(key in os.environ for key in REPORT_ENV) or (Path.home() / ".ncbirc").exists():
        raise SystemExit(f"unset {REPORT_ENV} and remove ~/.ncbirc before freezing")
    oracle = str(Path(args.oracle).resolve())
    version = subprocess.run([oracle, "-version"], capture_output=True, text=True).stdout.strip().replace("\n", "; ")
    rows = []
    for case_id, argv, extra, case_env, oracle_env in CASES:
        result = run_case(oracle, [], case_id, argv, oracle_env or case_env)
        (FIXTURES / f"{case_id}.out").write_bytes(result.stdout)
        err = FIXTURES / f"{case_id}.err"
        if result.stderr:
            err.write_bytes(result.stderr)
        elif err.exists():
            err.unlink()
        rows.append({"case_id": case_id, "argv": argv, "losat_extra": extra, "exit": str(result.returncode),
                     "stdout_sha256": sha256(result.stdout), "stdout_bytes": str(len(result.stdout)),
                     "stderr_sha256": sha256(result.stderr) if result.stderr else "", "env": case_env,
                     "oracle_env": oracle_env})
        print(f"{case_id}\texit {result.returncode}\t{result.stdout.count(b'\n')} lines\t{len(result.stderr)} stderr bytes", flush=True)
    with open(MANIFEST, "w", newline="") as handle:
        handle.write(f"# Frozen TBLASTX outputs of {version} (comparison oracle only; see tblastx_regression_fixtures.py).\n")
        handle.write("# Run from LOSAT/: env <env> tblastx <argv> (LOSAT: LOSAT tblastx <argv> <losat_extra>). Written by `freeze`; do not edit.\n")
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return 0


def check_one(losat: str, row: dict[str, str]) -> tuple[str, str]:
    result = run_case(losat, ["tblastx"], row["case_id"], row["argv"], row.get("env") or "", row["losat_extra"])
    expected_out = (FIXTURES / f"{row['case_id']}.out").read_bytes()
    err_path = FIXTURES / f"{row['case_id']}.err"
    expected_err = err_path.read_bytes() if err_path.exists() else b""
    problems = []
    if sha256(expected_out) != row["stdout_sha256"]:
        problems.append("frozen stdout file does not match the manifest")
    if str(result.returncode) != row["exit"]:
        problems.append(f"exit {result.returncode}, expected {row['exit']}")
    if result.stdout != expected_out:
        expected_lines, actual_lines = expected_out.split(b"\n"), result.stdout.split(b"\n")
        line = next((n for n, (a, b) in enumerate(zip(expected_lines, actual_lines), 1) if a != b),
                    min(len(expected_lines), len(actual_lines)))
        problems.append(f"stdout differs (expected {len(expected_lines) - 1} lines, got {len(actual_lines) - 1}; first difference at line {line})")
    if result.stderr != expected_err:
        problems.append(f"stderr differs: {result.stderr[:200]!r}")
    return row["case_id"], "; ".join(problems) or "same"


def command_check(args) -> int:
    losat = str(Path(args.losat).resolve())
    rows = read_manifest()
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        results = list(pool.map(lambda row: check_one(losat, row), rows))
    lines = [f"{case_id}\t{verdict}" for case_id, verdict in results]
    if args.out:
        Path(args.out).write_text("case_id\tresult\n" + "\n".join(lines) + "\n")
    print("\n".join(lines))
    differing = [case_id for case_id, verdict in results if verdict != "same"]
    print(f"{len(results)} cases, {len(differing)} differ")
    return 1 if differing else 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="action", required=True)
    sub.add_parser("generate")
    freeze = sub.add_parser("freeze")
    freeze.add_argument("--oracle", required=True)
    check = sub.add_parser("check")
    check.add_argument("--losat", required=True)
    check.add_argument("--jobs", type=int, default=os.cpu_count() or 4)
    check.add_argument("--out")
    args = parser.parse_args()
    return {"generate": command_generate, "freeze": command_freeze, "check": command_check}[args.action](args)


if __name__ == "__main__":
    sys.exit(main())
