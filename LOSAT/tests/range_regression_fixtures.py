#!/usr/bin/env python3
"""Frozen NCBI BLAST+ 2.17.0 outputs of `-query_loc` / `-subject_loc` searches that end with an
error or need an environment variable (Session S11, E2d), for all four programs.

The successful range searches are outfmt 0/6/7 fixtures of LOSAT/tests/outfmt0_manifest.tsv
(`e2d.*`). These cases cover what that manifest cannot hold: NCBI's range errors and their
order (exit 3 for a range NCBI cannot read, exit 1 for a subject record shorter than the
range start), query records that the range skips until a query batch is empty (`Empty
CBlastQueryVector`, exit 3, after the reports of the batches before), and a split ranged
BLASTN query (CHUNK_SIZE), whose chunk masks NCBI restricts in the searched coordinates.

The format is that of blastn_regression_fixtures.py (whose helpers it uses), with the program
as the first word of `argv`.

- freeze --bin-dir DIR: runs NCBI BLAST+ (comparison oracle only) from LOSAT/ and writes
  <case>.out, <case>.err and the hash columns of manifest.tsv.
- check --losat LOSAT: runs `LOSAT <argv> <losat_extra>` from LOSAT/ and compares stdout,
  stderr and the exit status with the frozen files.

Usage:
  range_regression_fixtures.py freeze --bin-dir /path/to/ncbi/bin
  range_regression_fixtures.py check --losat BIN [--jobs N] [--out TSV]
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
FIXTURES = ENGINE / "tests/fixtures/range_regression"
MANIFEST = FIXTURES / "manifest.tsv"

O = "tests/fasta/outfmt0"
N = f"-query {O}/e2d_n_q.fa -subject {O}/e2d_n_s_plus.fa"
P = f"-query {O}/e2d_p_q.faa -subject {O}/e2d_p_s3.faa"
T = f"-query {O}/e2d_p_q.faa -subject {O}/e2d_n_sub8k.fa"
X = f"-query {O}/e2d_x_q.fa -subject {O}/e2d_t_ts.fa"
# (case_id, NCBI argv (program first), extra LOSAT-only arguments[, environment])
_CASES = [
    # ParseSequenceRange (blast_input_aux.cpp:145-179): CBlastException, exit 3.
    ("blastn.parse_empty_range", f"blastn {N} -query_loc 5-5 -outfmt 6", ""),
    ("blastn.parse_format", f"blastn {N} -query_loc 10-20-30 -outfmt 6", ""),
    ("blastn.parse_blank", f"blastn {N} -query_loc '' -outfmt 6", ""),
    ("blastn.parse_subject_first", f"blastn {N} -query_loc 10-20-30 -subject_loc 2-1 -outfmt 6", ""),
    ("blastn.parse_subject_before_dust", f"blastn {N} -subject_loc 2-1 -dust foo -outfmt 6", ""),
    ("blastn.parse_query_after_dust", f"blastn {N} -query_loc 2-1 -dust foo -outfmt 6", ""),
    ("blastn.parse_query_before_validate", f"blastn {N} -query_loc 2-1 -evalue 0 -outfmt 6", ""),
    ("blastp.parse_zero", f"blastp {P} -subject_loc 0-10 -outfmt 6", ""),
    ("blastp.parse_hyphen_value", f"blastp {P} -query_loc=-5 -outfmt 6", ""),
    ("blastp.parse_open_end", f"blastp {P} -query_loc 10- -outfmt 0", ""),
    ("tblastn.parse_reversed", f"tblastn {T} -query_loc 20-10 -outfmt 6", ""),
    ("tblastn.parse_subject", f"tblastn {T} -subject_loc 10 -outfmt 7", ""),
    ("tblastx.parse_open_start", f"tblastx {X} -subject_loc -10 -outfmt 6", ""),
    ("tblastx.parse_query", f"tblastx {X} -query_loc 1-1 -outfmt 0", ""),
    # x_FastaToSeqLoc (blast_fasta_input.cpp:446-453): a subject record shorter than the
    # range start stops the subject reader, exit 1, before any output.
    ("blastn.subject_past_end", f"blastn -query {O}/e2d_n_qn.fa -subject {O}/e2d_n_sub_multi.fa -subject_loc 1001-2500 -outfmt 0", ""),
    ("blastp.subject_past_end", f"blastp -query {O}/e2d_p_qp75.faa -subject {O}/e2d_p_sp75m.faa -subject_loc 900-1000 -outfmt 7", ""),
    ("tblastn.subject_past_end", f"tblastn {T} -subject_loc 8003-9000 -outfmt 0", ""),
    ("tblastx.subject_past_end", f"tblastx -query {O}/e2d_n_qn.fa -subject {O}/e2d_n_sub_multi.fa -subject_loc 1001-2500 -outfmt 6", ""),
    # Every query record skipped (blast_input.cpp:144-155): the first batch is empty, NCBI
    # stops with Empty CBlastQueryVector after the outfmt 0 prolog.
    ("blastn.all_skipped", f"blastn {N} -query_loc 9000-9100 -outfmt 0", ""),
    ("blastn.all_skipped.7", f"blastn {N} -query_loc 9000-9100 -outfmt 7", ""),
    ("blastp.all_skipped", f"blastp {P} -query_loc 676-800 -outfmt 0", ""),
    ("tblastn.all_skipped", f"tblastn {T} -query_loc 676-800 -outfmt 0", ""),
    ("tblastx.all_skipped", f"tblastx {X} -query_loc 2000-3000 -outfmt 7", ""),
    ("tblastx.all_skipped.0", f"tblastx {X} -query_loc 2000-3000 -outfmt 0", ""),
    # A last batch of skipped records (whole-record batch lengths, blast_input.cpp:157-168):
    # the reports of the batches before, without the epilog, then the error.
    ("blastn.last_batch_skipped", f"blastn -query {O}/e2d_n_ab.fa -subject {O}/e2d_n_s_plus.fa -query_loc 502-1000 -outfmt 0", ""),
    ("blastn.last_batch_skipped.7", f"blastn -query {O}/e2d_n_ab.fa -subject {O}/e2d_n_s_plus.fa -query_loc 502-1000 -outfmt 7", ""),
    ("blastp.last_batch_skipped", f"blastp -query {O}/e2d_p_qabshort.faa -subject {O}/e2d_p_s3.faa -query_loc 100-300 -outfmt 0", ""),
    ("blastp.last_batch_skipped.7", f"blastp -query {O}/e2d_p_qabshort.faa -subject {O}/e2d_p_s3.faa -query_loc 100-300 -outfmt 7", ""),
    ("tblastn.last_batch_skipped", f"tblastn -query {O}/e2d_t_qabshort.faa -subject {O}/e2d_n_sub8k.fa -query_loc 100-300 -outfmt 0", ""),
    ("tblastn.last_batch_skipped.7", f"tblastn -query {O}/e2d_t_qabshort.faa -subject {O}/e2d_n_sub8k.fa -query_loc 100-300 -outfmt 7", ""),
    ("tblastx.last_batch_skipped", f"tblastx -query {O}/e2d_x_ab15.fa -subject {O}/e2d_n_s_plus.fa -query_loc 502-1000 -outfmt 0", ""),
    ("tblastx.last_batch_skipped.7", f"tblastx -query {O}/e2d_x_ab15.fa -subject {O}/e2d_n_s_plus.fa -query_loc 502-1000 -outfmt 7", ""),
    # A split ranged BLASTN query (split_query_cxx.cpp:256-281): the chunk masks are
    # restricted in the searched coordinates (NCBI's result, reproduced).
    ("blastn.split_masks", f"blastn -task blastn -query {O}/e2d_n_split_q.fa -subject {O}/e2d_n_split_s.fa -query_loc 10001-60000 -lcase_masking -outfmt 6", "", "CHUNK_SIZE=20000"),
    ("blastn.split_masks.0", f"blastn -task blastn -query {O}/e2d_n_split_q.fa -subject {O}/e2d_n_split_s.fa -query_loc 10001-60000 -lcase_masking -outfmt 0", "", "CHUNK_SIZE=20000"),
    ("blastn.split_nomask_range", f"blastn -task blastn -query {O}/e2d_n_split_q.fa -subject {O}/e2d_n_split_s.fa -query_loc 1-60000 -lcase_masking -outfmt 6", "", "CHUNK_SIZE=20000"),
]
CASES = [(case[0], case[1], case[2], case[3] if len(case) > 3 else "", "") for case in _CASES]


def read_manifest() -> list[dict[str, str]]:
    with open(MANIFEST, newline="") as handle:
        return list(csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t"))


def command_freeze(args) -> int:
    if any(key in os.environ for key in REPORT_ENV) or (Path.home() / ".ncbirc").exists():
        raise SystemExit(f"unset {REPORT_ENV} and remove ~/.ncbirc before freezing")
    bin_dir = Path(args.bin_dir).resolve()
    version = subprocess.run([str(bin_dir / "blastn"), "-version"], capture_output=True,
                             text=True).stdout.strip().replace("\n", "; ")
    FIXTURES.mkdir(parents=True, exist_ok=True)
    rows = []
    for case_id, argv, extra, case_env, oracle_env in CASES:
        program, rest = argv.split(" ", 1)
        result = run_case(str(bin_dir / program), [], case_id, rest, oracle_env or case_env)
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
        handle.write(f"# Frozen -query_loc/-subject_loc outputs of {version} (comparison oracle only; see range_regression_fixtures.py).\n")
        handle.write("# Run from LOSAT/: env <env> <argv> (LOSAT: LOSAT <argv> <losat_extra>). Written by `freeze`; do not edit.\n")
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return 0


def check_one(losat: str, row: dict[str, str]) -> tuple[str, str]:
    program, rest = row["argv"].split(" ", 1)
    result = run_case(losat, [program], row["case_id"], rest, row.get("env") or "", row["losat_extra"])
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
    freeze = sub.add_parser("freeze")
    freeze.add_argument("--bin-dir", required=True)
    check = sub.add_parser("check")
    check.add_argument("--losat", required=True)
    check.add_argument("--jobs", type=int, default=os.cpu_count() or 4)
    check.add_argument("--out")
    args = parser.parse_args()
    return {"freeze": command_freeze, "check": command_check}[args.action](args)


if __name__ == "__main__":
    sys.exit(main())
