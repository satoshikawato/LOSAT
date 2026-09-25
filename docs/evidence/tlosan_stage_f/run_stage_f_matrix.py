#!/usr/bin/env python3
"""Recheck pinned Stage E local TBLASTN oracle bytes across native thread counts."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
E = ROOT / "docs/evidence/tlosan_stage_e"
sys.path.insert(0, str(E))
from run_stage_e_cli_matrix import LOSAT, NCBI, NCBI_SHA, PIN, SOURCE, read_commands
from run_stage_e_codes import calibrate_pairwise

THREADS = (1, 2, 4, 8)
FORMATS = (0, 6, 7)
API_SHA = "d8e0f100143031db5121e13907ec34c31a84e258adbf3f92241976916f1536bb"
REPEAT_NAMES = ("multi_query", "multi_hsp", "tie", "seg_hard", "seg_soft",
                "natural_c_d_early", "natural_mode2", "result_order",
                "containment_run", "heap", "run_20260923")


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def records(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text().splitlines()]


def run(command: list[str], debug: bool = False) -> subprocess.CompletedProcess[bytes]:
    env = dict(os.environ)
    if debug:
        env["LOSAT_WASI_THREADS_DEBUG"] = "1"
    return subprocess.run(command, cwd=ROOT, capture_output=True, env=env)


def verify_environment(api: Path | None) -> dict:
    assert subprocess.check_output(["git", "-C", str(SOURCE), "rev-parse", "HEAD"]).decode().strip() == PIN
    assert sha(NCBI.read_bytes()) == NCBI_SHA
    assert LOSAT.is_file()
    if api is not None:
        assert api.is_file() and sha(api.read_bytes()) == API_SHA
    return {
        "ncbi_source_commit": PIN, "ncbi_executable": str(NCBI),
        "ncbi_executable_sha256": NCBI_SHA, "losat_executable": str(LOSAT),
        "losat_executable_sha256": sha(LOSAT.read_bytes()),
        "ncbi_api_oracle": str(api) if api else None,
        "ncbi_api_oracle_sha256": API_SHA if api else None,
        "threads": THREADS, "formats": FORMATS,
    }


def first_difference(expected: bytes, actual: bytes) -> int:
    return next((i for i, (a, b) in enumerate(zip(expected, actual)) if a != b),
                min(len(expected), len(actual)))


# NCBI c++/src/algo/blast/core/blast_hspstream.c:289-319:
# *hsp_list_out = hit_list->hsplist_array[last_hsplist_index];
# --hit_list->hsplist_count;
# Compare complete output bytes after indexed subject work enters that stream.
def compare_one(out: Path, log, name: str, query: Path, subject: Path,
                common: list[str], fmt: int, expected: bytes,
                authority: dict, registered: dict, repeats: bool) -> None:
    query_sha, subject_sha = sha(query.read_bytes()), sha(subject.read_bytes())
    assert query_sha == registered["query_sha256"], (name, "query")
    assert subject_sha == registered["subject_sha256"], (name, "subject")
    assert sha(expected) == authority["registered_sha256"], (name, "NCBI oracle")
    for n in THREADS:
        count = 3 if repeats and n > 1 else 1
        for repetition in range(count):
            args = common.copy()
            if "-num_threads" in args:
                args[args.index("-num_threads") + 1] = str(n)
            else:
                args.extend(["-num_threads", str(n)])
            command = [str(LOSAT), "tblastn", *args, "-outfmt", str(fmt)]
            actual = run(command, debug=True)
            activity = re.search(rb"\[losat-thread-activity\].*peak_active=(\d+)", actual.stderr)
            row = {
                "case": name, "outfmt": fmt, "target": "x86_64-unknown-linux-gnu",
                "threads": n, "repetition": repetition, "query": str(query),
                "subject": str(subject), "query_sha256": query_sha,
                "subject_sha256": subject_sha, "losat_command": command,
                "losat_executable_sha256": sha(LOSAT.read_bytes()),
                "ncbi_contract": authority, "expected_sha256": sha(expected),
                "actual_sha256": sha(actual.stdout),
                "stdout_bytes": len(actual.stdout), "stderr_sha256": sha(actual.stderr),
                "exit": actual.returncode, "peak_active": int(activity.group(1)) if activity else None,
                "activity_log": activity.group(0).decode() if activity else None,
                "equal": actual.returncode == 0 and actual.stdout == expected,
            }
            log.write(json.dumps(row, sort_keys=True) + "\n")
            log.flush()
            if not row["equal"]:
                (out / "first_failure.json").write_text(json.dumps(row, indent=2) + "\n")
                (out / "first_failure.expected").write_bytes(expected)
                (out / "first_failure.actual").write_bytes(actual.stdout)
                (out / "first_failure.stderr").write_bytes(actual.stderr)
                raise AssertionError(f"{name} fmt{fmt} n{n} repeat{repetition} first byte {first_difference(expected, actual.stdout)}")


def physical(out: Path, log) -> None:
    baseline = {row["case"]: row for row in records(E / "cli_matrix_20260925/comparison.jsonl")}
    cases = list(read_commands())
    assert len(cases) == 73 and len(baseline) == 219
    for family, case, base in cases:
        query = Path(base[base.index("-query") + 1])
        subject = Path(base[base.index("-subject") + 1])
        stem = f"{family.replace('/', '__')}__{case}"
        for fmt in FORMATS:
            name = f"{stem}__fmt{fmt}"
            saved = baseline[name]
            assert saved["equal"]
            oracle_cmd = [str(NCBI), *base[1:], "-outfmt", str(fmt)]
            oracle = run(oracle_cmd)
            assert oracle.returncode == 0 and sha(oracle.stdout) == saved["ncbi_stdout_sha256"], name
            authority = {
                "kind": "pinned_ncbi_cli_local_subject", "command": oracle_cmd,
                "executable_sha256": NCBI_SHA,
                "registered_sha256": saved["ncbi_stdout_sha256"],
                "current_stderr_sha256": sha(oracle.stderr),
            }
            compare_one(out, log, name, query, subject, base[1:], fmt, oracle.stdout,
                        authority, saved, any(key in name for key in REPEAT_NAMES))
            print("PASS", name, flush=True)


# NCBI c++/src/algo/blast/core/blast_engine.c:1452-1465:
# seq_arg.seq->gen_code_string = GenCodeSingletonFind(db_options->genetic_code);
# stat_length /= CODON_LENGTH;
# The API oracle selects FindGeneticCode(code); code 1 is calibrated to the CLI.
def codes(out: Path, log, api: Path) -> None:
    baseline = {(row["code"], row["outfmt"]): row for row in records(E / "all_codes_20260925/comparison.jsonl")}
    assert len(baseline) == 81
    fixtures = ROOT / "docs/evidence/tlosan_stage_d/all_codes_20260925/fixtures"
    for code in sorted({key[0] for key in baseline}):
        query, subject = fixtures / f"code{code}.faa", fixtures / f"code{code}.fna"
        for fmt in FORMATS:
            saved = baseline[(code, fmt)]
            assert saved["equal"] and saved["calibrated"]
            control_cli = run([str(NCBI), "-task", "tblastn", "-query", str(query),
                               "-subject", str(subject), "-num_threads", "1",
                               "-outfmt", str(fmt), "-db_gencode", "1"])
            control_api = run([str(api), str(query), str(subject), "1", str(fmt)])
            selected_cmd = [str(api), str(query), str(subject), str(code), str(fmt)]
            selected = run(selected_cmd)
            control_bytes = calibrate_pairwise(control_api.stdout) if fmt == 0 else control_api.stdout
            expected = calibrate_pairwise(selected.stdout) if fmt == 0 else selected.stdout
            assert control_cli.returncode == control_api.returncode == selected.returncode == 0
            assert control_cli.stdout == control_bytes
            assert sha(control_cli.stdout) == saved["cli_control_sha256"]
            assert sha(expected) == saved["api_selected_sha256"]
            authority = {
                "kind": "pinned_ncbi_cpp_api_FindGeneticCode_local_subject",
                "selected_command": selected_cmd,
                "code1_cli_command": [str(NCBI), "-task", "tblastn", "-query", str(query),
                    "-subject", str(subject), "-num_threads", "1", "-outfmt", str(fmt),
                    "-db_gencode", "1"],
                "code1_cli_sha256": sha(control_cli.stdout),
                "code1_api_calibrated_sha256": sha(control_bytes),
                "api_executable_sha256": API_SHA,
                "registered_sha256": saved["api_selected_sha256"],
                "calibration": "Stage E pairwise summary coverage column only" if fmt == 0 else "none",
            }
            base = ["-task", "tblastn", "-query", str(query), "-subject", str(subject),
                    "-db_gencode", str(code), "-num_threads", "1"]
            compare_one(out, log, f"code{code}", query, subject, base, fmt, expected,
                        authority, saved, code in (1, 4, 11, 32))
            print("PASS", f"code{code} fmt{fmt}", flush=True)


def codes_multi(out: Path, log, api: Path) -> None:
    # NCBI c++/src/algo/blast/core/blast_engine.c:1410-1416:
    # each FASTA subject supplies an independent seq_arg.oid to preliminary search.
    fixtures = ROOT / "docs/evidence/tlosan_stage_d/all_codes_20260925/fixtures"
    inputs = out / "generated_inputs"
    inputs.mkdir()
    for code in sorted(int(line.split("\t")[0]) for line in
                       (fixtures.parent / "fixtures.tsv").read_text().splitlines()[1:]):
        query = fixtures / f"code{code}.faa"
        original = (fixtures / f"code{code}.fna").read_text().splitlines()
        assert len(original) == 2 and original[0].startswith(">")
        subject = inputs / f"code{code}_eight_subjects.fna"
        subject.write_text("".join(f">s_code{code}_copy{i}\n{original[1]}\n" for i in range(8)))
        for fmt in FORMATS:
            control_cli_cmd = [str(NCBI), "-task", "tblastn", "-query", str(query),
                "-subject", str(subject), "-num_threads", "1", "-outfmt", str(fmt),
                "-db_gencode", "1"]
            control_cli = run(control_cli_cmd)
            control_api = run([str(api), str(query), str(subject), "1", str(fmt)])
            selected_cmd = [str(api), str(query), str(subject), str(code), str(fmt)]
            selected = run(selected_cmd)
            control_bytes = calibrate_pairwise(control_api.stdout) if fmt == 0 else control_api.stdout
            expected = calibrate_pairwise(selected.stdout) if fmt == 0 else selected.stdout
            assert control_cli.returncode == control_api.returncode == selected.returncode == 0, (code, fmt)
            assert control_cli.stdout == control_bytes, (code, fmt, "CLI/API calibration")
            authority = {
                "kind": "pinned_ncbi_cpp_api_FindGeneticCode_local_subject",
                "selected_command": selected_cmd,
                "code1_cli_command": control_cli_cmd,
                "code1_cli_sha256": sha(control_cli.stdout),
                "code1_api_calibrated_sha256": sha(control_bytes),
                "api_executable_sha256": API_SHA,
                "registered_sha256": sha(expected),
                "calibration": "Stage E pairwise summary coverage column only" if fmt == 0 else "none",
            }
            base = ["-task", "tblastn", "-query", str(query), "-subject", str(subject),
                    "-db_gencode", str(code), "-num_threads", "1"]
            registered = {"query_sha256": sha(query.read_bytes()),
                          "subject_sha256": sha(subject.read_bytes())}
            compare_one(out, log, f"code{code}_eight_subjects", query, subject, base,
                        fmt, expected, authority, registered, True)
            print("PASS", f"code{code} eight subjects fmt{fmt}", flush=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("family", choices=("physical", "codes", "codes_multi"))
    parser.add_argument("out", type=Path)
    parser.add_argument("--api-oracle", type=Path)
    args = parser.parse_args()
    if args.family in ("codes", "codes_multi") and args.api_oracle is None:
        parser.error("codes requires --api-oracle")
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    (out / "environment.json").write_text(json.dumps(verify_environment(args.api_oracle), indent=2) + "\n")
    with (out / "comparison.jsonl").open("w") as log:
        if args.family == "physical":
            physical(out, log)
        elif args.family == "codes":
            codes(out, log, args.api_oracle.resolve())
        else:
            codes_multi(out, log, args.api_oracle.resolve())
    rows = records(out / "comparison.jsonl")
    (out / "summary.json").write_text(json.dumps({
        "family": args.family, "comparisons": len(rows), "equal": sum(row["equal"] for row in rows),
        "unique_case_formats": len({(row["case"], row["outfmt"]) for row in rows}),
        "by_thread": {str(n): sum(row["threads"] == n for row in rows) for n in THREADS},
        "parallel_activity_max": {str(n): max((row["peak_active"] or 0 for row in rows if row["threads"] == n), default=0)
                                  for n in THREADS},
    }, indent=2) + "\n")
    print(f"{len(rows)}/{len(rows)} native thread outputs match registered NCBI bytes", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
