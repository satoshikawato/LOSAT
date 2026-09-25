#!/usr/bin/env python3
"""Compare complete real FASTA query sets with NCBI's selected-code oracle."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "docs/evidence/tlosan_stage_e"))
from run_stage_e_codes import calibrate_pairwise  # noqa: E402

LOSAT = ROOT / "LOSAT/target/release/LOSAT"
NCBI = Path("/home/kawato/micromamba/bin/tblastn")
API_SHA = "d8e0f100143031db5121e13907ec34c31a84e258adbf3f92241976916f1536bb"
BATCH_WRAPPER = Path(__file__).with_name("run_batch_api_oracle.py")
NCBI_SHA = "e3956f1e107a30439d56c8f72fae4267a7d62ebcc16e4f99d4baf2bdf00402e0"


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def count(path: Path) -> int:
    return sum(line.startswith(b">") for line in path.read_bytes().splitlines())


def run(command: list[str]) -> tuple[subprocess.CompletedProcess[bytes], float]:
    start = time.perf_counter()
    result = subprocess.run(command, cwd=ROOT, capture_output=True)
    return result, time.perf_counter() - start


def first_difference(expected: bytes, actual: bytes) -> int:
    return next((i for i, (a, b) in enumerate(zip(expected, actual)) if a != b),
                min(len(expected), len(actual)))


def main() -> int:
    out = Path(sys.argv[1]).resolve()
    api = Path(sys.argv[2]).resolve()
    out.mkdir(parents=True, exist_ok=False)
    assert sha(api.read_bytes()) == API_SHA and sha(NCBI.read_bytes()) == NCBI_SHA
    base = ROOT / "LOSAT/tests/fasta"
    cases = [
        ("av_self_code1", 1, base / "AvCLPV.faa", base / "AvCLPV.fasta", 120),
        ("av_ps_code11", 11, base / "AvCLPV.faa", base / "PsCLPV.fasta", 120),
        ("ap_self_code4", 4, base / "AP027131.faa", base / "AP027131.fasta", 595),
        ("ap_cross_code4", 4, base / "AP027131.faa", base / "AP027133.fasta", 595),
    ]
    environment = {"losat_sha256": sha(LOSAT.read_bytes()), "ncbi_sha256": NCBI_SHA,
                   "api_sha256": API_SHA, "batch_wrapper_sha256": sha(BATCH_WRAPPER.read_bytes()), "ncbi_source_commit": "598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4",
                   "cases": {name: {"queries": count(query), "query_sha256": sha(query.read_bytes()),
                                    "subject_sha256": sha(subject.read_bytes())}
                             for name, _, query, subject, _ in cases}}
    for name, _, query, subject, expected_count in cases:
        assert environment["cases"][name]["queries"] == expected_count
        assert query.is_file() and subject.is_file()
    (out / "environment.json").write_text(json.dumps(environment, indent=2) + "\n")
    rows = []
    with (out / "comparison.jsonl").open("w") as log:
        for name, code, query, subject, query_count in cases:
            for fmt in (0, 6, 7):
                print("ORACLE", name, "fmt", fmt, flush=True)
                common = ["-task", "tblastn", "-query", str(query), "-subject", str(subject),
                          "-outfmt", str(fmt)]
                cli_cmd = [str(NCBI), *common, "-db_gencode", "1"]
                cli, cli_s = run(cli_cmd)
                assert cli.returncode == 0, (name, fmt, cli.stderr.decode(errors="replace"))
                if code == 1:
                    selected = cli.stdout
                    selected_cmd = cli_cmd
                    raw_sha = sha(selected)
                    calibrated = False
                else:
                    api_control_cmd = [sys.executable, str(BATCH_WRAPPER), str(api), str(query), str(subject), "1", str(fmt)]
                    api_selected_cmd = [sys.executable, str(BATCH_WRAPPER), str(api), str(query), str(subject), str(code), str(fmt)]
                    api_control, api_control_s = run(api_control_cmd)
                    api_selected, api_selected_s = run(api_selected_cmd)
                    assert api_control.returncode == api_selected.returncode == 0, (name, fmt, api_control.stderr, api_selected.stderr)
                    control = calibrate_pairwise(api_control.stdout) if fmt == 0 else api_control.stdout
                    assert control == cli.stdout, (name, fmt, "code-1 CLI/API calibration", first_difference(cli.stdout, control))
                    selected = calibrate_pairwise(api_selected.stdout) if fmt == 0 else api_selected.stdout
                    selected_cmd = api_selected_cmd
                    raw_sha = sha(api_selected.stdout)
                    calibrated = fmt == 0
                oracle_path = out / f"{name}.fmt{fmt}.oracle"
                oracle_path.write_bytes(selected)
                for threads in (1, 2, 4, 8):
                    command = [str(LOSAT), "tblastn", *common, "-db_gencode", str(code),
                               "-num_threads", str(threads)]
                    print("LOSAT", name, "fmt", fmt, "threads", threads, flush=True)
                    result, seconds = run(command)
                    row = {"case": name, "code": code, "query_count": query_count,
                           "outfmt": fmt, "threads": threads, "command": command,
                           "cli_code1_command": cli_cmd, "cli_code1_seconds": cli_s,
                           "selected_oracle_command": selected_cmd,
                           "selected_raw_sha256": raw_sha, "pairwise_cover_column_calibrated": calibrated,
                           "query_sha256": sha(query.read_bytes()), "subject_sha256": sha(subject.read_bytes()),
                           "expected_sha256": sha(selected), "actual_sha256": sha(result.stdout),
                           "stderr_sha256": sha(result.stderr), "code1_cli_stderr_sha256": sha(cli.stderr),
                           "diagnostic_equal_to_code1_single_thread": result.stderr == cli.stderr,
                           "stdout_bytes": len(result.stdout), "elapsed_s": seconds, "exit": result.returncode,
                           "equal": result.returncode == 0 and result.stdout == selected
                                    and (result.stderr == cli.stderr)}
                    rows.append(row)
                    log.write(json.dumps(row, sort_keys=True) + "\n")
                    log.flush()
                    if not row["equal"]:
                        (out / "first_failure.json").write_text(json.dumps(row, indent=2) + "\n")
                        (out / "first_failure.actual").write_bytes(result.stdout)
                        (out / "first_failure.stderr").write_bytes(result.stderr)
                        (out / "first_failure.ncbi.stderr").write_bytes(cli.stderr)
                        raise AssertionError(f"{name} fmt{fmt} n{threads} byte {first_difference(selected, result.stdout)}")
                    print("PASS", name, fmt, threads, "seconds", round(seconds, 3), flush=True)
    (out / "summary.json").write_text(json.dumps({"comparisons": len(rows),
        "equal": sum(row["equal"] for row in rows), "query_counts": {name: n for name, _, _, _, n in cases}}, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
