#!/usr/bin/env python3
"""Comparison-only CLI-calibrated NCBI C++ API oracle for all 27 codes."""
from __future__ import annotations

import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
LOSAT = ROOT / "LOSAT/target/release/LOSAT"
NCBI = Path("/home/kawato/micromamba/bin/tblastn")
FIXTURES = ROOT / "docs/evidence/tlosan_stage_d/all_codes_20260925/fixtures"
NCBI_SHA = "e3956f1e107a30439d56c8f72fae4267a7d62ebcc16e4f99d4baf2bdf00402e0"


def sha(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def run(command: list[str]) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(command, cwd=ROOT, capture_output=True)


def first_diff(a: bytes, b: bytes) -> int:
    return next((index for index, (x, y) in enumerate(zip(a, b)) if x != y), min(len(a), len(b)))


def calibrate_pairwise(raw: bytes) -> bytes:
    """Apply the pinned CLI's default showdefline coverage flag to API output.

    NCBI blast_format.cpp:505-524 enables the cover column only when
    m_HitsSortOption >= 0. The oracle requests -1, but the separately
    distributed formatter library still emits that optional column.
    Calibration removes only that column; code-1 raw/adjusted checks below
    guard the entire search and formatter output on every fixture.
    """
    lines = raw.splitlines(keepends=True)
    score_header = b"                                                                      Score   Query    E\n"
    field_header = b"Sequences producing significant alignments:                          (Bits)  cover   Value\n"
    in_summary = False
    calibrated = []
    for line in lines:
        if line == score_header:
            calibrated.append(b"                                                                      Score     E\n")
            continue
        if line == field_header:
            calibrated.append(b"Sequences producing significant alignments:                          (Bits)  Value\n")
            in_summary = True
            continue
        if in_summary:
            if line == b"\n":
                calibrated.append(line)
                continue
            if line.startswith(b">"):
                in_summary = False
            elif len(line) > 78:
                tail = re.sub(rb"^\d+%\s+", b"", line[78:], count=1)
                if tail != line[78:]:
                    calibrated.append(line[:78] + tail)
                    continue
        calibrated.append(line)
    return b"".join(calibrated)


def main() -> int:
    oracle, out = Path(sys.argv[1]), Path(sys.argv[2])
    assert sha(NCBI.read_bytes()) == NCBI_SHA
    assert oracle.is_file() and LOSAT.is_file()
    out.mkdir(parents=True, exist_ok=False)
    codes = [int(line.split("\t")[0]) for line in (FIXTURES.parent / "fixtures.tsv").read_text().splitlines()[1:]]
    assert len(codes) == 27 and 32 in codes
    env = {
        "ncbi_executable_sha256": NCBI_SHA,
        "ncbi_api_oracle_sha256": sha(oracle.read_bytes()),
        "losat_executable_sha256": sha(LOSAT.read_bytes()),
        "ncbi_source_commit": "598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4",
    }
    (out / "environment.json").write_text(json.dumps(env, indent=2) + "\n")
    records = []
    for code in codes:
        query = FIXTURES / f"code{code}.faa"
        subject = FIXTURES / f"code{code}.fna"
        for fmt in (0, 6, 7):
            common = ["-task", "tblastn", "-query", str(query), "-subject", str(subject),
                      "-num_threads", "1", "-outfmt", str(fmt)]
            cli_control_cmd = [str(NCBI), *common, "-db_gencode", "1"]
            api_control_cmd = [str(oracle), str(query), str(subject), "1", str(fmt)]
            cli_control = run(cli_control_cmd)
            api_control = run(api_control_cmd)
            selected_cmd = [str(oracle), str(query), str(subject), str(code), str(fmt)]
            losat_cmd = [str(LOSAT), "tblastn", *common, "-db_gencode", str(code)]
            selected = run(selected_cmd)
            losat = run(losat_cmd)
            api_control_bytes = calibrate_pairwise(api_control.stdout) if fmt == 0 else api_control.stdout
            selected_bytes = calibrate_pairwise(selected.stdout) if fmt == 0 else selected.stdout
            record = {
                "code": code, "outfmt": fmt,
                "query_sha256": sha(query.read_bytes()),
                "subject_sha256": sha(subject.read_bytes()),
                "cli_control_command": cli_control_cmd,
                "api_control_command": api_control_cmd,
                "api_selected_command": selected_cmd,
                "losat_command": losat_cmd,
                "cli_control_exit": cli_control.returncode,
                "api_control_exit": api_control.returncode,
                "api_selected_exit": selected.returncode,
                "losat_exit": losat.returncode,
                "cli_control_sha256": sha(cli_control.stdout),
                "api_control_raw_sha256": sha(api_control.stdout),
                "api_selected_raw_sha256": sha(selected.stdout),
                "api_control_sha256": sha(api_control_bytes),
                "api_selected_sha256": sha(selected_bytes),
                "losat_sha256": sha(losat.stdout),
                "calibrated": cli_control.returncode == api_control.returncode == 0
                    and cli_control.stdout == api_control_bytes,
                "equal": selected.returncode == losat.returncode == 0
                    and selected_bytes == losat.stdout,
            }
            records.append(record)
            if not record["calibrated"] or not record["equal"]:
                (out / "first_failure.json").write_text(json.dumps(record, indent=2) + "\n")
                for name, result in (("cli_control", cli_control), ("api_control", api_control),
                                     ("api_selected", selected), ("losat", losat)):
                    (out / f"first_failure.{name}.out").write_bytes(result.stdout)
                    (out / f"first_failure.{name}.err").write_bytes(result.stderr)
                (out / "first_failure.api_control.calibrated.out").write_bytes(api_control_bytes)
                (out / "first_failure.api_selected.calibrated.out").write_bytes(selected_bytes)
                print(f"FAIL code={code} fmt={fmt}: calibration byte {first_diff(cli_control.stdout, api_control_bytes)}, selected byte {first_diff(selected_bytes, losat.stdout)}", flush=True)
                break
            print(f"PASS code={code} fmt={fmt}: {len(selected.stdout)} bytes", flush=True)
        if records and (not records[-1]["calibrated"] or not records[-1]["equal"]):
            break
    (out / "comparison.jsonl").write_text("".join(json.dumps(record, sort_keys=True) + "\n" for record in records))
    print(f"{sum(record['calibrated'] and record['equal'] for record in records)}/{len(records)} byte-identical", flush=True)
    return 0 if len(records) == 81 and all(record["calibrated"] and record["equal"] for record in records) else 1


if __name__ == "__main__":
    raise SystemExit(main())
