#!/usr/bin/env python3
"""Fresh real-subject and no-hit TBLASTN reports for codes 1, 4 and 11."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "docs/evidence/tlosan_stage_e"))
from run_stage_e_codes import calibrate_pairwise  # noqa: E402

LOSAT = ROOT / "LOSAT/target/release/LOSAT"
NCBI = Path("/home/kawato/micromamba/bin/tblastn")
API_SHA = "d8e0f100143031db5121e13907ec34c31a84e258adbf3f92241976916f1536bb"
NCBI_SHA = "e3956f1e107a30439d56c8f72fae4267a7d62ebcc16e4f99d4baf2bdf00402e0"


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def first_record(source: Path) -> bytes:
    lines = source.read_bytes().splitlines(keepends=True)
    assert lines and lines[0].startswith(b">")
    end = next((i for i, line in enumerate(lines[1:], 1) if line.startswith(b">")), len(lines))
    return b"".join(lines[:end])


def run(command: list[str]) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(command, cwd=ROOT, capture_output=True)


def first_difference(expected: bytes, actual: bytes) -> int:
    return next((i for i, (a, b) in enumerate(zip(expected, actual)) if a != b),
                min(len(expected), len(actual)))


def main() -> int:
    out, api = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
    out.mkdir(parents=True, exist_ok=False)
    assert sha(api.read_bytes()) == API_SHA and sha(NCBI.read_bytes()) == NCBI_SHA
    sources = {
        "av_code1": (1, ROOT / "LOSAT/tests/fasta/AvCLPV.faa", ROOT / "LOSAT/tests/fasta/AvCLPV.fasta"),
        "ap_code4": (4, ROOT / "LOSAT/tests/fasta/AP027131.faa", ROOT / "LOSAT/tests/fasta/AP027131.fasta"),
        "av_code11": (11, ROOT / "LOSAT/tests/fasta/AvCLPV.faa", ROOT / "LOSAT/tests/fasta/AvCLPV.fasta"),
    }
    cases = []
    for name, (code, query_source, subject) in sources.items():
        query = out / f"{name}.faa"
        query.write_bytes(first_record(query_source))
        cases.append((name, code, query, subject))
    nohit_query = out / "nohit.faa"
    nohit_subject = out / "nohit.fna"
    nohit_query.write_bytes(b">nohit_query\n" + b"W" * 60 + b"\n")
    nohit_subject.write_bytes(b">nohit_subject\n" + b"A" * 600 + b"\n")
    cases.append(("nohit_code1", 1, nohit_query, nohit_subject))
    valid_nohit_query = out / "valid_nohit.faa"
    valid_nohit_query.write_bytes(b">valid_nohit_query\n" + b"ARNDCQEGHILMFPSTWYV" * 3 + b"\n")
    cases.append(("valid_nohit_code1", 1, valid_nohit_query, nohit_subject))
    # NCBI c++/src/objects/seqfeat/gc.prt:260-266:
    # id 23; ncbieaa "FF*L..."; Base1/2/3 put TTA at the third position.
    # Replace code-1 L codons with TTA to exercise code 23's internal stop.
    code23_base = ROOT / "docs/evidence/tlosan_stage_d/all_codes_20260925/fixtures"
    code23_query = code23_base / "code23.faa"
    code23_aa = b"".join(code23_query.read_bytes().splitlines()[1:])
    code23_nt = b"".join((code23_base / "code23.fna").read_bytes().splitlines()[1:])
    assert len(code23_nt) == 3 * len(code23_aa)
    codons = [code23_nt[i:i + 3] for i in range(0, len(code23_nt), 3)]
    changed = sum(aa == ord("L") for aa in code23_aa)
    assert changed >= 3
    codons = [b"TTA" if aa == ord("L") else codon
              for aa, codon in zip(code23_aa, codons)]
    code23_subject = out / "code23_tta.fna"
    code23_subject.write_bytes(b">s_code23_tta\n" + b"".join(codons) + b"\n")
    cases.append(("code23_tta_stop", 23, code23_query, code23_subject))
    (out / "environment.json").write_text(json.dumps({
        "losat_executable_sha256": sha(LOSAT.read_bytes()),
        "ncbi_executable_sha256": NCBI_SHA, "api_executable_sha256": API_SHA,
        "ncbi_source_commit": "598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4",
        "formats": [0, 6, 7], "threads": [1, 2, 4, 8],
        "real_source_sha256": {name: {"query": sha(query_source.read_bytes()),
                                      "subject": sha(subject.read_bytes())}
                               for name, (_, query_source, subject) in sources.items()},
    }, indent=2) + "\n")
    count = 0
    with (out / "comparison.jsonl").open("w") as log:
        for name, code, query, subject in cases:
            for fmt in (0, 6, 7):
                local = ["-task", "tblastn", "-query", str(query), "-subject", str(subject),
                         "-db_gencode", str(code), "-outfmt", str(fmt)]
                control_cmd = [str(NCBI), "-task", "tblastn", "-query", str(query),
                               "-subject", str(subject), "-db_gencode", "1", "-outfmt", str(fmt)]
                control = run(control_cmd)
                assert control.returncode == 0, control.stderr.decode(errors="replace")
                if code == 1:
                    expected = control.stdout
                    selected_cmd = control_cmd
                    raw_selected_sha = sha(expected)
                    calibrated = False
                else:
                    api_control_cmd = [str(api), str(query), str(subject), "1", str(fmt)]
                    api_selected_cmd = [str(api), str(query), str(subject), str(code), str(fmt)]
                    api_control, selected = run(api_control_cmd), run(api_selected_cmd)
                    assert api_control.returncode == selected.returncode == 0
                    control_bytes = calibrate_pairwise(api_control.stdout) if fmt == 0 else api_control.stdout
                    assert control.stdout == control_bytes, (name, fmt, "CLI/API code-1 calibration")
                    expected = calibrate_pairwise(selected.stdout) if fmt == 0 else selected.stdout
                    if code == 23:
                        assert expected != control.stdout, "code-23 TTA fixture is not code sensitive"
                    selected_cmd = api_selected_cmd
                    raw_selected_sha = sha(selected.stdout)
                    calibrated = fmt == 0
                (out / f"{name}.fmt{fmt}.oracle").write_bytes(expected)
                for threads in (1, 2, 4, 8):
                    command = [str(LOSAT), "tblastn", *local, "-num_threads", str(threads)]
                    result = run(command)
                    row = {
                        "case": name, "code": code, "outfmt": fmt, "threads": threads,
                        "query": str(query), "subject": str(subject),
                        "query_sha256": sha(query.read_bytes()), "subject_sha256": sha(subject.read_bytes()),
                        "code1_cli_command": control_cmd, "code1_cli_sha256": sha(control.stdout),
                        "selected_oracle_command": selected_cmd,
                        "selected_raw_sha256": raw_selected_sha,
                        "pairwise_cover_column_calibrated": calibrated,
                        "expected_sha256": sha(expected), "actual_sha256": sha(result.stdout),
                        "code1_cli_stderr_sha256": sha(control.stderr),
                        "actual_stderr_sha256": sha(result.stderr),
                        "diagnostic_equal_to_code1_single_thread": result.stderr == control.stderr,
                        "stdout_bytes": len(result.stdout), "losat_command": command,
                        "exit": result.returncode,
                        "equal": result.returncode == 0 and result.stdout == expected
                                 and (result.stderr == control.stderr),
                    }
                    log.write(json.dumps(row, sort_keys=True) + "\n")
                    log.flush()
                    if not row["equal"]:
                        (out / "first_failure.json").write_text(json.dumps(row, indent=2) + "\n")
                        (out / "first_failure.actual").write_bytes(result.stdout)
                        (out / "first_failure.stderr").write_bytes(result.stderr)
                        (out / "first_failure.ncbi.stderr").write_bytes(control.stderr)
                        raise AssertionError(f"{name} fmt{fmt} n{threads}: byte {first_difference(expected, result.stdout)}")
                    count += 1
                print("PASS", name, f"fmt{fmt}", flush=True)
    assert count == 6 * 3 * 4
    (out / "summary.json").write_text(json.dumps({"comparisons": count, "equal": count,
                                                    "cases": list(sources) + ["nohit_code1", "valid_nohit_code1", "code23_tta_stop"]}, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
