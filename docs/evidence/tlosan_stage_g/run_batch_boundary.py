#!/usr/bin/env python3
"""Verify NCBI TBLASTN search-skip reporting at 20,000-residue batch boundaries."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[3]
NCBI = Path("/home/kawato/micromamba/bin/tblastn")
LOSAT = ROOT / "LOSAT/target/release/LOSAT"
SUBJECT = ROOT / "docs/evidence/tlosan_stage_d/all_codes_20260925/fixtures/code1.fna"
VALID = (ROOT / "docs/evidence/tlosan_stage_d/all_codes_20260925/fixtures/code1.faa").read_bytes().split(b"\n", 1)[1].strip()


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def run(command: list[str]) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(command, cwd=ROOT, capture_output=True)


def first_diff(a: bytes, b: bytes) -> int:
    return next((i for i, (x, y) in enumerate(zip(a, b)) if x != y), min(len(a), len(b)))


def main() -> int:
    out = Path(sys.argv[1]).resolve()
    out.mkdir(parents=True, exist_ok=False)
    cases = {
        "skipped_then_valid": [("invalid_long", b"W" * 20_001), ("valid", VALID)],
        "valid_then_skipped": [("valid", VALID), ("invalid_long", b"W" * 20_001),
                               ("invalid_short", b"W" * 100)],
    }
    rows = []
    for name, entries in cases.items():
        query = out / f"{name}.faa"
        query.write_bytes(b"".join(b">" + label.encode() + b"\n" + residues + b"\n"
                                   for label, residues in entries))
        for fmt in (0, 6, 7):
            common = ["-task", "tblastn", "-query", str(query), "-subject", str(SUBJECT),
                      "-outfmt", str(fmt), "-db_gencode", "1"]
            cli_cmd = [str(NCBI), *common, "-num_threads", "1"]
            cli = run(cli_cmd)
            assert cli.returncode == 0, (name, fmt, cli.stderr.decode(errors="replace"))
            (out / f"{name}.fmt{fmt}.ncbi").write_bytes(cli.stdout)
            (out / f"{name}.fmt{fmt}.ncbi.stderr").write_bytes(cli.stderr)
            for threads in (1, 4):
                command = [str(LOSAT), "tblastn", *common, "-num_threads", str(threads)]
                actual = run(command)
                row = {"case": name, "outfmt": fmt, "threads": threads, "query_sha256": sha(query.read_bytes()),
                       "subject_sha256": sha(SUBJECT.read_bytes()), "ncbi_sha256": sha(cli.stdout),
                       "losat_sha256": sha(actual.stdout), "ncbi_stderr_sha256": sha(cli.stderr),
                       "losat_stderr_sha256": sha(actual.stderr), "ncbi_command": cli_cmd,
                       "losat_command": command, "stdout_first_diff": first_diff(cli.stdout, actual.stdout),
                       "stderr_first_diff": first_diff(cli.stderr, actual.stderr),
                       "stdout_equal": actual.stdout == cli.stdout, "stderr_equal": actual.stderr == cli.stderr,
                       "exit": actual.returncode}
                rows.append(row)
                print(name, fmt, threads, "stdout", row["stdout_equal"], "stderr", row["stderr_equal"], flush=True)
                if not row["stdout_equal"] or not row["stderr_equal"] or actual.returncode:
                    (out / "first_failure.json").write_text(json.dumps(row, indent=2) + "\n")
                    (out / "first_failure.actual").write_bytes(actual.stdout)
                    (out / "first_failure.stderr").write_bytes(actual.stderr)
                    return 1
    (out / "summary.json").write_text(json.dumps({"rows": rows, "pass": len(rows)}, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
