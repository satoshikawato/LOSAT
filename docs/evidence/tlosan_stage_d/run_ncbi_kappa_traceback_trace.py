#!/usr/bin/env python3
"""Save pinned comparison-only TBLASTN Kappa traceback call and input state."""
from __future__ import annotations
from pathlib import Path
import subprocess
import sys
import tempfile
from run_ncbi_parameter_trace import HERE, command_from_manifest
from run_ncbi_parameter_trace import NCBI, NCBI_SHA256, NCBI_SOURCE, SOURCE_COMMIT, probe_run, sha

PROBE = HERE / "ncbi_kappa_traceback_trace.c"
CASES = ("seg_hard_query_20260924_default", "multi_query_20260924_default", "run_20260923_default")
PREFIX = b"K_TRACE_"


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: run_ncbi_kappa_traceback_trace.py NEW_OUTPUT_DIR")
    output = Path(sys.argv[1]).resolve()
    output.mkdir(parents=True, exist_ok=False)
    assert subprocess.check_output(["git", "-C", str(NCBI_SOURCE), "rev-parse", "HEAD"], text=True).strip() == SOURCE_COMMIT
    assert sha(NCBI) == NCBI_SHA256
    manifest = [f"NCBI source commit: {SOURCE_COMMIT}", f"NCBI binary SHA256: {NCBI_SHA256}", f"Probe SHA256: {sha(PROBE)}"]
    for case in CASES:
        command = command_from_manifest(case)
        plain = subprocess.run(command, capture_output=True, check=True)
        assert plain.stdout == (HERE / "run_20260924" / f"{case}.out").read_bytes()
        with tempfile.TemporaryDirectory(prefix="tlosan-d-kappa-traceback-") as tmp:
            traced = probe_run(command, PROBE, Path(tmp))
        selected = [line for line in traced.stderr.splitlines(keepends=True) if line.startswith(PREFIX)]
        assert traced.stdout == plain.stdout
        assert b"".join(line for line in traced.stderr.splitlines(keepends=True) if not line.startswith(PREFIX)) == plain.stderr
        assert any(line.startswith(b"K_TRACE_ENTER\t") for line in selected), case
        (output / f"{case}.tsv").write_bytes(b"".join(selected))
        manifest.append(f"{case} command: {command!r}")
    (output / "manifest.txt").write_text("\n".join(manifest) + "\n")
    files = sorted(path for path in output.iterdir() if path.is_file())
    (output / "outputs.sha256").write_text("".join(f"{sha(path)}  {path.name}\n" for path in files))

if __name__ == "__main__":
    main()
