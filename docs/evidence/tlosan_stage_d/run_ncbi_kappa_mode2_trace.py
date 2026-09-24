#!/usr/bin/env python3
"""Save pinned NCBI composition-mode-2 translated-subject Kappa call states."""
from __future__ import annotations
from pathlib import Path
import subprocess
import sys
import tempfile
from run_ncbi_parameter_trace import HERE, command_from_manifest
from run_ncbi_parameter_trace import NCBI, NCBI_SHA256, NCBI_SOURCE, SOURCE_COMMIT, probe_run, sha

PROBE = HERE / "ncbi_kappa_mode2_trace.c"
CASES = ("seg_hard_query_20260924_default", "multi_query_20260924_default", "run_20260923_default")
PREFIXES = (b"K_CALL\t", b"K_ALIGN\t", b"K_TRANSLATED\t")


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: run_ncbi_kappa_mode2_trace.py NEW_OUTPUT_DIR")
    output = Path(sys.argv[1]).resolve()
    output.mkdir(parents=True, exist_ok=False)
    assert subprocess.check_output(["git", "-C", str(NCBI_SOURCE), "rev-parse", "HEAD"], text=True).strip() == SOURCE_COMMIT
    assert sha(NCBI) == NCBI_SHA256
    manifest = [f"NCBI source commit: {SOURCE_COMMIT}", f"NCBI binary SHA256: {NCBI_SHA256}", f"Probe SHA256: {sha(PROBE)}"]
    for case in CASES:
        command = command_from_manifest(case)
        plain = subprocess.run(command, capture_output=True, check=True)
        assert plain.stdout == (HERE / "run_20260924" / f"{case}.out").read_bytes()
        with tempfile.TemporaryDirectory(prefix="tlosan-d-kappa-") as tmp:
            traced = probe_run(command, PROBE, Path(tmp))
        lines = traced.stderr.splitlines(keepends=True)
        selected = [line for line in lines if line.startswith(PREFIXES)]
        assert traced.stdout == plain.stdout
        assert b"".join(line for line in lines if not line.startswith(PREFIXES)) == plain.stderr
        assert any(b"\tredo_enter\t" in line for line in selected), case
        (output / f"{case}.tsv").write_bytes(b"".join(selected))
        manifest.append(f"{case} command: {command!r}")
    (output / "manifest.txt").write_text("\n".join(manifest) + "\n")
    files = sorted(path for path in output.iterdir() if path.is_file())
    (output / "outputs.sha256").write_text("".join(f"{sha(path)}  {path.name}\n" for path in files))

if __name__ == "__main__":
    main()
