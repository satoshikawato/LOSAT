#!/usr/bin/env python3
"""Trace the actual local-subject OneSubjectUpdate gate in pinned NCBI."""
from __future__ import annotations

from pathlib import Path
import subprocess
import sys
import tempfile

from run_ncbi_parameter_trace import CASES, HERE, command_from_manifest
from run_ncbi_parameter_trace import NCBI, NCBI_SHA256, NCBI_SOURCE, SOURCE_COMMIT, probe_run, sha

PROBE = HERE / "ncbi_seqsrc_callstate_trace.c"
PREFIX = b"D_SEQSRC_"
UPDATE_PREFIX = b"D_ONE_SUBJECT_"


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: run_ncbi_seqsrc_callstate.py NEW_OUTPUT_DIR")
    output = Path(sys.argv[1]).resolve()
    output.mkdir(parents=True, exist_ok=False)
    assert subprocess.check_output(
        ["git", "-C", str(NCBI_SOURCE), "rev-parse", "HEAD"], text=True
    ).strip() == SOURCE_COMMIT
    assert sha(NCBI) == NCBI_SHA256
    manifest = [
        f"NCBI source commit: {SOURCE_COMMIT}",
        f"NCBI binary SHA256: {NCBI_SHA256}",
        f"Probe SHA256: {sha(PROBE)}",
    ]
    for case in CASES:
        command = command_from_manifest(case)
        plain = subprocess.run(command, capture_output=True, check=True)
        assert plain.stdout == (HERE / "run_20260924" / f"{case}.out").read_bytes()
        with tempfile.TemporaryDirectory(prefix="tlosan-d-seqsrc-") as tmp:
            traced = probe_run(command, PROBE, Path(tmp))
        lines = traced.stderr.splitlines(keepends=True)
        selected = [line for line in lines if line.startswith((PREFIX, UPDATE_PREFIX))]
        original = b"".join(line for line in lines if not line.startswith((PREFIX, UPDATE_PREFIX)))
        assert traced.stdout == plain.stdout
        assert original == plain.stderr
        lengths = [
            int(line.split(b"\t", 1)[1])
            for line in selected
            if line.startswith(b"D_SEQSRC_TOTLEN\t")
        ]
        assert lengths and all(length > 0 for length in lengths), (case, lengths)
        assert not any(line.startswith(UPDATE_PREFIX) for line in selected), case
        (output / f"{case}.tsv").write_bytes(b"".join(selected))
        manifest.append(f"{case} command: {command!r}")
    (output / "manifest.txt").write_text("\n".join(manifest) + "\n")
    files = sorted(path for path in output.iterdir() if path.is_file())
    (output / "outputs.sha256").write_text(
        "".join(f"{sha(path)}  {path.name}\n" for path in files)
    )


if __name__ == "__main__":
    main()
