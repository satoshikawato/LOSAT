#!/usr/bin/env python3
"""Replay comparison-only NCBI local-subject Kappa heap rejection fixture."""
from __future__ import annotations
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
STAGE_D = HERE.parent
sys.path.insert(0, str(STAGE_D))
from run_ncbi_parameter_trace import NCBI, NCBI_SHA256, NCBI_SOURCE, SOURCE_COMMIT, probe_run, sha

PROBE = STAGE_D / "ncbi_kappa_traceback_trace.c"
PREFIX = b"K_TRACE_"


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: run_ncbi_heap_rejection_trace.py NEW_OUTPUT_DIR")
    out = Path(sys.argv[1]).resolve()
    out.mkdir(parents=True, exist_ok=False)
    assert subprocess.check_output(
        ["git", "-C", str(NCBI_SOURCE), "rev-parse", "HEAD"], text=True
    ).strip() == SOURCE_COMMIT
    assert sha(NCBI) == NCBI_SHA256
    command = [
        str(NCBI), "-task", "tblastn",
        "-query", str(HERE / "query.faa"),
        "-subject", str(HERE / "subjects.fna"),
        "-db_gencode", "1", "-num_threads", "1",
        "-evalue", "10", "-max_target_seqs", "2",
        "-outfmt", "6 qseqid sseqid score bitscore evalue qstart qend sstart send sframe length",
    ]
    plain = subprocess.run(command, capture_output=True, check=True)
    with tempfile.TemporaryDirectory(prefix="tlosan-d-heap-rejection-") as tmp:
        traced = probe_run(command, PROBE, Path(tmp))
    selected = [
        line for line in traced.stderr.splitlines(keepends=True)
        if line.startswith(PREFIX)
    ]
    assert traced.stdout == plain.stdout
    assert b"".join(
        line for line in traced.stderr.splitlines(keepends=True)
        if not line.startswith(PREFIX)
    ) == plain.stderr
    assert sum(line.startswith(b"K_TRACE_HEAP_WOULD\t") for line in selected) == 11
    assert sum(line.startswith(b"K_TRACE_HEAP_INSERT\t") for line in selected) == 10
    assert sum(line.startswith(b"K_TRACE_HEAP_POP\t") for line in selected) >= 11
    (out / "ncbi.out").write_bytes(plain.stdout)
    (out / "ncbi.stderr").write_bytes(plain.stderr)
    (out / "ncbi.trace").write_bytes(b"".join(selected))
    (out / "manifest.txt").write_text(
        f"NCBI source commit: {SOURCE_COMMIT}\n"
        f"NCBI binary SHA256: {NCBI_SHA256}\n"
        f"Probe SHA256: {sha(PROBE)}\n"
        f"Query SHA256: {sha(HERE / 'query.faa')}\n"
        f"Subject SHA256: {sha(HERE / 'subjects.fna')}\n"
        f"Command: {command!r}\n"
        "Comparison only; local -subject; code 1; one thread.\n"
    )
    files = sorted(path for path in out.iterdir() if path.is_file())
    (out / "outputs.sha256").write_text(
        "".join(f"{sha(path)}  {path.name}\n" for path in files)
    )


if __name__ == "__main__":
    main()
