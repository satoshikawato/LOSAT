#!/usr/bin/env python3
"""Trace pinned NCBI local max-target-seqs=2 preliminary calls and parameters."""
from __future__ import annotations

from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
STAGE_D = HERE.parent
sys.path.insert(0, str(STAGE_D))
from run_ncbi_parameter_trace import (  # noqa: E402
    NCBI, NCBI_SHA256, NCBI_SOURCE, SOURCE_COMMIT, probe_run, sha,
)

PROBES = (
    (STAGE_D / "ncbi_parameter_trace.c", b"D_PARAM_", "parameters"),
    (STAGE_D / "ncbi_d_call_trace.c", b"D_", "calls"),
    (HERE / "ncbi_kappa_early_termination_trace.c", b"K_EARLY_", "early"),
)


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: run_ncbi_natural_c_d_trace.py NEW_OUTPUT_DIR")
    output = Path(sys.argv[1]).resolve()
    output.mkdir(parents=True, exist_ok=False)
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
    assert plain.stdout == (HERE / "result_order_20260925" / "ncbi.out").read_bytes()
    manifest = [
        f"NCBI source commit: {SOURCE_COMMIT}",
        f"NCBI binary SHA256: {NCBI_SHA256}",
        f"Query SHA256: {sha(HERE / 'query.faa')}",
        f"Subject SHA256: {sha(HERE / 'subjects.fna')}",
        f"Command: {command!r}",
    ]
    for probe, prefix, label in PROBES:
        with tempfile.TemporaryDirectory(prefix="tlosan-d-natural-cd-") as tmp:
            traced = probe_run(command, probe, Path(tmp))
        assert traced.stdout == plain.stdout
        selected = [line for line in traced.stderr.splitlines(keepends=True)
                    if line.startswith(prefix)]
        assert selected, label
        assert b"".join(line for line in traced.stderr.splitlines(keepends=True)
                        if not line.startswith(prefix)) == plain.stderr
        (output / f"{label}.tsv").write_bytes(b"".join(selected))
        manifest.append(f"{label} probe SHA256: {sha(probe)}")
    (output / "manifest.txt").write_text("\n".join(manifest) + "\n")
    files = sorted(path for path in output.iterdir() if path.is_file())
    (output / "outputs.sha256").write_text(
        "".join(f"{sha(path)}  {path.name}\n" for path in files)
    )


if __name__ == "__main__":
    main()
