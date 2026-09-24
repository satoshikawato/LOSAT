#!/usr/bin/env python3
"""Trace pinned NCBI mode-2 redo calls for the local 112-subject fixture."""
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
    (STAGE_D / "ncbi_kappa_mode2_trace.c",
     (b"K_CALL\t", b"K_ALIGN\t", b"K_TRANSLATED\t"), "mode2"),
    (STAGE_D / "ncbi_kappa_composition_trace.c",
     (b"K_COMP_CALL\t", b"K_MATRIX\t", b"K_COMP\t",
      b"K_COMP_RESULT\t", b"K_ADJUSTED\t"), "composition"),
)


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: run_ncbi_natural_mode2_trace.py NEW_OUTPUT_DIR")
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
        with tempfile.TemporaryDirectory(prefix="tlosan-d-natural-mode2-") as tmp:
            traced = probe_run(command, probe, Path(tmp))
        assert traced.stdout == plain.stdout
        selected = [line for line in traced.stderr.splitlines(keepends=True)
                    if line.startswith(prefix)]
        assert b"".join(line for line in traced.stderr.splitlines(keepends=True)
                        if not line.startswith(prefix)) == plain.stderr
        if label == "mode2":
            assert sum(b"\tredo_enter\t" in line for line in selected) == 11
        else:
            assert any(line.startswith(b"K_ADJUSTED\t") for line in selected)
        (output / f"{label}.tsv").write_bytes(b"".join(selected))
        manifest.append(f"{label} probe SHA256: {sha(probe)}")
    (output / "manifest.txt").write_text("\n".join(manifest) + "\n")
    files = sorted(path for path in output.iterdir() if path.is_file())
    (output / "outputs.sha256").write_text(
        "".join(f"{sha(path)}  {path.name}\n" for path in files)
    )


if __name__ == "__main__":
    main()
