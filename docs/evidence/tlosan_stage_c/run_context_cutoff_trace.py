#!/usr/bin/env python3
"""Retain every pinned NCBI BLAST_GetGappedScore query-context cutoff."""
from __future__ import annotations

import ast
from pathlib import Path
import subprocess
import sys
import tempfile

from run_multi_query_trace import HERE, NCBI_SHA256, SOURCE_COMMIT, probe_run, sha


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: run_context_cutoff_trace.py NEW_OUTPUT_DIR")
    output = Path(sys.argv[1]).resolve()
    output.mkdir(parents=True, exist_ok=False)
    source = HERE / "ncbi_context_cutoff_trace.c"
    for profile, folder in (
        ("blosum62_word3", "multi_query_20260924"),
        ("blosum45_word2", "alternate_matrix_word2_20260924"),
    ):
        fixture = HERE / folder
        command = ast.literal_eval(
            (fixture / "manifest.txt").read_text().split("Command: ", 1)[1].splitlines()[0]
        )
        assert sha(Path(command[0])) == NCBI_SHA256
        plain = subprocess.run(command, capture_output=True, check=True)
        assert plain.stdout == (fixture / "ncbi_output.out").read_bytes()
        with tempfile.TemporaryDirectory(prefix="tlosan-context-cutoffs-") as tmp:
            traced = probe_run(command, source, Path(tmp))
        assert traced.stdout == plain.stdout
        rows = [
            line for line in traced.stderr.decode().splitlines()
            if line.startswith("GAPPED_CONTEXT_CUTOFF\t")
        ]
        assert len(rows) == 18
        (output / f"{profile}.tsv").write_text("\n".join(rows) + "\n")
    (output / "manifest.txt").write_text(
        "Comparison-only pinned NCBI BLAST_GetGappedScore context cutoff probe.\n"
        f"NCBI source commit: {SOURCE_COMMIT}\n"
        f"NCBI binary SHA256: {NCBI_SHA256}\n"
        f"Base probe source SHA256: {sha(HERE / 'ncbi_gapped_trace.c')}\n"
        f"Extension source SHA256: {sha(source)}\n"
        "Both extended-probe outputs match saved unprobed NCBI output byte for byte.\n"
        "Commands: original commands in multi_query_20260924 and alternate_matrix_word2_20260924 manifests.\n"
    )
    names = ("blosum62_word3.tsv", "blosum45_word2.tsv", "manifest.txt")
    (output / "retained.sha256").write_text(
        "".join(f"{sha(output / name)}  {name}\n" for name in names)
    )


if __name__ == "__main__":
    main()
