#!/usr/bin/env python3
"""Exercise NCBI's real traceback HSPTest deletion with a valid API option."""
from __future__ import annotations

import gzip
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from run_long_chunk_trace import (
    HERE, NCBI, NCBI_SHA256, NCBI_SOURCE, SOURCE_COMMIT, sha,
)


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: run_hsp_test_deletion_trace.py NEW_OUTPUT_DIR")
    output = Path(sys.argv[1]).resolve()
    output.mkdir(parents=True, exist_ok=False)
    assert subprocess.check_output(
        ["git", "-C", str(NCBI_SOURCE), "rev-parse", "HEAD"], text=True
    ).strip() == SOURCE_COMMIT
    assert sha(NCBI) == NCBI_SHA256
    fixture = HERE / "multi_query_20260924"
    command = [
        str(NCBI), "-task", "tblastn", "-query", str(fixture / "query.faa"),
        "-subject", str(fixture / "subjects.fna"), "-db_gencode", "1",
        "-matrix", "BLOSUM62", "-word_size", "3", "-threshold", "13",
        "-window_size", "40", "-gapopen", "11", "-gapextend", "1",
        "-evalue", "10000", "-num_threads", "1", "-comp_based_stats", "0",
        "-seg", "no", "-sum_stats", "false", "-outfmt",
        "6 qseqid sseqid score qstart qend sstart send sframe qseq sseq",
    ]
    plain = subprocess.run(command, capture_output=True, check=True)
    assert plain.stdout == (fixture / "ncbi_output.out").read_bytes()
    sources = (HERE / "ncbi_hsp_test_inject.c", HERE / "ncbi_gapped_trace.c")
    with tempfile.TemporaryDirectory(prefix="tlosan-c-hsp-delete-") as tmp:
        libraries = []
        for source in sources:
            library = Path(tmp) / (source.stem + ".so")
            subprocess.run([
                "gcc", "-shared", "-fPIC", "-std=c11", "-O2", "-Wall",
                "-Wextra", "-o", str(library), str(source), "-ldl",
            ], check=True)
            libraries.append(library)
        env = os.environ.copy()
        env["LD_PRELOAD"] = ":".join(map(str, libraries))
        traced = subprocess.run(command, capture_output=True, check=True, env=env)
    # Pinned blast_traceback.c:585-605: retain the complete probe stream
    # byte for byte, plus the ordered function/event rows used by Rust.
    (output / "trace.stderr.gz").write_bytes(gzip.compress(traced.stderr, mtime=0))
    (output / "ncbi_output.out").write_bytes(traced.stdout)
    lines = traced.stderr.decode().splitlines()
    events = [line for line in lines if line.startswith((
        "TRACEBACK_", "HSP_TEST\t", "INJECT_HSP_TEST\t", "CONTAINS\t",
        "ENDPOINT_", "TARGET_TRANSLATION\t",
    ))]
    (output / "trace.tsv").write_text("\n".join(events) + "\n")
    injected = [line for line in lines if line.startswith("INJECT_HSP_TEST\t")]
    tested = [line for line in lines if line.startswith("HSP_TEST\t")]
    assert len(injected) == len(tested) and len(tested) > 0
    positive = sum(line.startswith("HSP_TEST\t1\t") for line in tested)
    assert positive > 0
    assert traced.stdout != plain.stdout
    (output / "manifest.txt").write_text(
        "Comparison-only NCBI real traceback deletion fixture; no LOSAT runtime/build dependency.\n"
        f"NCBI source commit: {SOURCE_COMMIT}\n"
        f"NCBI binary SHA256: {NCBI_SHA256}\n"
        "Pinned source: blast_traceback.c:585-605; blast_hits.c:993-1001.\n"
        "Intervention: set BlastHitSavingOptions.percent_identity to valid 100.0 "
        "only while each actual Blast_HSPTest call executes, then restore it.\n"
        f"Original calls: {len(tested)}; positive deletions: {positive}.\n"
        f"Input SHA256: query={sha(fixture / 'query.faa')} "
        f"subject={sha(fixture / 'subjects.fna')}\n"
        f"Probe source SHA256: inject={sha(sources[0])} gapped={sha(sources[1])}\n"
        f"Baseline output SHA256: {sha(fixture / 'ncbi_output.out')}\n"
        f"Intervened output SHA256: {sha(output / 'ncbi_output.out')}\n"
        f"Command: {command!r}\n"
    )
    names = ("trace.stderr.gz", "trace.tsv", "ncbi_output.out", "manifest.txt")
    (output / "retained.sha256").write_text(
        "".join(f"{sha(output / name)}  {name}\n" for name in names)
    )
    print("HSPTest calls", len(tested), "positive deletions", positive)


if __name__ == "__main__":
    main()
