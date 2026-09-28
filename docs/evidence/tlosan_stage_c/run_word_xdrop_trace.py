#!/usr/bin/env python3
"""Pin NCBI real WordFinder behavior with a distinct query-context x-drop."""
from __future__ import annotations

import gzip
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from run_long_chunk_trace import HERE, NCBI, NCBI_SHA256, NCBI_SOURCE, SOURCE_COMMIT, sha


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: run_word_xdrop_trace.py NEW_OUTPUT_DIR")
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
    sources = (HERE / "ncbi_word_xdrop_inject.c", HERE / "ncbi_wordfinder_trace.c")
    with tempfile.TemporaryDirectory(prefix="tlosan-c-word-xdrop-") as tmp:
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
    lines = traced.stderr.decode().splitlines()
    injected = [line for line in lines if line.startswith("WORD_XDROP_INJECT\t")]
    assert injected == [f"WORD_XDROP_INJECT\t{call}\t16\t1" for call in range(6)]
    events = [line for line in lines if line.startswith((
        "WORD_XDROP_INJECT\t", "FRAME_CHUNK\t", "WORD_FINDER\t", "INIT\t",
    ))]
    original = [line for line in (fixture / "wordfinder.stderr").read_text().splitlines()
                if line.startswith("INIT\t")]
    modified = [line for line in events if line.startswith("INIT\t")]
    assert len(modified) == len(original) == 32 and modified != original
    (output / "trace.stderr.gz").write_bytes(gzip.compress(traced.stderr, mtime=0))
    (output / "trace.tsv").write_text("\n".join(events) + "\n")
    (output / "ncbi_output.out").write_bytes(traced.stdout)
    (output / "manifest.txt").write_text(
        "Comparison-only NCBI real WordFinder per-context x-drop intervention.\n"
        f"NCBI source commit: {SOURCE_COMMIT}\n"
        f"NCBI binary SHA256: {NCBI_SHA256}\n"
        "Pinned source: aa_ungapped.c:562-590.\n"
        "Artificial input: each of six actual WordFinder calls receives "
        "context-1 x-drop 1 instead of 16; restored after each call.\n"
        f"Input SHA256: query={sha(fixture / 'query.faa')} "
        f"subject={sha(fixture / 'subjects.fna')}\n"
        f"Probe source SHA256: inject={sha(sources[0])} wordfinder={sha(sources[1])}\n"
        f"Baseline output SHA256: {sha(fixture / 'ncbi_output.out')}\n"
        f"Intervened output SHA256: {sha(output / 'ncbi_output.out')}\n"
        f"Command: {command!r}\n"
    )
    names = ("trace.stderr.gz", "trace.tsv", "ncbi_output.out", "manifest.txt")
    (output / "retained.sha256").write_text(
        "".join(f"{sha(output / name)}  {name}\n" for name in names)
    )
    print("six injected calls;", len(modified), "init HSP rows; output bytes", len(traced.stdout))


if __name__ == "__main__":
    main()
