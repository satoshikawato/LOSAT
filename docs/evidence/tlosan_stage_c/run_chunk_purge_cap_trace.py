#!/usr/bin/env python3
"""Pin positive chunk endpoint removal and frame append cap at NCBI calls."""
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
        raise SystemExit("usage: run_chunk_purge_cap_trace.py NEW_OUTPUT_DIR")
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
    sources = (HERE / "ncbi_chunk_purge_cap_inject.c", HERE / "ncbi_gapped_trace.c")
    with tempfile.TemporaryDirectory(prefix="tlosan-c-purge-cap-") as tmp:
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
    inject = [line for line in lines if line.startswith("INJECT_ENDPOINT\t")]
    assert inject == ["INJECT_ENDPOINT\t0\t20\t0\t200\t3"], inject
    append = [line for line in lines if line.startswith("APPEND_INPUT\t")]
    assert len(append) == 6 and all(line.split("\t")[2] == "3" for line in append)
    end_input = [line for line in lines if line.startswith("ENDPOINT_INPUT\t")]
    end_output = [line for line in lines if line.startswith("ENDPOINT_OUTPUT\t")]
    assert end_input[0] == "ENDPOINT_INPUT\t1\t3"
    assert end_output[0] == "ENDPOINT_OUTPUT\t2\t2"
    assert any(line.startswith("APPEND_OUTPUT\t") and line.endswith("\t3") for line in lines)
    (output / "trace.stderr.gz").write_bytes(gzip.compress(traced.stderr, mtime=0))
    events = [line for line in lines if line.startswith((
        "INJECT_", "GAPPED_", "ENDPOINT_", "MERGE_", "APPEND_",
        "TRACEBACK_", "QUERY_CONTEXT\t",
    ))]
    (output / "trace.tsv").write_text("\n".join(events) + "\n")
    (output / "ncbi_output.out").write_bytes(traced.stdout)
    (output / "manifest.txt").write_text(
        "Comparison-only NCBI positive chunk-local endpoint purge and frame append cap fixture.\n"
        f"NCBI source commit: {SOURCE_COMMIT}\n"
        f"NCBI binary SHA256: {NCBI_SHA256}\n"
        "Pinned source: blast_engine.c:539-552,840-850; blast_hits.c:2455-2537,2809-2864.\n"
        "Artificial intervention: on first 3-HSP purge input, lower-scoring "
        "context-0 HSP receives top HSP's query and subject start; first "
        "purge removes it. Each actual frame append call receives cap 3 "
        "instead of kHspNumMax.\n"
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
    print("purge 3 -> 2; frame append calls", len(append), "output bytes", len(traced.stdout))


if __name__ == "__main__":
    main()
