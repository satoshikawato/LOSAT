#!/usr/bin/env python3
"""Exercise pinned NCBI start-offset failure or success in real traceback."""
from __future__ import annotations

import gzip
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from run_long_chunk_trace import HERE, NCBI, NCBI_SHA256, NCBI_SOURCE, SOURCE_COMMIT, sha


def main() -> None:
    if len(sys.argv) not in (2, 3) or (len(sys.argv) == 3 and sys.argv[2] != "positive"):
        raise SystemExit("usage: run_start_failure_trace.py NEW_OUTPUT_DIR [positive]")
    positive = len(sys.argv) == 3
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
    sources = (HERE / "ncbi_start_failure_inject.c", HERE / "ncbi_gapped_trace.c")
    with tempfile.TemporaryDirectory(prefix="tlosan-c-start-failure-") as tmp:
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
        if positive:
            env["TLOSAN_START_POSITIVE"] = "1"
        else:
            env.pop("TLOSAN_START_POSITIVE", None)
        traced = subprocess.run(command, capture_output=True, check=True, env=env)
    lines = traced.stderr.decode().splitlines()
    results = [line for line in lines if line.startswith("START_RESULT\t")]
    expected_prefix = ("START_RESULT\t1\t0\t120\t200\t320\t34\t234"
                       if positive else "START_RESULT\t0\t0\t20\t0\t20\t-1\t-1")
    assert len(results) == 2 and all(line == expected_prefix for line in results)
    if positive:
        assert traced.stdout == plain.stdout
    else:
        assert traced.stdout != plain.stdout
    # Pinned blast_traceback.c:436-445: retain raw probe bytes and ordered
    # function rows so the two failures and final HSP deletion are auditable.
    (output / "trace.stderr.gz").write_bytes(gzip.compress(traced.stderr, mtime=0))
    events = [line for line in lines if line.startswith((
        "INJECT_START\t", "START_RESULT\t", "AFTER_START_HSP\t", "TRACEBACK_", "CONTAINS\t",
        "ENDPOINT_", "TARGET_TRANSLATION\t", "HSP_TEST\t",
    ))]
    (output / "trace.tsv").write_text("\n".join(events) + "\n")
    (output / "ncbi_output.out").write_bytes(traced.stdout)
    observation = (
        "Observed start acquisition: two TRUE returns (34,234), and "
        "successful gapped-start writeback after full retry.\n"
        if positive else
        "Observed start acquisition: two FALSE returns, then deletion before "
        "HSPTest, once on the partial pass and once after full retry.\n"
    )
    (output / "manifest.txt").write_text(
        f"Comparison-only NCBI real traceback start-offset {'success' if positive else 'failure'} fixture.\n"
        f"NCBI source commit: {SOURCE_COMMIT}\n"
        f"NCBI binary SHA256: {NCBI_SHA256}\n"
        "Pinned source: blast_traceback.c:436-445; blast_gapalign.c:3248-3321.\n"
        "Intervention: at query-index 0 traceback entry, replace first HSP "
        f"with query 0..{120 if positive else 20}, subject +1 "
        f"{200 if positive else 0}..{320 if positive else 20}, both gapped "
        "starts 0; preserve its pre-sorted score. This HSP input is artificial.\n"
        f"{observation}"
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
    print("start results", len(results), "positive", positive, "output bytes", len(traced.stdout))


if __name__ == "__main__":
    main()
