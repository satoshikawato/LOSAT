#!/usr/bin/env python3
"""Build and run a comparison-only NCBI WordFinder trace probe."""
from __future__ import annotations

import argparse
import hashlib
import os
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
NCBI = Path("/home/kawato/micromamba/bin/tblastn")
NCBI_SHA256 = "e3956f1e107a30439d56c8f72fae4267a7d62ebcc16e4f99d4baf2bdf00402e0"
FRAMES = (1, 2, 3, -1, -2, -3)
FIELDS = "6 qseqid sseqid score qstart qend sstart send sframe qseq sseq"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("fixture_dir", type=Path)
    args = parser.parse_args()
    fixture = args.fixture_dir.resolve()
    assert hashlib.sha256(NCBI.read_bytes()).hexdigest() == NCBI_SHA256
    names = [line[1:].split()[0] for line in
             (fixture / "subjects.fna").read_text().splitlines() if line.startswith(">")]
    with tempfile.TemporaryDirectory(prefix="tlosan-ncbi-trace-") as temp:
        probe = Path(temp) / "ncbi_wordfinder_trace.so"
        subprocess.run(["gcc", "-shared", "-fPIC", "-std=c11", "-O2",
                        "-o", str(probe), str(HERE / "ncbi_wordfinder_trace.c"),
                        "-ldl"], check=True)
        command = [str(NCBI), "-task", "tblastn", "-query", str(fixture / "query.faa"),
                   "-subject", str(fixture / "subjects.fna"), "-db_gencode", "1",
                   "-matrix", "BLOSUM62", "-word_size", "3", "-threshold", "13",
                   "-window_size", "40", "-gapopen", "11", "-gapextend", "1",
                   "-evalue", "10000", "-num_threads", "1", "-comp_based_stats", "0",
                   "-seg", "no", "-sum_stats", "false", "-outfmt", FIELDS]
        env = os.environ.copy()
        env["LD_PRELOAD"] = str(probe)
        result = subprocess.run(command, env=env, capture_output=True, check=True)
        assert result.stdout == (fixture / "raw_isolation_fields.out").read_bytes()
        (fixture / "raw_isolation_wordfinder_trace.out").write_bytes(result.stdout)
        (fixture / "raw_isolation_wordfinder_trace.stderr").write_bytes(result.stderr)
        lines = result.stderr.decode().splitlines()
        calls = [line for line in lines if line.startswith("WORD_FINDER\t")]
        assert len(calls) == 6 * len(names), (len(calls), len(names))
        for index, line in enumerate(calls):
            assert int(line.split("\t")[1]) == index
        rows = ["call\tsubject\tframe\tinit_index\tq_seed\ts_seed\tq_start\ts_start\tlength\traw_score\n"]
        for line in lines:
            if not line.startswith("INIT\t"):
                continue
            parts = line.split("\t")
            assert len(parts) == 9
            call = int(parts[1])
            rows.append(f"{call}\t{names[call // 6]}\t{FRAMES[call % 6]}\t" +
                        "\t".join(parts[2:]) + "\n")
        (fixture / "ncbi_init_trace.tsv").write_text("".join(rows))
        (fixture / "trace_manifest.txt").write_text(
            "Comparison-only LD_PRELOAD probe; never used in LOSAT build or runtime.\n"
            f"NCBI SHA256: {NCBI_SHA256}\n"
            f"Probe C SHA256: {hashlib.sha256((HERE / 'ncbi_wordfinder_trace.c').read_bytes()).hexdigest()}\n"
            "Source commit: 598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4\n"
            "NCBI final outfmt 6 bytes with probe equal the unprobed fixture.\n"
            f"WordFinder calls: {len(calls)}; subjects: {len(names)}; frames per subject: 6\n"
            + "Command: " + repr(command) + "\n")
    paths = ["raw_isolation_wordfinder_trace.out",
             "raw_isolation_wordfinder_trace.stderr",
             "ncbi_init_trace.tsv", "trace_manifest.txt"]
    (fixture / "trace_outputs.sha256").write_text(
        "".join(f"{hashlib.sha256((fixture / name).read_bytes()).hexdigest()}  {name}\n"
                for name in paths))
    print(f"NCBI WordFinder calls {len(calls)}; ungapped init HSPs {len(rows)-1}; output bytes unchanged")


if __name__ == "__main__":
    main()
