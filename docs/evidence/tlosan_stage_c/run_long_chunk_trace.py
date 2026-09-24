#!/usr/bin/env python3
"""Trace NCBI local TBLASTN over the translated 5,000,000-residue chunk edge."""
from __future__ import annotations

import hashlib
import os
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
NCBI_SOURCE = Path("/mnt/c/Users/genom/GitHub/ncbi-blast")
SOURCE_COMMIT = "598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4"
NCBI = Path("/home/kawato/micromamba/bin/tblastn")
NCBI_SHA256 = "e3956f1e107a30439d56c8f72fae4267a7d62ebcc16e4f99d4baf2bdf00402e0"
PREFIX_CODONS = 4_999_950
SUFFIX_CODONS = 250


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_probe(command: list[str], source: Path, temp: Path) -> subprocess.CompletedProcess[bytes]:
    library = temp / (source.stem + ".so")
    subprocess.run(
        ["gcc", "-shared", "-fPIC", "-std=c11", "-O2", "-Wall", "-Wextra",
         "-o", str(library), str(source), "-ldl"],
        check=True,
    )
    env = os.environ.copy()
    env["LD_PRELOAD"] = str(library)
    return subprocess.run(command, capture_output=True, check=True, env=env)


def main() -> None:
    if len(sys.argv) not in (2, 3) or (len(sys.argv) == 3 and sys.argv[2] != "--mask-boundary"):
        raise SystemExit("usage: run_long_chunk_trace.py NEW_OUTPUT_DIR [--mask-boundary]")
    mask_boundary = len(sys.argv) == 3
    output = Path(sys.argv[1]).resolve()
    output.mkdir(parents=True, exist_ok=False)
    assert subprocess.check_output(
        ["git", "-C", str(NCBI_SOURCE), "rev-parse", "HEAD"], text=True
    ).strip() == SOURCE_COMMIT
    assert sha(NCBI) == NCBI_SHA256

    original = (HERE / "run_20260923/subjects.fna").read_text()
    records = {}
    for block in original.split(">"):
        if block.strip():
            name, *lines = block.splitlines()
            records[name] = "".join(lines)
    insert = records["plus1"]
    assert len(insert) == 362
    assert PREFIX_CODONS < 5_000_000 < PREFIX_CODONS + len(insert) // 3
    query = (HERE / "run_20260923/query.faa").read_bytes()
    subject = "ATG" * PREFIX_CODONS + insert + "ATG" * SUFFIX_CODONS
    assert len(subject) // 3 > 5_000_000
    if mask_boundary:
        # Pinned blast_engine.c:283-307: the first soft range's right end
        # equals the second chunk's start at translated offset 4,999,900.
        left, right = 4_999_900 * 3, 5_000_000 * 3
        subject = subject[:left] + subject[left:right].lower() + subject[right:]
    (output / "query.faa").write_bytes(query)
    (output / "subjects.fna").write_text(">chunk_edge\n" + subject + "\n")
    command = [
        str(NCBI), "-task", "tblastn", "-query", str(output / "query.faa"),
        "-subject", str(output / "subjects.fna"), "-db_gencode", "1",
        "-matrix", "BLOSUM62", "-word_size", "3", "-threshold", "13",
        "-window_size", "40", "-gapopen", "11", "-gapextend", "1",
        "-evalue", "10000", "-num_threads", "1", "-comp_based_stats", "0",
        "-seg", "no", "-sum_stats", "false", "-outfmt",
        "6 qseqid sseqid score qstart qend sstart send sframe qseq sseq",
    ]
    if mask_boundary:
        command.append("-lcase_masking")
    plain = subprocess.run(command, capture_output=True, check=True)
    (output / "ncbi_output.out").write_bytes(plain.stdout)
    sources = {
        "wordfinder": HERE / "ncbi_wordfinder_trace.c",
        "gapped": HERE / "ncbi_gapped_trace.c",
        "ranges": HERE / "ncbi_chunk_ranges_trace.c",
    }
    with tempfile.TemporaryDirectory(prefix="tlosan-c-longchunk-") as tmp:
        for name, source in sources.items():
            traced = run_probe(command, source, Path(tmp))
            assert traced.stdout == plain.stdout, f"{name} probe changed NCBI bytes"
            (output / f"{name}.stderr").write_bytes(traced.stderr)
    for name, starts in {
        "frame_chunks.tsv": ("FRAME_CHUNK\t", "WORD_FINDER\t", "INIT\t"),
        "gapped_events.tsv": ("GAPPED_", "MERGE_", "APPEND_", "QUERY_CONTEXT\t",
                              "TRACEBACK_", "TARGET_TRANSLATION\t"),
    }.items():
        lines = (output / ("wordfinder.stderr" if name == "frame_chunks.tsv"
                           else "gapped.stderr")).read_text().splitlines()
        (output / name).write_text(
            "\n".join(line for line in lines if line.startswith(starts)) + "\n"
        )
    range_lines = [line for line in (output / "ranges.stderr").read_text().splitlines()
                   if line.startswith(("RANGE_CALL\t", "RANGE\t", "SET_RANGES\t", "SET_RANGE\t"))]
    (output / "chunk_ranges.tsv").write_text("\n".join(range_lines) + "\n")
    mask_line = (
        "Lowercase mask: nucleotide [14999700,15000000), "
        "translated +1 right edge 4999900.\n"
        if mask_boundary else ""
    )
    (output / "manifest.txt").write_text(
        "Comparison-only pinned NCBI local -subject TBLASTN chunk fixture.\n"
        f"NCBI source commit: {SOURCE_COMMIT}\n"
        f"NCBI binary SHA256: {NCBI_SHA256}\n"
        f"Prefix: ATG x {PREFIX_CODONS} codons; insert: saved plus1 362 nt; "
        f"suffix: ATG x {SUFFIX_CODONS} codons.\n"
        f"{mask_line}"
        f"Subject length: {len(subject)} nt; first frame: {len(subject) // 3} residues.\n"
        f"Input SHA256: query={sha(output / 'query.faa')} "
        f"subject={sha(output / 'subjects.fna')}\n"
        f"Probe source SHA256: {', '.join(f'{name}={sha(source)}' for name, source in sources.items())}\n"
        "Unprobed output equals all three probe outputs byte for byte: yes\n"
        f"Command: {command!r}\n"
    )
    files = ["query.faa", "subjects.fna", "ncbi_output.out", "manifest.txt",
             "wordfinder.stderr", "gapped.stderr", "frame_chunks.tsv",
             "gapped_events.tsv", "ranges.stderr", "chunk_ranges.tsv"]
    (output / "outputs.sha256").write_text(
        "".join(f"{sha(output / name)}  {name}\n" for name in files)
    )
    print("nt", len(subject), "frame chunks",
          sum(line.startswith("FRAME_CHUNK\t") for line in
              (output / "frame_chunks.tsv").read_text().splitlines()),
          "output lines", len(plain.stdout.splitlines()))


if __name__ == "__main__":
    main()
