#!/usr/bin/env python3
"""Pin NCBI SUBJECT_SPLIT_NO_RANGE timing for a fully masked middle chunk."""
from __future__ import annotations

from pathlib import Path
import subprocess
import sys
import tempfile

from run_long_chunk_trace import (
    HERE, NCBI, NCBI_SHA256, NCBI_SOURCE, SOURCE_COMMIT, run_probe, sha,
)


def main() -> None:
    if len(sys.argv) not in (2, 3) or (len(sys.argv) == 3 and sys.argv[2] != "--early-hit"):
        raise SystemExit("usage: run_no_range_middle_trace.py NEW_OUTPUT_DIR [--early-hit]")
    early_hit = len(sys.argv) == 3
    output = Path(sys.argv[1]).resolve()
    output.mkdir(parents=True, exist_ok=False)
    assert subprocess.check_output(
        ["git", "-C", str(NCBI_SOURCE), "rev-parse", "HEAD"], text=True
    ).strip() == SOURCE_COMMIT
    assert sha(NCBI) == NCBI_SHA256
    query = (HERE / "run_20260923/query.faa").read_bytes()
    records = {}
    for block in (HERE / "run_20260923/subjects.fna").read_text().split(">"):
        if block.strip():
            name, *lines = block.splitlines()
            records[name] = "".join(lines)
    insert = records["plus1"].encode()
    assert len(insert) == 362
    if early_hit:
        first = 4_998_450
        second = 10_000_050
        subject = bytearray(
            b"ATG" * first + insert + b"A"
            + b"ATG" * (second - first - 121) + insert + b"ATG" * 250
        )
    else:
        subject = bytearray(b"ATG" * 10_000_050 + insert + b"ATG" * 250)
    left, right = 4_999_000 * 3, 10_000_000 * 3
    subject[left:right] = subject[left:right].lower()
    (output / "query.faa").write_bytes(query)
    (output / "subjects.fna").write_bytes(b">masked_middle\n" + subject + b"\n")
    command = [
        str(NCBI), "-task", "tblastn", "-query", str(output / "query.faa"),
        "-subject", str(output / "subjects.fna"), "-db_gencode", "1",
        "-matrix", "BLOSUM62", "-word_size", "3", "-threshold", "13",
        "-window_size", "40", "-gapopen", "11", "-gapextend", "1",
        "-evalue", "10000", "-num_threads", "1", "-comp_based_stats", "0",
        "-seg", "no", "-sum_stats", "false", "-outfmt", "6",
        "-lcase_masking",
    ]
    plain = subprocess.run(command, capture_output=True, check=True)
    (output / "ncbi_output.out").write_bytes(plain.stdout)
    sources = {
        "wordfinder": HERE / "ncbi_wordfinder_trace.c",
        "gapped": HERE / "ncbi_gapped_trace.c",
        "ranges": HERE / "ncbi_chunk_ranges_trace.c",
        # Pinned aa_ungapped.c:478-505 scans every surviving chunk range.
        "candidate": HERE / "ncbi_candidate_trace.c",
    }
    with tempfile.TemporaryDirectory(prefix="tlosan-no-range-") as tmp:
        traced = {name: run_probe(command, source, Path(tmp))
                  for name, source in sources.items()}
    assert all(result.stdout == plain.stdout for result in traced.values())
    (output / "candidate.stderr").write_bytes(traced["candidate"].stderr)
    lines = [line for line in traced["wordfinder"].stderr.decode().splitlines()
             if line.startswith(("FRAME_CHUNK\t", "WORD_FINDER\t", "INIT\t"))]
    (output / "frame_chunks.tsv").write_text("\n".join(lines) + "\n")
    range_lines = [line for line in traced["ranges"].stderr.decode().splitlines()
                   if line.startswith(("RANGE_CALL\t", "RANGE\t", "SET_RANGES\t", "SET_RANGE\t"))]
    (output / "chunk_ranges.tsv").write_text("\n".join(range_lines) + "\n")
    gapped_lines = [line for line in traced["gapped"].stderr.decode().splitlines()
                    if line.startswith(("GAPPED_", "MERGE_", "APPEND_", "QUERY_CONTEXT\t"))]
    (output / "gapped_events.tsv").write_text("\n".join(gapped_lines) + "\n")
    (output / "manifest.txt").write_text(
        "Comparison-only pinned NCBI TBLASTN masked middle-chunk oracle.\n"
        f"NCBI source commit: {SOURCE_COMMIT}\n"
        f"NCBI binary SHA256: {NCBI_SHA256}\n"
        + ("Subject: saved plus1 362 nt + A pad at codon 4998450, second saved plus1 at "
         "codon 10000050, ATG elsewhere; " if early_hit else
         "Subject: ATG x 10000050 codons, saved plus1 362 nt, ATG x 250; ")
        + "lowercase nucleotide [14997000,30000000), translated +1 [4999000,10000000).\n"
        f"Input SHA256: query={sha(output / 'query.faa')} "
        f"subject={sha(output / 'subjects.fna')}\n"
        f"Probe source SHA256: {', '.join(f'{name}={sha(source)}' for name, source in sources.items())}\n"
        "All four probed final outputs equal unprobed bytes: yes\n"
        f"Command: {command!r}\n"
    )
    names = ("query.faa", "ncbi_output.out", "frame_chunks.tsv", "chunk_ranges.tsv", "gapped_events.tsv", "candidate.stderr", "manifest.txt")
    (output / "retained.sha256").write_text(
        "".join(f"{sha(output / name)}  {name}\n" for name in names)
    )
    print("subject nt", len(subject), "WordFinder calls",
          sum(line.startswith("WORD_FINDER\t") for line in lines))


if __name__ == "__main__":
    main()
