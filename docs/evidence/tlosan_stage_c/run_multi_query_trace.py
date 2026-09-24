#!/usr/bin/env python3
"""Pin NCBI TBLASTN query contexts, candidates, gapped HSPs, and append order."""
from __future__ import annotations

import hashlib
import os
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
NCBI_SOURCE = Path("/mnt/c/Users/genom/GitHub/ncbi-blast")
SOURCE_COMMIT = "598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4"
NCBI = Path("/home/kawato/micromamba/bin/tblastn")
NCBI_SHA256 = "e3956f1e107a30439d56c8f72fae4267a7d62ebcc16e4f99d4baf2bdf00402e0"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def probe_run(command: list[str], source: Path, directory: Path) -> subprocess.CompletedProcess[bytes]:
    library = directory / (source.stem + ".so")
    subprocess.run(
        ["gcc", "-shared", "-fPIC", "-std=c11", "-O2", "-Wall", "-Wextra",
         "-o", str(library), str(source), "-ldl"],
        check=True,
    )
    env = os.environ.copy()
    env["LD_PRELOAD"] = str(library)
    return subprocess.run(command, capture_output=True, check=True, env=env)


def main() -> None:
    if len(sys.argv) not in (2, 3):
        raise SystemExit("usage: run_multi_query_trace.py NEW_OUTPUT_DIR [--blosum45-word2]")
    mode = sys.argv[2] if len(sys.argv) == 3 else ""
    assert mode in ("", "--blosum45-word2")
    output = Path(sys.argv[1]).resolve()
    output.mkdir(parents=True, exist_ok=False)
    assert subprocess.check_output(
        ["git", "-C", str(NCBI_SOURCE), "rev-parse", "HEAD"], text=True
    ).strip() == SOURCE_COMMIT
    assert sha(NCBI) == NCBI_SHA256

    original = (HERE / "run_20260923/query.faa").read_text()
    protein = "".join(original.splitlines()[1:])
    assert len(protein) == 120
    # Pinned blast_query_info.c:68-96 gives one query context per
    # protein query; offset/length are captured rather than inferred.
    queries = [("full", protein), ("internal70", protein[25:95]), ("no_hit_x", "X" * 120)]
    (output / "query.faa").write_text(
        "".join(f">{name}\n{sequence}\n" for name, sequence in queries)
    )
    if mode:
        records = {}
        for block in (HERE / "run_20260923/subjects.fna").read_text().split(">"):
            if block.strip():
                name, *lines = block.splitlines()
                records[name] = "".join(lines)
        (output / "subjects.fna").write_text(">plus1\n" + records["plus1"] + "\n")
        matrix, word_size, threshold, window, gap_open, gap_extend = (
            "BLOSUM45", "2", "16", "60", "14", "2"
        )
        subject_label = "saved plus1"
    else:
        (output / "subjects.fna").write_bytes(
            (HERE / "six_frame_merged_20260924/subjects.fna").read_bytes()
        )
        matrix, word_size, threshold, window, gap_open, gap_extend = (
            "BLOSUM62", "3", "13", "40", "11", "1"
        )
        subject_label = "saved six-frame sequence"
    command = [
        str(NCBI), "-task", "tblastn", "-query", str(output / "query.faa"),
        "-subject", str(output / "subjects.fna"), "-db_gencode", "1",
        "-matrix", matrix, "-word_size", word_size, "-threshold", threshold,
        "-window_size", window, "-gapopen", gap_open, "-gapextend", gap_extend,
        "-evalue", "10000", "-num_threads", "1", "-comp_based_stats", "0",
        "-seg", "no", "-sum_stats", "false", "-outfmt",
        "6 qseqid sseqid score qstart qend sstart send sframe qseq sseq",
    ]
    plain = subprocess.run(command, capture_output=True, check=True)
    (output / "ncbi_output.out").write_bytes(plain.stdout)
    with tempfile.TemporaryDirectory(prefix="tlosan-c-multiquery-") as tmp:
        temp = Path(tmp)
        sources = {
            "gapped": HERE / "ncbi_gapped_trace.c",
            "wordfinder": HERE / "ncbi_wordfinder_trace.c",
            "candidate": HERE / "ncbi_candidate_trace.c",
        }
        for name, source in sources.items():
            traced = probe_run(command, source, temp)
            assert traced.stdout == plain.stdout, f"{name} probe changed NCBI output"
            (output / f"{name}.stderr").write_bytes(traced.stderr)

    gapped_lines = (output / "gapped.stderr").read_text().splitlines()
    selected = {
        "query_context.tsv": ("QUERY_INFO\t", "QUERY_CONTEXT\t"),
        "append_events.tsv": ("APPEND_", "MERGE_"),
        "gapped_events.tsv": ("GAPPED_",),
        "traceback_events.tsv": ("TRACEBACK_", "TARGET_TRANSLATION\t",
                                 "CONTAINS\t", "HSP_TEST\t", "ENDPOINT_"),
    }
    for filename, prefixes in selected.items():
        lines = [line for line in gapped_lines if line.startswith(prefixes)]
        (output / filename).write_text("\n".join(lines) + "\n")
    (output / "manifest.txt").write_text(
        "Comparison-only pinned NCBI TBLASTN local -subject oracle.\n"
        f"NCBI source commit: {SOURCE_COMMIT}\n"
        f"NCBI binary SHA256: {NCBI_SHA256}\n"
        f"LOSAT HEAD: {subprocess.check_output(['git', '-C', str(ROOT), 'rev-parse', 'HEAD'], text=True).strip()}\n"
        "Queries: full=120 aa, internal70=70 aa (full[25:95]), no_hit_x=120 aa.\n"
        f"Subject: {subject_label}, {sum(len(line) for line in (output / 'subjects.fna').read_text().splitlines()[1:])} nt.\n"
        f"Probe source SHA256: {', '.join(f'{name}={sha(source)}' for name, source in sources.items())}\n"
        "Unprobed output equals all three probe outputs byte for byte: yes\n"
        f"Command: {command!r}\n"
    )
    files = ["query.faa", "subjects.fna", "ncbi_output.out", "manifest.txt"]
    files += [f"{name}.stderr" for name in sources]
    files += list(selected)
    (output / "outputs.sha256").write_text(
        "".join(f"{sha(output / name)}  {name}\n" for name in files)
    )
    print("queries", len(queries), "contexts",
          sum(line.startswith("QUERY_CONTEXT\t") for line in gapped_lines),
          "append calls", sum(line.startswith("APPEND_INPUT\t") for line in gapped_lines),
          "gapped HSPs", sum(line.startswith("GAPPED_HSP\t") for line in gapped_lines))


if __name__ == "__main__":
    main()
