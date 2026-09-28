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
        raise SystemExit("usage: run_multi_query_trace.py NEW_OUTPUT_DIR [--blosum45-word2|--seg-hard|--seg-cross|--seg-soft|--lcase-query|--lcase-soft|--seg-lcase-overlap|--seg-lcase-overlap-soft]")
    mode = sys.argv[2] if len(sys.argv) == 3 else ""
    assert mode in ("", "--blosum45-word2", "--seg-hard", "--seg-cross", "--seg-soft", "--lcase-query", "--lcase-soft", "--seg-lcase-overlap", "--seg-lcase-overlap-soft")
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
    if mode == "--seg-hard":
        queries = [("masked_prefix", "K" * 40 + protein)]
    elif mode in ("--seg-lcase-overlap", "--seg-lcase-overlap-soft"):
        sequence = "K" * 40 + protein
        queries = [("seg_lcase_overlap", sequence[:35] + sequence[35:45].lower() + sequence[45:])]
    elif mode in ("--seg-cross", "--seg-soft"):
        queries = [("masked_internal", protein[:60] + "K" * 18 + protein[60:])]
    elif mode in ("--lcase-query", "--lcase-soft"):
        queries = [("lowercase_query", protein[:40] + protein[40:60].lower() + protein[60:])]
    else:
        queries = [("full", protein), ("internal70", protein[25:95]),
                   ("no_hit_x", "X" * 120)]
    (output / "query.faa").write_text(
        "".join(f">{name}\n{sequence}\n" for name, sequence in queries)
    )
    if mode in ("--blosum45-word2", "--seg-hard", "--seg-cross", "--seg-soft", "--lcase-query", "--lcase-soft", "--seg-lcase-overlap", "--seg-lcase-overlap-soft"):
        records = {}
        for block in (HERE / "run_20260923/subjects.fna").read_text().split(">"):
            if block.strip():
                name, *lines = block.splitlines()
                records[name] = "".join(lines)
        subject_sequence = records["plus1"]
        if mode in ("--seg-cross", "--seg-soft"):
            assert len(subject_sequence) == 362
            subject_sequence = subject_sequence[:180] + "AAA" * 18 + subject_sequence[180:]
        (output / "subjects.fna").write_text(">plus1\n" + subject_sequence + "\n")
        matrix, word_size, threshold, window, gap_open, gap_extend = (
            ("BLOSUM45", "2", "16", "60", "14", "2")
            if mode == "--blosum45-word2" else
            ("BLOSUM62", "3", "13", "40", "11", "1")
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
        "-seg", "12 2.2 2.5" if mode in ("--seg-hard", "--seg-cross", "--seg-soft", "--seg-lcase-overlap", "--seg-lcase-overlap-soft") else "no",
        "-sum_stats", "false", "-outfmt",
        "6 qseqid sseqid score qstart qend sstart send sframe qseq sseq",
    ]
    if mode in ("--seg-soft", "--lcase-soft", "--seg-lcase-overlap-soft"):
        command.extend(["-soft_masking", "true"])
    if mode in ("--lcase-query", "--lcase-soft", "--seg-lcase-overlap", "--seg-lcase-overlap-soft"):
        command.append("-lcase_masking")
    plain = subprocess.run(command, capture_output=True, check=True)
    (output / "ncbi_output.out").write_bytes(plain.stdout)
    with tempfile.TemporaryDirectory(prefix="tlosan-c-multiquery-") as tmp:
        temp = Path(tmp)
        sources = {
            "gapped": HERE / "ncbi_gapped_trace.c",
            "wordfinder": HERE / "ncbi_wordfinder_trace.c",
            "candidate": HERE / "ncbi_candidate_trace.c",
        }
        if mode in ("--seg-hard", "--seg-cross", "--seg-soft", "--lcase-query", "--lcase-soft", "--seg-lcase-overlap", "--seg-lcase-overlap-soft"):
            sources["seg_state"] = HERE / "ncbi_seg_query_state_trace.c"
            sources["context_cutoff"] = HERE / "ncbi_wordfinder_context_cutoff_trace.c"
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
        f"Queries: {[(name, len(seq)) for name, seq in queries]}; SEG={mode in ('--seg-hard', '--seg-cross', '--seg-soft', '--seg-lcase-overlap', '--seg-lcase-overlap-soft')}.\n"
        f"Subject: {subject_label}, {sum(len(line) for line in (output / 'subjects.fna').read_text().splitlines()[1:])} nt.\n"
        f"Probe source SHA256: {', '.join(f'{name}={sha(source)}' for name, source in sources.items())}\n"
        f"Unprobed output equals all {len(sources)} probe outputs byte for byte: yes\n"
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
