#!/usr/bin/env python3
"""Generate fixed TBLASTN subjects and fresh, comparison-only NCBI output.

NCBI c++/src/algo/blast/core/blast_engine.c:804-855:
    for (context=first_context; context<=last_context; context++) {
        subject->frame = BLAST_ContextToFrame(eBlastTypeBlastx, context);
        status = s_BlastSearchEngineOneContext(...);
        Blast_HSPListAppend(&hsp_list_for_chunks, &hsp_list_out, kHspNumMax);
    }
"""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[3]
STAGE_A = ROOT / "docs/evidence/tlosan_stage_a/fixtures"
NCBI = Path("/home/kawato/micromamba/bin/tblastn")
PINNED_SHA256 = "e3956f1e107a30439d56c8f72fae4267a7d62ebcc16e4f99d4baf2bdf00402e0"


def fasta_sequence(path: Path) -> str:
    return "".join(line.strip() for line in path.read_text().splitlines()
                   if line and not line.startswith(">"))


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    path.write_text("".join(f">{name}\n{sequence}\n" for name, sequence in records))


def run() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path, help="new directory; existing paths are refused")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    output = args.output.resolve()
    assert hashlib.sha256(NCBI.read_bytes()).hexdigest() == PINNED_SHA256

    core = fasta_sequence(STAGE_A / "subject_code1.fna")
    assert len(core) == 360
    query = fasta_sequence(STAGE_A / "query.faa")
    assert len(query) == 120
    write_fasta(output / "query.faa", [("q1", query)])

    # The offset before the known 360-nt ORF selects the frame. Reverse
    # complements of the same constructs exercise contexts 3, 4, and 5.
    # NCBI blast_util.c:1080-1101 translates context order +1,+2,+3,-1,-2,-3.
    plus = [(f"plus{frame}", "A" * (frame - 1) + core + "C" * (3 - frame))
            for frame in (1, 2, 3)]
    minus = [(f"minus{frame}", reverse_complement(sequence))
             for frame, (_, sequence) in enumerate(plus, start=1)]
    ambiguous = core[:150] + "N" + core[151:]
    stopped = core[:120] + "TAA" + core[123:]
    records = plus + minus + [
        ("partial_codon", core + "AC"),
        ("ambiguous", ambiguous),
        ("internal_stop", stopped),
        ("low_complexity", "AAA" * 120),
        ("no_hit", "TAA" * 120),
        ("tie_a", core),
        ("tie_b", core),
    ]
    write_fasta(output / "subjects.fna", records)
    write_fasta(output / "query_low.faa", [("q_low", "K" * 120)])
    write_fasta(output / "low_subject.fna", [("low_complexity", "AAA" * 120)])
    (output / "cases.tsv").write_text(
        "subject\tfeature\tlength_nt\n" +
        "".join(f"{name}\t{('frame ' + name[-1]) if name.startswith(('plus', 'minus')) else name}\t{len(seq)}\n"
                for name, seq in records))

    fields = "6 qseqid sseqid score qstart qend sstart send sframe qseq sseq"
    common = [str(NCBI), "-task", "tblastn", "-query", str(output / "query.faa"),
              "-subject", str(output / "subjects.fna"), "-db_gencode", "1",
              "-matrix", "BLOSUM62", "-word_size", "3", "-threshold", "13",
              "-window_size", "40", "-gapopen", "11", "-gapextend", "1",
              "-evalue", "10000", "-num_threads", "1"]
    runs = [("raw_isolation", ["-comp_based_stats", "0", "-seg", "no",
                               "-sum_stats", "false"]),
            ("default_profile", ["-comp_based_stats", "2", "-seg", "12 2.2 2.5",
                                 "-sum_stats", "true"])]
    manifest = [f"NCBI binary: {NCBI}", f"NCBI sha256: {PINNED_SHA256}",
                "NCBI source commit: 598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4",
                "Target: local -subject, native Linux, one thread",
                "Output 6 includes raw score and sframe; output 0 is retained separately."]
    for name, options in runs:
        for label, fmt in (("fields", fields), ("pairwise", "0")):
            cmd = common + options + ["-outfmt", fmt]
            result = subprocess.run(cmd, capture_output=True, check=False)
            (output / f"{name}_{label}.out").write_bytes(result.stdout)
            (output / f"{name}_{label}.stderr").write_bytes(result.stderr)
            if result.returncode:
                raise RuntimeError(f"NCBI {name}/{label} exited {result.returncode}: {result.stderr.decode()}")
            manifest.append(" ".join(cmd))
    # NCBI blast_filter.c:337-370 filters the protein query before lookup.
    # The same low-complexity pair probes SEG's seed suppression separately.
    for label, seg_value in (("seg_on", "12 2.2 2.5"), ("seg_off", "no")):
        cmd = [str(NCBI), "-task", "tblastn", "-query", str(output / "query_low.faa"),
               "-subject", str(output / "low_subject.fna"), "-db_gencode", "1",
               "-comp_based_stats", "0", "-seg", seg_value, "-sum_stats", "false",
               "-evalue", "10000", "-outfmt", fields]
        result = subprocess.run(cmd, capture_output=True, check=False)
        (output / f"low_query_{label}.out").write_bytes(result.stdout)
        (output / f"low_query_{label}.stderr").write_bytes(result.stderr)
        if result.returncode:
            raise RuntimeError(f"NCBI low-complexity {label} exited {result.returncode}: {result.stderr.decode()}")
        manifest.append(" ".join(cmd))
    (output / "manifest.txt").write_text("\n".join(manifest) + "\n")
    files = sorted(path for path in output.iterdir() if path.is_file())
    (output / "inputs_outputs.sha256").write_text(
        "".join(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}\n"
                for path in files))
    print(f"NCBI six-frame oracle generated in {output}")


if __name__ == "__main__":
    run()
