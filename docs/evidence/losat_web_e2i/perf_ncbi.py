#!/usr/bin/env python3
"""Time LOSAT's dc-megablast and blastn-short against NCBI BLAST+ 2.17.0 (Session SD).

The tasks are new in LOSAT, so there is no before/after comparison; this records the
ratio of LOSAT's wall time to NCBI's. As AGENTS.md's benchmark protocol requires, NCBI
searches the subject as a database (`makeblastdb -dbtype nucl` before the timing; its
command, version and time are recorded separately) and LOSAT searches it with -subject;
each case has one untimed warmup and three timed runs per side, alternating, and the
median and the full range are reported. The outputs are not compared here (the fixtures
and sweeps do that); the HSP lines of LOSAT's outfmt 6 are counted.

Usage: perf_ncbi.py --losat LOSAT --bin-dir NCBI_BIN --work DIR --out JSON [--repeat 3]
"""
from __future__ import annotations

import argparse
import json
import statistics
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
FASTA = REPO / "LOSAT/tests/fasta"
# (name, task, query, subject, threads, extra arguments)
CASES = [
    ("dc-edl933-sakai", "dc-megablast", "EDL933.fna", "Sakai.fna", 1, []),
    ("dc-edl933-sakai-t4", "dc-megablast", "EDL933.fna", "Sakai.fna", 4, []),
    ("dc-lc738874-lc738875", "dc-megablast", "LC738874.fasta", "LC738875.fasta", 1, []),
    ("dc-two21-lc738874-lc738875", "dc-megablast", "LC738874.fasta", "LC738875.fasta", 1,
     ["-template_type", "coding_and_optimal", "-template_length", "21"]),
    ("short-lc738874-lc738875", "blastn-short", "LC738874.fasta", "LC738875.fasta", 1, []),
    ("short-primers-lc738875", "blastn-short", "outfmt0/short_query.fasta", "LC738875.fasta", 1, []),
]


def timed(argv: list[str], cwd: Path) -> tuple[float, bytes]:
    start = time.perf_counter()
    run = subprocess.run(argv, cwd=cwd, capture_output=True, check=True)
    return time.perf_counter() - start, run.stdout


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--losat", type=Path, required=True)
    parser.add_argument("--bin-dir", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--repeat", type=int, default=3)
    args = parser.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    makeblastdb = args.bin_dir / "makeblastdb"
    record = {"losat": str(args.losat), "ncbi": str(args.bin_dir / "blastn"),
              "ncbi_version": subprocess.run([args.bin_dir / "blastn", "-version"], capture_output=True,
                                             text=True).stdout.strip(), "databases": {}, "cases": []}
    for subject in sorted({case[3] for case in CASES}):
        db = args.work / Path(subject).name
        start = time.perf_counter()
        command = [str(makeblastdb), "-in", str(FASTA / subject), "-dbtype", "nucl", "-out", str(db)]
        subprocess.run(command, capture_output=True, check=True)
        record["databases"][subject] = {"command": command, "seconds": round(time.perf_counter() - start, 3)}
    for name, task, query, subject, threads, extra in CASES:
        common = ["-task", task, "-query", str(FASTA / query), *extra, "-outfmt", "6", "-num_threads", str(threads)]
        losat = [str(args.losat), "blastn", *common, "-subject", str(FASTA / subject)]
        ncbi = [str(args.bin_dir / "blastn"), *common, "-db", str(args.work / Path(subject).name)]
        timed(losat, args.work)
        timed(ncbi, args.work)
        samples = {"losat": [], "ncbi": []}
        lines = None
        for rep in range(args.repeat):
            for side in (("losat", "ncbi") if rep % 2 == 0 else ("ncbi", "losat")):
                seconds, out = timed(losat if side == "losat" else ncbi, args.work)
                samples[side].append(round(seconds, 3))
                if side == "losat":
                    lines = out.count(b"\n")
        medians = {side: statistics.median(values) for side, values in samples.items()}
        row = {"case": name, "threads": threads, "losat_hsp_lines": lines, "samples": samples, "median": medians,
               "ratio_losat_to_ncbi": round(medians["losat"] / medians["ncbi"], 3)}
        record["cases"].append(row)
        print(f"{name}: LOSAT {medians['losat']:.3f}s {samples['losat']}, NCBI {medians['ncbi']:.3f}s "
              f"{samples['ncbi']} (x{row['ratio_losat_to_ncbi']}), {lines} HSP lines", flush=True)
    args.out.write_text(json.dumps(record, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
