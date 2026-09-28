#!/usr/bin/env python3
"""Non-regression timing for the LOSAT Web engine refactors (plan §6.2, V-PERF).

For each named case (one fixture and output format) and each execution mode (native,
serial command-WASI, threaded command-WASI), runs the binaries built before a change
and after it with one untimed warmup each, then exactly three timed repetitions each
(the AGENTS.md benchmark protocol). The two builds alternate on every repetition
(before, after, then after, before, ...), so a drift of the machine over the run
affects both alike. Records every sample, the median and full range per side, the
output SHA-256 of both builds, and the start and end time of the run.

`check` passes when every after median is within +5% of its before median and every
pair of output hashes is equal.

Usage:
  measure_perf.py run --before NATIVE,SERIAL,THREADED --after NATIVE,SERIAL,THREADED \\
      --out FILE.json [--cases blastp,...] [--repeat N]
  measure_perf.py check FILE.json

Use --repeat above 3 only for an investigation of inconclusive results (AGENTS.md).
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
TESTS = REPO / "LOSAT" / "tests"
FASTA = TESTS / "fasta"
TBLASTN_BENCHMARK = REPO / "docs" / "evidence" / "tlosan_stage_g" / "benchmark"
BLASTP = ["blastp", "-query", str(FASTA / "SicyWSV.faa"), "-subject", str(FASTA / "PajaWSV.faa"),
          "-max_hsps", "1"]
# The Stage G TBLASTN benchmark fixture (docs/evidence/tlosan_stage_g/benchmark/environment.json).
TBLASTN = ["tblastn", "-task", "tblastn", "-query", str(TBLASTN_BENCHMARK / "first_AvCLPV_protein.faa"),
           "-subject", str(FASTA / "AvCLPV.fasta")]
# Named cases: one representative fixture per program, in the tabular format and, where
# the outfmt 0 path changed, in outfmt 0. The argv is the same in every mode.
FIXTURES = {
    "blastp": BLASTP + ["-outfmt", "6"],
    "blastp-fmt0": BLASTP + ["-outfmt", "0"],
    "tblastn": TBLASTN + ["-outfmt", "6"],
    "tblastn-fmt0": TBLASTN + ["-outfmt", "0"],
    "blastn": ["blastn", "-task", "megablast", "-query", str(FASTA / "AP027152.fasta"),
               "-subject", str(FASTA / "LC738884.fasta"), "-outfmt", "6"],
    "tblastx": ["tblastx", "-query", str(FASTA / "LC738884.fasta"), "-subject", str(FASTA / "LC741431.fasta"),
                "-outfmt", "6"],
}
# (mode, index in the NATIVE,SERIAL,THREADED list, -num_threads)
MODES = (("native", 0, 1), ("serial-wasi", 1, 1), ("threaded-wasi", 2, 4))
SIDES = ("before", "after")
LIMIT = 1.05


def command(mode: str, binary: str, argv: list[str], threads: int, out: Path) -> list[str]:
    tail = argv + ["-num_threads", str(threads), "-out", str(out)]
    if mode == "native":
        return [binary] + tail
    runner = "run_losat_wasi.js" if mode == "serial-wasi" else "run_losat_wasi_threads.js"
    return ["node", str(TESTS / runner), binary] + tail


def measure(argv: list[str]) -> float:
    start = time.perf_counter()
    subprocess.run(argv, check=True, capture_output=True)
    return time.perf_counter() - start


def now() -> str:
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")


def command_run(args: argparse.Namespace) -> int:
    builds = {"before": args.before.split(","), "after": args.after.split(",")}
    for side, paths in builds.items():
        if len(paths) != 3:
            raise SystemExit(f"--{side} needs NATIVE,SERIAL,THREADED")
    results = {
        "binaries": {side: [hashlib.sha256(Path(path).read_bytes()).hexdigest() for path in paths]
                     for side, paths in builds.items()},
        "repeat": args.repeat, "started_at": now(), "cases": [],
    }
    with tempfile.TemporaryDirectory() as tmp:
        for case in args.cases.split(","):
            fixture = FIXTURES[case]
            for mode, index, threads in MODES:
                outs = {side: Path(tmp) / f"{case}.{mode}.{side}.out" for side in SIDES}
                argvs = {side: command(mode, builds[side][index], fixture, threads, outs[side])
                         for side in SIDES}
                for side in SIDES:
                    measure(argvs[side])  # untimed warmup
                samples: dict[str, list[float]] = {side: [] for side in SIDES}
                for repetition in range(args.repeat):
                    for side in (SIDES if repetition % 2 == 0 else SIDES[::-1]):
                        samples[side].append(measure(argvs[side]))
                entry = {"case": case, "mode": mode, "threads": threads,
                         "argv": argvs["after"][argvs["after"].index(fixture[0]):]}
                for side in SIDES:
                    entry[side] = {
                        "samples_s": samples[side], "median_s": statistics.median(samples[side]),
                        "min_s": min(samples[side]), "max_s": max(samples[side]),
                        "output_sha256": hashlib.sha256(outs[side].read_bytes()).hexdigest(),
                    }
                results["cases"].append(entry)
                print(verdict(entry), flush=True)
    results["finished_at"] = now()
    Path(args.out).write_text(json.dumps(results, indent=2) + "\n")
    return 0


def verdict(entry: dict) -> str:
    before, after = entry["before"], entry["after"]
    ratio = after["median_s"] / before["median_s"]
    same = after["output_sha256"] == before["output_sha256"]
    word = "ok" if ratio <= LIMIT and same else "FAIL"
    return (f"{word} {entry['case']} {entry['mode']}: median {before['median_s']:.3f}s "
            f"[{before['min_s']:.3f}, {before['max_s']:.3f}] -> {after['median_s']:.3f}s "
            f"[{after['min_s']:.3f}, {after['max_s']:.3f}] (x{ratio:.3f}), "
            f"output {'same' if same else 'DIFFERENT'}")


def command_check(args: argparse.Namespace) -> int:
    lines = [verdict(entry) for entry in json.loads(Path(args.file).read_text())["cases"]]
    print("\n".join(lines))
    return 1 if any(line.startswith("FAIL") for line in lines) else 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="action", required=True)
    run = sub.add_parser("run")
    run.add_argument("--before", required=True)
    run.add_argument("--after", required=True)
    run.add_argument("--out", required=True)
    run.add_argument("--cases", default="blastp,blastp-fmt0")
    run.add_argument("--repeat", type=int, default=3)
    check = sub.add_parser("check")
    check.add_argument("file")
    args = parser.parse_args()
    return command_run(args) if args.action == "run" else command_check(args)


if __name__ == "__main__":
    sys.exit(main())
