#!/usr/bin/env python3
"""Non-regression timing for the LOSAT Web engine refactors (plan §6.2, V-PERF).

For each named case (one fixture and output format), runs native LOSAT and the
serial and threaded command-WASI modules with one untimed warmup and exactly three
timed repetitions (the AGENTS.md benchmark protocol), and records the median, the
full range, the output SHA-256 and the start and end time of the run. Run it with the
binaries built before a change and after it, alternately, then `compare` each pair of
JSON files: every median must stay within +5% of the baseline median, and every
output hash must be unchanged.

Usage:
  measure_perf.py run --native BIN --serial WASM --threaded WASM --out FILE.json [--cases blastp,...]
  measure_perf.py compare BEFORE.json AFTER.json
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
MODES = (("native", 1), ("serial-wasi", 1), ("threaded-wasi", 4))
LIMIT = 1.05


def command(mode: str, binaries: dict[str, str], argv: list[str], threads: int, out: Path) -> list[str]:
    tail = argv + ["-num_threads", str(threads), "-out", str(out)]
    if mode == "native":
        return [binaries["native"]] + tail
    runner = "run_losat_wasi.js" if mode == "serial-wasi" else "run_losat_wasi_threads.js"
    wasm = binaries["serial"] if mode == "serial-wasi" else binaries["threaded"]
    return ["node", str(TESTS / runner), wasm] + tail


def measure(argv: list[str]) -> float:
    start = time.perf_counter()
    subprocess.run(argv, check=True, capture_output=True)
    return time.perf_counter() - start


def command_run(args: argparse.Namespace) -> int:
    binaries = {"native": args.native, "serial": args.serial, "threaded": args.threaded}
    cases = args.cases.split(",")
    now = lambda: datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
    results = {"binaries": {name: hashlib.sha256(Path(path).read_bytes()).hexdigest()
                            for name, path in binaries.items()},
               "started_at": now(), "cases": []}
    with tempfile.TemporaryDirectory() as tmp:
        for case in cases:
            fixture = FIXTURES[case]
            for mode, threads in MODES:
                out = Path(tmp) / f"{case}.{mode}.out"
                argv = command(mode, binaries, fixture, threads, out)
                measure(argv)  # untimed warmup
                samples = [measure(argv) for _ in range(3)]
                results["cases"].append({
                    "case": case, "mode": mode, "threads": threads,
                    "argv": argv[argv.index(fixture[0]):], "samples_s": samples,
                    "median_s": statistics.median(samples), "min_s": min(samples), "max_s": max(samples),
                    "output_sha256": hashlib.sha256(out.read_bytes()).hexdigest(),
                })
                print(f"{case} {mode}: median {statistics.median(samples):.3f}s "
                      f"[{min(samples):.3f}, {max(samples):.3f}]")
    results["finished_at"] = now()
    Path(args.out).write_text(json.dumps(results, indent=2) + "\n")
    return 0


def command_compare(args: argparse.Namespace) -> int:
    before = {(c["case"], c["mode"]): c for c in json.loads(Path(args.before).read_text())["cases"]}
    after = {(c["case"], c["mode"]): c for c in json.loads(Path(args.after).read_text())["cases"]}
    failures = 0
    for key in sorted(before):
        if key not in after:
            print(f"missing after: {key}")
            failures += 1
            continue
        ratio = after[key]["median_s"] / before[key]["median_s"]
        same = after[key]["output_sha256"] == before[key]["output_sha256"]
        verdict = "ok" if ratio <= LIMIT and same else "FAIL"
        failures += verdict == "FAIL"
        print(f"{verdict} {key[0]} {key[1]}: median {before[key]['median_s']:.3f}s -> "
              f"{after[key]['median_s']:.3f}s (x{ratio:.3f}), output {'same' if same else 'DIFFERENT'}")
    return 1 if failures else 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="action", required=True)
    run = sub.add_parser("run")
    run.add_argument("--native", required=True)
    run.add_argument("--serial", required=True)
    run.add_argument("--threaded", required=True)
    run.add_argument("--out", required=True)
    run.add_argument("--cases", default="blastp,blastp-fmt0")
    compare = sub.add_parser("compare")
    compare.add_argument("before")
    compare.add_argument("after")
    args = parser.parse_args()
    return command_run(args) if args.action == "run" else command_compare(args)


if __name__ == "__main__":
    sys.exit(main())
