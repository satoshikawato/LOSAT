#!/usr/bin/env python3
"""A/B of the threaded Wasm command built with the 1 GiB maximum (A, out-base) and with
65,535 pages (B, out-new), on the same Node host files, plus the same searches with the
heap above 2.4 GiB (B with LOSAT_WASM_MEMORY_INITIAL_PAGES=40000) and natively.

Usage: ab_time.py timing|high [--rounds N] [--case NAME]
  timing: 1 warm-up and N alternating rounds per case and thread count; wall and CPU seconds.
  high:   one run each of native, A, B, B-high; the outputs must have one SHA-256 per case.
Results: logs/ab_timing.tsv, logs/ab_timing_summary.txt, logs/high_memory.tsv
"""
import argparse
import hashlib
import os
import resource
import statistics
import subprocess
import sys
import time
from pathlib import Path

T = Path("/home/kawato/losat-baselines/wasm-mem-20261012")
WT = Path("/home/kawato/losat-work/.worktrees/wasm-mem")
BR = Path("/home/kawato/.cache/losat-work/wasm-mem")
NODE = "/home/kawato/.local/bin/node"
RUNNER = str(WT / "LOSAT/tests/run_losat_wasi_threads.js")
NATIVE = str(BR / "target-native/release/LOSAT")
WASM = {"A": str(BR / "out-base/engine/losat-threaded-command.wasm"),
        "B": str(BR / "out-new/engine/losat-threaded-command.wasm")}
F = "tests/fasta"
CASES = {
    "p03": ["tblastx", "-query", f"{F}/MelaMJNV.fasta", "-subject", f"{F}/PemoMJNVA.fasta", "-outfmt", "6",
            "-query_gencode", "1", "-db_gencode", "1"],
    "d04": ["tblastx", "-query", f"{F}/AP027131.fasta", "-subject", f"{F}/AP027133.fasta", "-outfmt", "6",
            "-query_gencode", "4", "-db_gencode", "4"],
    "blastp_ap": ["blastp", "-query", f"{F}/AP027078.faa", "-subject", f"{F}/AP027133.faa", "-outfmt", "6"],
    "tblastn_ap": ["tblastn", "-query", f"{F}/AP027131.faa", "-subject", f"{F}/AP027133.fasta", "-outfmt", "6",
                   "-db_gencode", "4"],
    "blastn_big": ["blastn", "-task", "blastn", "-query", f"{F}/AP027202.fasta", "-subject", f"{F}/LC738875.fasta",
                   "-outfmt", "6"],
}
# The repository's Node flags for comparisons (LOSAT/tests/README.md): TurboFan from the start.
NODE_FLAGS = ["--no-liftoff", "--no-wasm-tier-up"]
OUT = T / "tmp" / "ab"


def environment(extra=None):
    env = {k: v for k, v in os.environ.items() if not k.startswith(("LOSAT_", "RAYON_", "NODE_"))}
    env.update(NODE_NO_WARNINGS="1", LC_ALL="C", **(extra or {}))
    return env


def run(variant, case, threads, extra=None, tag=""):
    out = OUT / f"{case}.{variant}{tag}.n{threads}.out"
    argv = [*CASES[case], "-num_threads", str(threads), "-out", str(out)]
    command = [NATIVE, *argv] if variant == "native" else [NODE, *NODE_FLAGS, RUNNER, WASM[variant], *argv]
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    start = time.monotonic()
    done = subprocess.run(command, cwd=WT / "LOSAT", env=environment(extra), capture_output=True, text=True)
    wall = time.monotonic() - start
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    cpu = (after.ru_utime - before.ru_utime) + (after.ru_stime - before.ru_stime)
    digest = hashlib.sha256(out.read_bytes()).hexdigest() if out.exists() else "missing"
    return dict(rc=done.returncode, wall=wall, cpu=cpu, sha=digest, stderr=done.stderr, rows=out.read_bytes().count(b"\n") if out.exists() else 0)


PLAN = {"p03": ((1, 8), 5), "d04": ((8,), 3), "blastp_ap": ((1, 8), 5)}


def timing(rounds, only=None):
    """One case per call when `only` is given, so that the machine-wide lock is held for minutes, not for the whole plan.
    d04 takes 150 s at one thread on this tree (no round-3 switches): 8 threads only, 3 rounds."""
    tsv, text = T / "logs/ab_timing.tsv", T / "logs/ab_timing_summary.txt"
    if not tsv.exists():
        tsv.write_text("case\tthreads\tround\tvariant\trc\twall_s\tcpu_s\tsha256_16\tload1\n")
    for case, (counts, plan_rounds) in PLAN.items():
        if only and case != only:
            continue
        case_rounds = min(rounds, plan_rounds)
        for threads in counts:
            for variant in ("A", "B"):
                run(variant, case, threads)  # warm-up (file cache, guard)
            rows, ratios_wall, ratios_cpu, walls = [], [], [], {"A": [], "B": []}
            shas = set()
            for index in range(case_rounds):
                order = ("A", "B") if index % 2 == 0 else ("B", "A")
                result = {v: run(v, case, threads) for v in order}
                for v in order:
                    r = result[v]
                    rows.append((case, threads, index, v, r["rc"], f"{r['wall']:.3f}", f"{r['cpu']:.3f}", r["sha"][:16],
                                 f"{os.getloadavg()[0]:.1f}"))
                    walls[v].append(r["wall"])
                    shas.add((r["rc"], r["sha"]))
                ratios_wall.append(result["B"]["wall"] / result["A"]["wall"])
                ratios_cpu.append(result["B"]["cpu"] / result["A"]["cpu"])
            line = (f"{case} n{threads} ({case_rounds} rounds): wall median A {statistics.median(walls['A']):.3f} s, B {statistics.median(walls['B']):.3f} s; "
                    f"B/A per round: wall median {statistics.median(ratios_wall):.3f} (min {min(ratios_wall):.3f}, max {max(ratios_wall):.3f}), "
                    f"CPU median {statistics.median(ratios_cpu):.3f} (min {min(ratios_cpu):.3f}, max {max(ratios_cpu):.3f}); "
                    f"outputs {'one SHA-256, rc 0' if len(shas) == 1 and next(iter(shas))[0] == 0 else 'DIFFER ' + str(shas)}")
            print(line, flush=True)
            with tsv.open("a") as handle:
                handle.write("".join("\t".join(map(str, row)) + "\n" for row in rows))
            with text.open("a") as handle:
                handle.write(line + "\n")


def high():
    lines, bad = [], 0
    for case in CASES:
        for threads in (1, 8):
            results = {
                "native": run("native", case, threads),
                "A": run("A", case, threads),
                "B": run("B", case, threads),
                "B-high": run("B", case, threads, {"LOSAT_WASM_MEMORY_INITIAL_PAGES": "40000", "LOSAT_WASI_THREADS_DEBUG": "1"}, ".high"),
            }
            memory = [l for l in results["B-high"]["stderr"].splitlines() if "shared memory:" in l]
            shas = {r["sha"] for r in results.values()}
            ok = len(shas) == 1 and all(r["rc"] == 0 for r in results.values()) and "missing" not in shas
            bad += not ok
            line = "\t".join([case, str(threads), "same" if ok else "DIFFERENT", results["B"]["sha"][:16], str(results["B"]["rows"]),
                              " ".join(f"{k}:rc{r['rc']}:{r['sha'][:8]}:{r['wall']:.1f}s" for k, r in results.items()),
                              memory[0].split("] ")[-1] if memory else "no memory line"])
            print(line, flush=True)
            lines.append(line)
    (T / "logs/high_memory.tsv").write_text("case\tthreads\tverdict\tsha256_16\trows\truns\thigh_memory\n" + "\n".join(lines) + "\n")
    return bad


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("mode", choices=["timing", "high"])
    parser.add_argument("--rounds", type=int, default=5)
    parser.add_argument("--case", choices=sorted(PLAN), help="timing: only this case")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if args.mode == "timing":
        timing(args.rounds, args.case)
        return 0
    return 1 if high() else 0


if __name__ == "__main__":
    sys.exit(main())
