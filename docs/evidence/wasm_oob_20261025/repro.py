#!/usr/bin/env python3
"""Run the n16 threaded-Wasm blastx command of check_wasm_threading.py many times, arms interleaved.

Usage: repro.py --wasm W --out DIR --arms default,prespawn0,refresh --runs N [--threads 16]
The command is check_wasm_threading.py's `format/blastx-blastx-query-2/threaded-n16` run
(prefix "threaded", outfmt "6 std qframe sframe btop"), with its environment.
"""
import argparse, json, os, subprocess, sys, time
from pathlib import Path

ARMS = {"default": {}, "prespawn0": {"LOSAT_WASI_PRESPAWN": "0"}, "refresh": {"OOB_REFRESH": "1"}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--wasm", type=Path, required=True); ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--arms", default="default"); ap.add_argument("--runs", type=int, default=100)
    ap.add_argument("--threads", type=int, default=16); ap.add_argument("--node", default="node")
    a = ap.parse_args()
    root = Path.cwd(); sys.path.insert(0, str(root / "LOSAT" / "tests"))
    from check_wasm_threading import fixtures
    a.out.mkdir(parents=True, exist_ok=True)
    inputs = fixtures(a.out / "fixtures")
    env0 = {k: v for k, v in os.environ.items() if not k.startswith(("LOSAT_", "RAYON_")) and k != "BL2SEQ_LEGACY"}
    env0.update(NODE_NO_WARNINGS="1", LC_ALL="C")
    arms = a.arms.split(","); stats = {arm: {"runs": 0, "fail": 0, "oob": 0, "other": [], "probes": []} for arm in arms}
    rows = open(a.out / "runs.tsv", "w"); rows.write("arm\trun\texit\twall_s\toob\n")
    result = a.out / "result.out"
    for i in range(a.runs):
        for arm in arms:
            argv = [a.node, "LOSAT/tests/run_losat_wasi_threads.js", str(a.wasm.resolve()), "blastx",
                    "-query", str(inputs["blastx-query"]), "-subject", str(inputs["aa3"]),
                    "-outfmt", "6 std qframe sframe btop", "-out", str(result), "-num_threads", str(a.threads)]
            t = time.monotonic()
            try:
                r = subprocess.run(argv, cwd=root, env={**env0, **ARMS[arm]}, capture_output=True, timeout=120)
                code, err = r.returncode, r.stderr.decode(errors="replace")
            except subprocess.TimeoutExpired as e:
                code, err = "timeout", (e.stderr or b"").decode(errors="replace")
            wall = time.monotonic() - t; oob = "memory access out of bounds" in err
            s = stats[arm]; s["runs"] += 1
            if code != 0:
                s["fail"] += 1; s["oob"] += oob
                (a.out / f"fail-{arm}-{i}.txt").write_text(err)
                s["probes"] += [l[6:] for l in err.splitlines() if l.startswith("PROBE ")]
                if not oob: s["other"].append({"run": i, "exit": code, "head": err[:300]})
            rows.write(f"{arm}\t{i}\t{code}\t{wall:.2f}\t{int(oob)}\n"); rows.flush()
    summary = {"node": subprocess.run([a.node, "--version"], capture_output=True, text=True).stdout.strip(),
               "nproc": os.cpu_count(), "threads": a.threads, "arms": stats}
    (a.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


main()
