#!/usr/bin/env python3
"""Line WO: trap counts of the threaded command, arms in turn, several cases.
usage: arm_rate.py --arm NAME=HOSTDIR=WASM [--arm ...] --out DIR --runs N [--threads 16] [--node node] [--cases bp,bn]
Run from the repository root. Keeps the first stderr lines of every failed run."""
import argparse, collections, hashlib, json, os, re, subprocess, sys, time
from pathlib import Path

F = "LOSAT/tests/fasta"
CASES = {
    "bn": ["blastn", "-task", "blastn", "-query", f"{F}/AP027131.fasta", "-subject", f"{F}/NZ_CP006932.fasta"],
    "mb": ["blastn", "-task", "megablast", "-query", f"{F}/AP027131.fasta", "-subject", f"{F}/AP027132.fasta"],
    "dc": ["blastn", "-task", "dc-megablast", "-query", f"{F}/AP027131.fasta", "-subject", f"{F}/AP027132.fasta"],
    "bp": ["blastp", "-query", f"{F}/AP027131.faa", "-subject", f"{F}/AP027132.faa"],
}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", action="append", required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--runs", type=int, default=10); ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--node", default="node"); ap.add_argument("--cases", default=",".join(CASES))
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    env0 = {k: v for k, v in os.environ.items() if not k.startswith(("LOSAT_", "RAYON_"))}
    env0.update(NODE_NO_WARNINGS="1", LC_ALL="C")
    rows = open(a.out / "runs.tsv", "w")
    rows.write("case\tarm\trun\texit\twall_s\tmaxrss_kb\toob\tsha16\n")
    stats = collections.defaultdict(lambda: {"runs": 0, "fail": 0, "oob": 0, "shas": collections.Counter(), "frames": collections.Counter(), "rss": [], "wall": []})
    result = a.out / "result.out"; timing = a.out / "time.txt"
    specs = {name: (Path(host), Path(wasm)) for name, host, wasm in (x.split("=") for x in a.arm)}
    arms = list(specs)
    for i in range(a.runs):
        for case in a.cases.split(","):
            for arm in arms[i % len(arms):] + arms[:i % len(arms)]:
                argv = ["/usr/bin/time", "-o", str(timing), "-f", "%M", a.node, str(specs[arm][0] / "run_losat_wasi_threads.js"), str(specs[arm][1].resolve()),
                        *CASES[case], "-outfmt", "6", "-out", str(result), "-num_threads", str(a.threads)]
                result.unlink(missing_ok=True)
                t = time.monotonic()
                try:
                    r = subprocess.run(argv, env={**env0, "XDG_CACHE_HOME": str(a.out / f"cache-{arm}")}, capture_output=True, timeout=600)
                    code, err = r.returncode, r.stderr.decode(errors="replace")
                except subprocess.TimeoutExpired as e:
                    code, err = "timeout", (e.stderr or b"").decode(errors="replace")
                wall = time.monotonic() - t
                oob = "memory access out of bounds" in err
                sha = hashlib.sha256(result.read_bytes()).hexdigest()[:16] if result.exists() else "-"
                try: rss = int(timing.read_text().split()[-1])
                except (OSError, ValueError, IndexError): rss = 0
                s = stats[f"{case}/{arm}"]
                s["runs"] += 1; s["fail"] += code != 0; s["oob"] += oob
                if code == 0: s["shas"][sha] += 1; s["rss"].append(rss); s["wall"].append(round(wall, 2))
                if code != 0:
                    frames = re.findall(r"^\s+at \S*\.wasm\.(\S+) ", err, re.M)
                    s["frames"][frames[0][:70] if frames else "no wasm frame"] += 1
                    (a.out / f"fail-{case}-{arm}-{i}.txt").write_text(err[:6000])
                rows.write(f"{case}\t{arm}\t{i}\t{code}\t{wall:.3f}\t{rss}\t{int(oob)}\t{sha}\n"); rows.flush()
        (a.out / "summary.json").write_text(json.dumps(stats, indent=1, default=dict))
        print(f"round {i}", {k: (v["runs"], v["oob"], v["fail"]) for k, v in stats.items()}, flush=True)
    for k, v in sorted(stats.items()):
        print(k, "runs", v["runs"], "oob", v["oob"], "fail", v["fail"], "shas", dict(v["shas"]), "first frames", dict(v["frames"]),
              "median rss MiB", sorted(v["rss"])[len(v["rss"]) // 2] // 1024 if v["rss"] else None, "median wall", sorted(v["wall"])[len(v["wall"]) // 2] if v["wall"] else None)

main()
