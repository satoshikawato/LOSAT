#!/usr/bin/env python3
"""Peak RSS and wall of the Node runner: A = dev host + raw module (guard-cache index hit), B = dev host + build-guarded
module, C = this branch's host + build-guarded module. Interleaved, order rotated per round; one line per run."""
import hashlib, json, os, subprocess, sys, tempfile
from pathlib import Path
B = Path.home() / "losat-baselines/wg-wasm-guard-20261101"
W = Path.home() / "losat-work/.worktrees/wg-guard/LOSAT/tests"
L = Path.home() / "losat-baselines"
T = B / "timing/data"
RAW = Path.home() / ".cache/losat-work/wt-wg-guard/threaded-command/wasm32-wasip1-threads/release/losat.wasm"
ART = B / "art/losat-threaded-command.wasm"
ARMS = {"A": (B / "mem/hostB", RAW), "B": (B / "mem/hostB", ART), "C": (W, ART), "N": (B / "mem/hostB", ART)}  # N = B again: the noise
sp4f, sp6 = T / "sp4_first3000.fa", L / "perf-of-20261013/ws/full7_a3/Species6.fa"
cdna, dna3m = L / "v4-A10/wasm_in/dmel_cdna_first31168.fa", L / "v4-A10/wasm_in/dmel_dna_3M.fa"
CASES = {
    "blastx": ["blastx", "-query", L / "v4-W/tmp/bxdata/g2l50k.fa", "-db", L / "v4-W/tmp/bxdata/dmel_pep", "-evalue", "1e-5", "-outfmt", "6"],
    "dc-megablast": ["blastn", "-task", "dc-megablast", "-query", cdna, "-subject", dna3m, "-outfmt", "6"],
    "blastp": ["blastp", "-query", sp4f, "-subject", sp6, "-outfmt", "6", "-evalue", "0.001"],
    "tblastx": ["tblastx", "-query", T / "txB_q.fa", "-subject", T / "txB_s.fa", "-outfmt", "6", "-query_gencode", "4", "-db_gencode", "4"],
}
rounds, threads, out = int(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
names = sys.argv[4].split(",") if len(sys.argv) > 4 else ["A", "B", "C"]
CASES = {c: CASES[c] for c in (sys.argv[5].split(",") if len(sys.argv) > 5 else CASES)}
env = {k: v for k, v in os.environ.items() if not k.startswith("LOSAT_") and k not in ("GLIBC_TUNABLES", "MALLOC_ARENA_MAX")}
env.update(NODE_NO_WARNINGS="1", TMPDIR=str(B / "tmp"), LOSAT_WASI_GUARD_CACHE=str(B / "mem/cache"))
with out.open("a") as log, tempfile.TemporaryDirectory(dir=B / "tmp") as d:
    def run(case, arm, r):
        host, module = ARMS[arm]
        res, txt = Path(d) / "out", Path(d) / "time"
        argv = ["/usr/bin/time", "-f", "%e %M %U %S", "-o", str(txt), "node", str(host / "run_losat_wasi_threads.js"), str(module),
                *map(str, CASES[case]), "-num_threads", threads, "-out", str(res)]
        p = subprocess.run(argv, env=env, capture_output=True, text=True)
        wall, rss, user, syst = txt.read_text().split()[-4:]
        sha = hashlib.sha256(res.read_bytes()).hexdigest() if res.exists() else None
        log.write(json.dumps(dict(case=case, arm=arm, round=r, status=p.returncode, wall_s=float(wall), rss_kib=int(rss), cpu_s=float(user) + float(syst), out_sha256=sha, load1=float(open("/proc/loadavg").read().split()[0]), err=p.stderr[-160:] if p.returncode else "")) + "\n"); log.flush()
    for case in CASES:
        for arm in names: run(case, arm, -1)   # discarded: the index of A, cold inputs
    for r in range(rounds):
        order = names[r % len(names):] + names[:r % len(names)]
        for case in CASES:
            for arm in order: run(case, arm, r)
