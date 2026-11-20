"""PROBE (not for merge): A/B of two pylosat wheels on one runner, runs alternating.

Usage: ab_probe.py <python A> <python B> <fasta dir> <out.tsv>
Each run is a fresh process of the arm's venv Python doing one workload; its CPU time
(user + system) and max RSS come from os.wait4. Reports must be the same in both arms.
"""
import hashlib, os, statistics, subprocess, sys, time

WL = r'''
import sys, hashlib, pylosat
W, reps, D = sys.argv[1], int(sys.argv[2]), sys.argv[3]
if W == "W1":
    s = pylosat.Search("tblastx", num_threads=1); q, sj = "AP027131.fasta", "NZ_CP006932.fasta"
elif W == "W2":
    s = pylosat.Search("blastn", task="blastn", num_threads=1); q, sj = "AP027131.fasta", "NZ_CP006932.fasta"
else:
    s = pylosat.Search("blastp", num_threads=1); q, sj = "SicyWSV.faa", "PajaWSV.faa"
q = pylosat.Sequences.from_fasta(D + "/" + q); sj = pylosat.Sequences.from_fasta(D + "/" + sj)
hs = set()
for _ in range(reps):
    b = s.run(q, sj, outfmt="6").reports[0].bytes
    assert b
    hs.add(hashlib.sha256(b).hexdigest())
assert len(hs) == 1
print(hs.pop())
'''
REPS = {"W1": 1, "W2": 5, "W3": 25}
ROUNDS = 7


def run(py, w, d):
    t0 = time.perf_counter()
    p = subprocess.Popen([py, "-I", "-c", WL, w, str(REPS[w]), d], stdout=subprocess.PIPE)
    out = p.stdout.read().decode().strip()
    _, status, ru = os.wait4(p.pid, 0)
    wall = time.perf_counter() - t0
    if status != 0:
        raise SystemExit(f"{py} {w}: status {status}")
    return ru.ru_utime + ru.ru_stime, ru.ru_maxrss, wall, out


def main():
    pa, pb, d, out = sys.argv[1:5]
    rows = []
    for w in ("W1", "W2", "W3"):
        for r in range(1, ROUNDS + 1):
            for arm, py in (("gnu", pa), ("zig", pb)) if r % 2 else (("zig", pb), ("gnu", pa)):
                cpu, rss, wall, h = run(py, w, d)
                rows.append((w, arm, r, cpu, rss, wall, h))
                print(w, arm, r, f"{cpu:.3f}", rss, f"{wall:.3f}", h[:12], flush=True)
    with open(out, "w") as f:
        f.write("workload\tarm\tround\tcpu_s\tmax_rss_kb\twall_s\tsha256\n")
        for row in rows:
            f.write("\t".join(str(x) for x in row) + "\n")
    print("| W | cpu gnu (med) | cpu zig (med) | zig/gnu | rss gnu | rss zig | zig/gnu | same report |")
    print("|---|---|---|---|---|---|---|---|")
    for w in ("W1", "W2", "W3"):
        g = [x for x in rows if x[0] == w and x[1] == "gnu"]; z = [x for x in rows if x[0] == w and x[1] == "zig"]
        cg = statistics.median(x[3] for x in g); cz = statistics.median(x[3] for x in z)
        rg = statistics.median(x[4] for x in g); rz = statistics.median(x[4] for x in z)
        same = len({x[6] for x in g + z}) == 1
        print(f"| {w} | {cg:.3f} | {cz:.3f} | {cz/cg:.4f} | {rg} | {rz} | {rz/rg:.4f} | {same} |")


main()
