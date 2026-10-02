"""D3/D4 (E2g T9/T14): merged stdout+stderr order and outfmt 0 write failure, NCBI vs LOSAT."""
import subprocess, sys, json, os, itertools
from pathlib import Path
NCBI = "/home/kawato/micromamba/bin/blastn"; LOSAT = sys.argv[1]
W = Path("/home/kawato/.cache/losat-web-gui-target/e2g-audit/b/w")
env = {k: v for k, v in os.environ.items() if k not in ("BATCH_SIZE", "CHUNK_SIZE", "OVERLAP_CHUNK_SIZE")}
inputs = sorted({p.name[:-5] for p in W.glob("T89p*_q.fa")})
rows = []
def run(binary, args, mode, extra_env=None):
    argv = [binary] + ([] if binary == NCBI else ["blastn"]) + args
    e = dict(env, **(extra_env or {}))
    if mode == "merged":
        p = subprocess.run(argv, cwd=W, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=e)
        return p.returncode, p.stdout, b""
    if mode == "devfull":
        with open("/dev/full", "wb") as f:
            p = subprocess.run(argv, cwd=W, stdout=f, stderr=subprocess.PIPE, env=e)
        return p.returncode, b"", p.stderr
    if mode == "out_devfull":
        p = subprocess.run(argv + ["-out", "/dev/full"], cwd=W, capture_output=True, env=e)
        return p.returncode, p.stdout, p.stderr
    p = subprocess.run(argv, cwd=W, capture_output=True, env=e)
    return p.returncode, p.stdout, p.stderr
for name, fmt, mode, batch in itertools.product(inputs, ("0", "6", "7"), ("merged", "separate", "devfull", "out_devfull"), ("", "1", "200")):
    if mode in ("devfull", "out_devfull") and fmt != "0":
        continue
    args = ["-query", f"{name}_q.fa", "-subject", f"{name}_s.fa", "-outfmt", fmt]
    extra = {"BATCH_SIZE": batch} if batch else {}
    n = run(NCBI, args, mode, extra); l = run(LOSAT, args, mode, extra)
    rows.append(dict(input=name, fmt=fmt, mode=mode, batch=batch, same=n == l, ncbi_exit=n[0], losat_exit=l[0]))
    if n != l:
        print("DIFF", json.dumps(rows[-1]), flush=True)
json.dump(rows, open(Path(__file__).with_name("results.json"), "w"), indent=1)
print(len(rows), "runs", sum(r["same"] for r in rows), "same")
