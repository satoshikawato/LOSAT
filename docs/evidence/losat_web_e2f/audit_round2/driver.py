#!/usr/bin/env python3
import subprocess, sys, json, os, importlib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
HERE = Path(__file__).resolve().parent
NCBI = "/home/kawato/micromamba/bin/blastn"
OLD = "/home/kawato/.cache/losat-web-gui-target/audit_s07pp_round2/old_target/release/LOSAT"
LOSAT = "/home/kawato/.cache/losat-web-gui-target/s07pg-native/release/LOSAT"
WARN = b"Warning: [blastn] 'num_threads' is currently ignored when 'subject' is specified.\n"
sys.path.insert(0, str(HERE))
import cases as C

def run(argv, env=None):
    p = subprocess.run(argv, cwd=HERE, capture_output=True, env=env)
    return p.returncode, p.stdout, p.stderr

def strip_opt(args, key):
    out, i = [], 0
    while i < len(args):
        if args[i] == key: i += 2
        else: out.append(args[i]); i += 1
    return out

def one(c):
    a = ["-query", c["query"], "-subject", c["subject"], *c["args"]]
    if "-outfmt" not in c["args"]: a += ["-outfmt", "6"]
    n = run([NCBI, *a]); l = run([LOSAT, "blastn", *a])
    n = (n[0], n[1], n[2].replace(WARN, b""))
    same = n == l
    o = run([OLD, "blastn", *a]); o = o
    old_same = n == o
    rejected = (not same) and l[0] == 2 and l[2].startswith(b"error: the NCBI BLAST+ option")
    # discriminator metrics from NCBI
    ua = ["-query", c["query"], "-subject", c["subject"], *strip_opt(strip_opt(c["args"], "-max_target_seqs"), "-outfmt"), "-outfmt", "6", "-max_target_seqs", "1000000"]
    u = run([NCBI, *ua])
    hit_subj = len({ln.split(b"\t")[1] for ln in u[1].splitlines()}) if u[0] == 0 else -1
    fmt = "6"
    if "-outfmt" in c["args"]: fmt = c["args"][c["args"].index("-outfmt") + 1]
    out_subj = len({ln.split(b"\t")[1] for ln in n[1].splitlines() if ln and not ln.startswith(b"#")}) if fmt in ("6", "7") else -1
    if not same and not rejected:
        d = HERE / "diff" / c["name"]; d.mkdir(parents=True, exist_ok=True)
        (d / "ncbi.out").write_bytes(n[1]); (d / "ncbi.err").write_bytes(n[2]); (d / "ncbi.rc").write_text(str(n[0]))
        (d / "losat.out").write_bytes(l[1]); (d / "losat.err").write_bytes(l[2]); (d / "losat.rc").write_text(str(l[0]))
        (d / "argv").write_text(" ".join(a))
    return c["name"], c["cat"], ("rejected" if rejected else same), n[0], l[0], len(n[1]), hit_subj, out_subj, " ".join(a), old_same

if __name__ == "__main__":
    only = sys.argv[1] if len(sys.argv) > 1 else ""
    sel = [c for c in C.cases if c["name"].startswith(only)]
    with ThreadPoolExecutor(6) as ex: res = list(ex.map(one, sel))
    out = sys.argv[2] if len(sys.argv) > 2 else "results.tsv"
    with open(out, "w") as f:
        f.write("name\tcat\tresult\tncbi_rc\tlosat_rc\tncbi_bytes\tuncapped_hit_subjects\toutput_subjects\tprefix_binary_matches_ncbi\targv\n")
        for r in res: f.write(f"{r[0]}\t{r[1]}\t{'rejected' if r[2]=='rejected' else ('same' if r[2] else 'DIFFERENT')}\t{r[3]}\t{r[4]}\t{r[5]}\t{r[6]}\t{r[7]}\t{'yes' if r[9] else 'NO'}\t{r[8]}\n")
    print(len(res), "cases;", sum(r[2] is True for r in res), "identical;", sum(r[2] is False for r in res), "different;", sum(r[2]=="rejected" for r in res), "rejected;", sum(not r[9] for r in res), "prefix-binary differs")
