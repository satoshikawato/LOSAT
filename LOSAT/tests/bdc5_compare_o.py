#!/usr/bin/env python3
"""BDC5: compare LOSAT's probe O outputs with NCBI's of the same platform.
usage: bdc5_compare_o.py --ncbi DIR --losat DIR [--crlf]
--crlf: NCBI's outputs are Windows text-mode bytes; each '\\r' before a '\\n' is taken away
(LOSAT writes '\\n' on every platform). NCBI's exit status 255 (its exception that names its
build path) must be LOSAT's refusal: exit 1, empty stdout, 'LOSAT's BLASTDBCMD' in stderr.
Prints one row per case (match / refuse / FAIL) and a total; exit 1 on any FAIL."""
import argparse, sys
from pathlib import Path

def manifest(d):
    rows = {}
    lines = (d / "manifest.tsv").read_text(encoding="utf-8").splitlines()
    assert lines[-1] == "RUNDONE", f"{d}: probe did not finish"
    for line in lines[1:-1]:
        f = line.split("\t")
        rows[f[0]] = f
    return rows

def uncrlf(b):
    out = bytearray()
    for i, c in enumerate(b):
        if not (c == 13 and i + 1 < len(b) and b[i + 1] == 10):
            out.append(c)
    return bytes(out)

p = argparse.ArgumentParser()
p.add_argument("--ncbi", type=Path, required=True)
p.add_argument("--losat", type=Path, required=True)
p.add_argument("--crlf", action="store_true")
a = p.parse_args()
n, l = manifest(a.ncbi), manifest(a.losat)
assert n.keys() == l.keys(), (sorted(n.keys() ^ l.keys()))
count = {"match": 0, "refuse": 0, "FAIL": 0}
for cid in n:
    so = (a.losat / f"{cid}.stdout").read_bytes()
    se = (a.losat / f"{cid}.stderr").read_bytes()
    if n[cid][1] == "255":
        ok = l[cid][1] == "1" and so == b"" and b"LOSAT's BLASTDBCMD" in se
        kind = "refuse" if ok else "FAIL"
    else:
        no = (a.ncbi / f"{cid}.stdout").read_bytes()
        ne = (a.ncbi / f"{cid}.stderr").read_bytes()
        if a.crlf:
            no, ne = uncrlf(no), uncrlf(ne)
        ok = l[cid][1] == n[cid][1] and so == no and se == ne
        kind = "match" if ok else "FAIL"
    count[kind] += 1
    print(f"{cid}\t{kind}\tncbi_rc={n[cid][1]}\tlosat_rc={l[cid][1]}")
print("total " + " ".join(f"{k} {v}" for k, v in count.items()))
sys.exit(1 if count["FAIL"] else 0)
