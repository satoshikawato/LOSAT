#!/usr/bin/env python3
"""NCBI's report (stderr) on a damaged nest against LOSAT's, line by line (line NP M (b), 2026-11-23).

usage: compare_reports.py NCBI_ERR LOSAT_ERR
NCBI writes "Error: NCBI C++ Exception:" and one `T0 "<file>", line <n>: Error: (<class>) <function> - <text>` line
per rethrow; LOSAT (crate::wording) writes "Error: <text>" and then each frame's text on a line of its own. The texts
(the message and every frame) must be the same, in the same order.
"""
import re
import sys

n = open(sys.argv[1], encoding="latin-1").read().split("\n")
lo = open(sys.argv[2], encoding="latin-1").read().split("\n")
assert n[0] == "Error: NCBI C++ Exception:", n[0]
pat = re.compile(r'^\s+T0 ".*", line \d+: Error: \(\S+\) \S+ - (.*)$')
nm = []
for x in n[1:]:
    if x == "":
        nm.append("")
        continue
    m = pat.match(x)
    assert m, x[:200]
    nm.append(m.group(1).strip())
assert lo[0].startswith("Error: "), lo[0][:80]
lm = [lo[0][len("Error: "):].strip()] + [x.strip() for x in lo[1:]]
print(f"NCBI lines {len(n)}, LOSAT lines {len(lo)}, texts compared {len(nm)}")
if nm == lm:
    print("REPORT TEXTS EQUAL")
else:
    for i, (a, b) in enumerate(zip(nm, lm)):
        if a != b:
            print("REPORT TEXTS DIFFER at", i, repr(a[:160]), repr(b[:160]))
            break
    else:
        print("REPORT TEXTS DIFFER in length", len(nm), len(lm))
