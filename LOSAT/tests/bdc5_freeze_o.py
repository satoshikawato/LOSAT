#!/usr/bin/env python3
"""BDC5: freeze probe O's NCBI outputs of one platform into LOSAT/tests/blastdbcmd/eol/<name>/.
usage: freeze_o.py <probe out dir (holds manifest.tsv)> <eol dir of the platform>
Copies manifest.tsv and the stdout/stderr of every case whose exit status is not 255 (NCBI's
exception with its build path, which LOSAT refuses; the test reads nothing of those)."""
import shutil, sys
from pathlib import Path
src, dst = Path(sys.argv[1]), Path(sys.argv[2])
if dst.exists():
    shutil.rmtree(dst)
dst.mkdir(parents=True)
lines = (src / "manifest.tsv").read_text(encoding="utf-8").splitlines()
assert lines[-1] == "RUNDONE", "probe did not finish"
n = 0
for line in lines[1:-1]:
    cid, rc = line.split("\t")[:2]
    if rc == "255":
        continue
    for k in ("stdout", "stderr"):
        shutil.copyfile(src / f"{cid}.{k}", dst / f"{cid}.{k}")
    n += 1
shutil.copyfile(src / "manifest.tsv", dst / "manifest.tsv")
print(f"{dst}: {n} cases, {len(lines) - 2} in the manifest")
