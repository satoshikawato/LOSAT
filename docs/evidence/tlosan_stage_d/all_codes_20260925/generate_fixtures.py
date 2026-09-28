#!/usr/bin/env python3
"""Generate code-changing TBLASTN local-subject inputs from pinned NCBI gc.prt."""
import hashlib
from pathlib import Path
import re
import subprocess

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
COMMIT = "598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4"
NCBI = Path("/mnt/c/Users/genom/GitHub/ncbi-blast")
assert subprocess.check_output(["git", "-C", str(NCBI), "rev-parse", "HEAD"], text=True).strip() == COMMIT
gc = subprocess.check_output(["git", "-C", str(NCBI), "show",
    f"{COMMIT}:c++/src/objects/seqfeat/gc.prt"], text=True)
tables = {int(code): aa for code, aa in re.findall(
    r'\bid\s+(\d+)\s*,\s*ncbieaa\s+"([A-Z*]+)"', gc)}
assert len(tables) == 27 and all(len(t) == 64 for t in tables.values())
codons = [a+b+c for a in "TCAG" for b in "TCAG" for c in "TCAG"]
stage_a = REPO / "docs/evidence/tlosan_stage_a/fixtures"
read = lambda path: "".join(line.strip() for line in path.read_text().splitlines()
                            if not line.startswith(">"))
base_query = read(stage_a / "query.faa")
base_subject = read(stage_a / "subject_code1.fna")
assert len(base_query) == 120 and len(base_subject) == 360
positions = (15,35,55,75,95)
rows = ["code\tcodon\tcode1_aa\tselected_aa\tpositions\tnote"]
for code, table in tables.items():
    changed = [(codons[i], tables[1][i], table[i]) for i in range(64)
               if table[i] != tables[1][i] and table[i] != "*"]
    codon, standard, selected = changed[0] if changed else ("TGG", "W", "W")
    q = list(base_query)
    s = list(base_subject)
    for pos in positions:
        q[pos] = selected
        s[3*pos:3*pos+3] = codon
    query = ">q1\n" + "".join(q) + "\n"
    subject = ">s_code%d\n" % code + "".join(s) + "\n"
    (HERE / "fixtures" / f"code{code}.faa").write_text(query)
    (HERE / "fixtures" / f"code{code}.fna").write_text(subject)
    note = "changed amino acid" if changed else "identical internal translation"
    rows.append(f"{code}\t{codon}\t{standard}\t{selected}\t{','.join(map(str,positions))}\t{note}")
(HERE / "fixtures.tsv").write_text("\n".join(rows) + "\n")
files = sorted((HERE / "fixtures").glob("*"))
(HERE / "fixtures.sha256").write_text("".join(
    f"{hashlib.sha256(path.read_bytes()).hexdigest()}  fixtures/{path.name}\n" for path in files))
(HERE / "generation.manifest.txt").write_text(
    f"NCBI gc.prt source: {COMMIT}\n"
    f"NCBI gc.prt SHA256: {hashlib.sha256(gc.encode()).hexdigest()}\n"
    f"Source query SHA256: {hashlib.sha256((stage_a/'query.faa').read_bytes()).hexdigest()}\n"
    f"Source subject SHA256: {hashlib.sha256((stage_a/'subject_code1.fna').read_bytes()).hexdigest()}\n"
    "Code 11 and 23 have unchanged ncbieaa translation; TGG checks their identical translated path.\n"
)
