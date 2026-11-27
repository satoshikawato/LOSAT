#!/usr/bin/env python3
"""Write the inputs of the `db.sw.*` rows of line PS-SW's fix round 1 (PSIBLAST `-use_sw_tback`, 2026-11-27).

NCBI BLAST+ 2.17.0 makeblastdb is the comparison oracle only. From tks_W (`tests/fasta/rps/tks_w.faa`), tks_M
(`tests/fasta/rps/tks_m.faa`) and sfail (`tests/fasta/x4/sw_subjects.faa`: tks_M with D first, twice), and the tk PSSM
`tests/fixtures/pssm/ps4_tk_paja.asn.gz`, it writes:
- `tests/fasta/pssw/sfail.faa`: sfail alone, a subject whose Smith-Waterman redo fails (two HSPs, and the query's first
  letter scores D at 0 or less), so that a search whose last redone match is sfail's ends in `Traceback failed`.
- `tests/fasta/pssw/smm.faa`: tks_M twice, a subject whose redo never ends in NCBI (two HSPs, M against the query's first
  letter above 0).
- the version 5 databases (`-parse_seqids`) `tests/db_search/pssw/psw` (sx_frag, whose only line has a ti, then sfail)
  and `pswk` (sexact, tks_W itself, then sx_frag and sfail), searched with `tests/db_search/x4excl/neg.txt`: the negative
  list keeps sx_frag without a line, and NCBI's fetch of the redo excludes it, which sets the redo's status back to 0.
  The redo runs the matches in OID order at `-comp_based_stats 0` and in e-value order at 1 to 3, so the excluded
  subject comes before sfail's failing match in the one and after it in the other.
- `tests/fixtures/pssm/ps_sw_tk_scores.asn.gz`: the tk PSSM without the Pssm member `intermediateData` (scores only, no
  frequency ratios), gzipped without a name or time.
The databases are written with the creation date in their index files, so the frozen outputs need the committed files,
not a new run of this script.

usage: make_pssw.py --ncbi-bin DIR"""
import gzip
import os
import subprocess
import sys

args = sys.argv[1:]
if len(args) != 2 or args[0] != "--ncbi-bin":
    sys.exit(__doc__)
NCBI = args[1]
D = os.path.dirname(os.path.abspath(__file__))
TESTS = os.path.dirname(os.path.dirname(D))


def fasta(path):
    """The records of `path` as {first word of the title: residues}."""
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = ""
        elif name:
            seqs[name] += line
    return seqs


tks_w = fasta(f"{TESTS}/fasta/rps/tks_w.faa")["tks_W"]
tks_m = fasta(f"{TESTS}/fasta/rps/tks_m.faa")["tks_M"]
sfail = fasta(f"{TESTS}/fasta/x4/sw_subjects.faa")["sfail"]
assert sfail == "D" + tks_m[1:] + "D" + tks_m[1:], "sfail is tks_M with D first, twice"
assert "MAGRVE" in tks_w

os.makedirs(f"{TESTS}/fasta/pssw", exist_ok=True)
with open(f"{TESTS}/fasta/pssw/sfail.faa", "w") as f:
    f.write(f">sfail tks_M with D first, twice\n{sfail}\n")
with open(f"{TESTS}/fasta/pssw/smm.faa", "w") as f:
    f.write(f">smm tks_M twice\n{tks_m}{tks_m}\n")

dbs = {
    "psw": [
        ("gnl|ti|7204", "sx_frag six letters of tks_W", "MAGRVD"),
        ("dbj|BAA200302.1|", "sfail tks_M with D first, twice", sfail),
    ],
    "pswk": [
        ("dbj|BAA200401.1|", "sexact tks_W itself", tks_w),
        ("gnl|ti|7304", "sx_frag six letters of tks_W", "MAGRVD"),
        ("dbj|BAA200403.1|", "sfail tks_M with D first, twice", sfail),
    ],
}
for name, records in dbs.items():
    fa = f"{D}/{name}.fa"
    with open(fa, "w") as f:
        f.write("".join(f">{i} {t}\n{r}\n" for i, t, r in records))
    subprocess.run([f"{NCBI}/makeblastdb", "-in", fa, "-dbtype", "prot", "-parse_seqids", "-out", f"{D}/{name}",
                    "-title", name, "-blastdb_version", "5"], check=True, stdout=subprocess.DEVNULL)
    os.remove(fa)

tk = gzip.open(f"{TESTS}/fixtures/pssm/ps4_tk_paja.asn.gz", "rt").read()
start, end = tk.index("    intermediateData {\n"), tk.index("    finalData {\n")
with open(f"{TESTS}/fixtures/pssm/ps_sw_tk_scores.asn.gz", "wb") as raw:
    with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as f:
        f.write((tk[:start] + tk[end:]).encode())
