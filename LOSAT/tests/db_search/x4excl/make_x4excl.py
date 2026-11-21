#!/usr/bin/env python3
"""Build the databases of the `db.x4excl.*` rows (line RG item X4-excl, 2026-11-20).

NCBI BLAST+ 2.17.0 makeblastdb is the comparison oracle only. BLASTP `-use_sw_tback` with a negative seqid list: a
subject that the list keeps without a definition line (its only line has a ti, as in `ll2/prot_v5`) is excluded at the
fetch of the redo, after the early termination test of its match. The databases (version 5, `-parse_seqids`) hold
the failing pair of the `db.x4.sw.fail.*` rows (tks_W of `tests/fasta/rps/tks_w.faa` against sfail of
`tests/fasta/x4/sw_subjects.faa`) and such a subject:
- `fwd`: sexact (tks_W itself, one HSP: its redo succeeds), sfail, sok1 (`tests/fasta/x4/sw_subjects.faa`, the
  subject of qok1), and sx_frag (ti only: six letters of tks_W, MAGRVD, whose match comes after sfail's in the
  stream of tks_W, with a best e-value above 5 x 0.002 at `-searchsp 1000`).
- `rev`: sx_strong (ti only: tks_W itself, a match before sfail's), sfail.
`neg.txt` lists an id that no database holds; `neg_exact.txt` also lists sexact's accession (the OID mask removes
sexact). Queries: `tests/fasta/rps/tks_w.faa` and `tests/fasta/x4/sw_fail_last.faa`. The databases are written with
the creation date in their index files, so the frozen outputs need the committed files, not a new run of this script.
Version 5 because on version 4 a string list masks a sequence whose only line has a gi (NCBI reference (598d8ae6):
c++/src/objtools/blast/seqdb_reader/seqdboidlist.cpp:375-398): such a subject is never searched and never reaches the
redo (seen on ws with a gdb trace of NCBI's fetch, 2026-11-20).

usage: make_x4excl.py --ncbi-bin DIR [--out LOSAT/tests/db_search/x4excl]"""
import os, subprocess, sys

args = sys.argv[1:]
if len(args) < 2 or args[0] != "--ncbi-bin":
    sys.exit(__doc__)
NCBI = args[1]
D = args[3] if len(args) >= 4 and args[2] == "--out" else os.path.dirname(os.path.abspath(__file__))
TESTS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


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
sw = fasta(f"{TESTS}/fasta/x4/sw_subjects.faa")
assert "MAGRVE" in tks_w
dbs = {
    "fwd": [
        ("dbj|BAA200001.1|", "sexact tks_W itself", tks_w),
        ("dbj|BAA200002.1|", "sfail tks_M with D first, twice", sw["sfail"]),
        ("dbj|BAA200003.1|", "sok1 BDT63546.1 of PajaWSV", sw["sok1"]),
        ("gnl|ti|7004", "sx_frag six letters of tks_W", "MAGRVD"),
    ],
    "rev": [
        ("gnl|ti|7100", "sx_strong tks_W itself", tks_w),
        ("dbj|BAA200102.1|", "sfail tks_M with D first, twice", sw["sfail"]),
    ],
}
os.makedirs(D, exist_ok=True)
for name, records in dbs.items():
    fa = f"{D}/{name}.fa"
    with open(fa, "w") as f:
        for seqid, title, residues in records:
            f.write(f">{seqid} {title}\n{residues}\n")
    cmd = [f"{NCBI}/makeblastdb", "-in", fa, "-dbtype", "prot", "-parse_seqids", "-out", f"{D}/{name}",
           "-title", name, "-blastdb_version", "5"]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
    os.remove(fa)
with open(f"{D}/neg.txt", "w") as f:
    f.write("XX999999.1\n")
with open(f"{D}/neg_exact.txt", "w") as f:
    f.write("BAA200001.1\nXX999999.1\n")
print("ok", {name: [r[1].split()[0] for r in records] for name, records in dbs.items()})
