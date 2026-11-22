#!/usr/bin/env python3
"""Build the databases of the `db.x4b.*` rows that search a database (line RG item X4b-1 fix rounds 1 and 2,
2026-11-20 and 21).

NCBI BLAST+ 2.17.0 makeblastdb is the comparison oracle only. The subjects are X4's crafted proteins back-translated
with one codon per amino acid (`tests/fasta/x4b/ok1_fail_rc.fna`): sok1 in frame +1 and sfail in frame -1, which
make tks_W's Smith-Waterman redo fail (`-use_sw_tback`, PD-LOSAT-NCBI-DEFECTS 1.72).

- `fail`: sok1_nt (OID 0) and sfail_rc (OID 1), each with an accession line (`dbj|AB30000N.1|`): the `-db` form of
  the `-subject` failure rows (`db.x4b.dbfail.*`, `db.x4b.cut0.db.*`, `db.x4b.sp.db.fail.*`, `db.x4b.sp.two.db.*`).
- `excl`: the same two and, at OID 2, `sx_rc`: sfail's first half with 8 substitutions, in frame -1, whose only line
  has a ti (`gnl|ti|9003`, as `tests/db_search/ll2`'s lineless sequences). tks_W's preliminary match on it comes
  after the one on sfail_rc in both of NCBI's stream orders (a worse best e-value with composition-based statistics,
  a higher OID without them). `neg.txt` lists only an id that no database holds: a negative seqid list with ids
  keeps sx_rc and removes its line, and NCBI's fetch with `check_oid_exclusion` excludes it at the redo
  (`seqsrc_seqdb.cpp:306-314`), after its early termination test; NCBI then sets the batch's status back to 0
  (`blast_kappa.c:3560-3575`), which clears the failure of the match before (rows `db.x4b.excl.*`,
  `db.x4b.sp.db.excl.*`).
- `polycg`: `tests/fasta/x4b/polycg.fna` (repeats of CCCCGGGGAA, no parsed Seq-ids), on which a PSSM finds no hit
  (rows `db.x4b.pssm.none.db.*`: `-in_pssm -use_sw_tback` ends as NCBI's).

The databases are written with the creation date in their index files, so the frozen outputs need the committed
files, not a new run of this script.

usage: make_x4b.py --ncbi-bin DIR [--out LOSAT/tests/db_search/x4b]"""
import os, subprocess, sys

CODON = dict(A="GCT", R="CGT", N="AAT", D="GAT", C="TGT", Q="CAA", E="GAA", G="GGT", H="CAT", I="ATT",
             L="CTG", K="AAA", M="ATG", F="TTT", P="CCG", S="TCT", T="ACT", W="TGG", Y="TAT", V="GTT")
COMP = str.maketrans("ACGT", "TGCA")
BODY = "QLILSHHLTMAGRVELVTGPMFAGKSTYLKNIYQQENGGNKHCLFVKHSLETRYGCGTGTIVTHAGEVIEGCTTVSSIK"
SOK1 = ("MATAPAALPFGKANAAACGNSNNTFDVLEGLVKRARVVLVFITALSIGLIVVAIMNLRMGMGNSAQGVPLLDEHGQPVLDGNDDPQIEEPNNTYPNKTRDYLI"
        "LMASILSFVIGAASSYTTRLK")
SFAIL = "D" + BODY + "D" + BODY


def nt(protein):
    return "".join(CODON[a] for a in protein)


def rc(seq):
    return seq.translate(COMP)[::-1]


def weaker(protein):
    """sfail's first half with 8 substitutions (every tenth letter from the 6th becomes A, or G where it is A)."""
    letters = list(protein)
    for at in range(5, len(letters), 10):
        letters[at] = "G" if letters[at] == "A" else "A"
    return "".join(letters)


def fasta(records):
    out = []
    for name, seq in records:
        out.append(">" + name)
        out += [seq[i:i + 70] for i in range(0, len(seq), 70)]
    return "\n".join(out) + "\n"


args = sys.argv[1:]
if len(args) < 2 or args[0] != "--ncbi-bin":
    sys.exit(__doc__)
NCBI = args[1]
D = args[3] if len(args) >= 4 and args[2] == "--out" else os.path.dirname(os.path.abspath(__file__))
SX = weaker("D" + BODY)
base = [("dbj|AB300001.1| sok1_nt sok1 in frame +1", nt(SOK1)), ("dbj|AB300002.1| sfail_rc sfail in frame -1", rc(nt(SFAIL)))]
dbs = {"fail": base, "excl": base + [("gnl|ti|9003 sx_rc sfail's first half with 8 substitutions in frame -1", rc(nt(SX)))]}
os.makedirs(D, exist_ok=True)
for name, records in dbs.items():
    fa = f"{D}/{name}.fna"
    with open(fa, "w") as f:
        f.write(fasta(records))
    subprocess.run([f"{NCBI}/makeblastdb", "-in", fa, "-dbtype", "nucl", "-parse_seqids", "-out", f"{D}/{name}",
                    "-title", name, "-blastdb_version", "5"], check=True, stdout=subprocess.DEVNULL)
    os.remove(fa)
polycg = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "fasta", "x4b", "polycg.fna")
subprocess.run([f"{NCBI}/makeblastdb", "-in", polycg, "-dbtype", "nucl", "-out", f"{D}/polycg", "-title", "polycg",
                "-blastdb_version", "5"], check=True, stdout=subprocess.DEVNULL)
with open(f"{D}/neg.txt", "w") as f:
    f.write("XX999999.1\n")
print("ok", {name: len(records) for name, records in dbs.items()}, SX)
