#!/usr/bin/env python3
"""Build the databases of the `db.ll2.*` rows (line DBID item lineless2, 2026-11-15).

NCBI BLAST+ 2.17.0 makeblastdb is the comparison oracle only. The databases hold 14 subjects made from one
random query of 300 residues (the nucleotide ones from its codons) with 0, 2, 4, ... substitutions, so that
their ranks differ: r0 (the query itself), r1, the equal pair r2a and r2b, r3 to r12. The only line of r0, r2a
and r6 has a gi (version 4: `gi|N`) or a ti (version 5: `gnl|ti|N`; a gi-only line makes makeblastdb 2.17.0
abort writing the LMDB file): a `-negative_seqidlist` keeps those sequences and removes their lines, and NCBI
excludes them at the traceback's, the composition-based redo's and the ungapped re-evaluation's fetches
(`check_oid_exclusion`). `neg.txt` lists r12's accession (the OID mask removes it) and an id that no database
holds; `neggi.txt` lists a gi only (empty on version 5 once the gis are dropped: no id set, no line rule).
Queries: `q.faa`, `q.fna` (the query and its codons), `q2.faa`, `q2.fna` (the same and its first 100
residues, 300 bases). The databases are written with the creation date in their index files, so the frozen
outputs need the committed files, not a new run of this script.

`prot_v4bad` and `nucl_v4bad`: copies of the version 4 databases whose header of OID 8 (r7, which keeps its line)
starts with a wrong byte: under the negative list NCBI's fetch reads that header and the search fails although r7 is
not reported at `-max_target_seqs 3`; LOSAT refuses it (rows `db.ll2.*.bad.6`). `--damaged-only` writes these two from
the committed version 4 files.
`prot_v4u8`: a copy of `prot_v4` whose header of OID 1 (r0, lineless under `neg.txt`) holds its title as a UTF8String where
a VisibleString is expected: NCBI's warning `(803.10)` where it first reads that header (rows `db.ll2.bp.*u8*`);
`--damaged-only` writes it too.

usage: make_ll2.py --ncbi-bin DIR [--out LOSAT/tests/db_search/ll2] | make_ll2.py --damaged-only [--out DIR]"""
import os, random, shutil, struct, subprocess, sys


def damage(d):
    """`prot_v4bad`, `nucl_v4bad`: the version 4 databases of `d` with the first byte of OID 8's header changed."""
    for base, ext in (("prot_v4", "p"), ("nucl_v4", "n")):
        bad = base + "bad"
        for suffix in ("in", "hr", "sq", "si", "sd", "ni", "nd", "og", "js"):
            if os.path.exists(f"{d}/{base}.{ext}{suffix}"):
                shutil.copy(f"{d}/{base}.{ext}{suffix}", f"{d}/{bad}.{ext}{suffix}")
        pin = open(f"{d}/{bad}.{ext}in", "rb").read()
        at = 8
        at += 4 + struct.unpack(">I", pin[at:at + 4])[0]
        at += 4 + struct.unpack(">I", pin[at:at + 4])[0]
        n = struct.unpack(">I", pin[at:at + 4])[0]
        at += 4 + 8 + 4
        offsets = struct.unpack(f">{n + 1}I", pin[at:at + 4 * (n + 1)])
        phr = bytearray(open(f"{d}/{bad}.{ext}hr", "rb").read())
        assert phr[offsets[8]] == 0x30
        phr[offsets[8]] = 0x31
        open(f"{d}/{bad}.{ext}hr", "wb").write(bytes(phr))


def utf8_title(d):
    """`prot_v4u8`: a copy of `prot_v4` whose header of OID 1 (r0, its only line a gi) holds its title as a UTF8String
    (tag 0x0c) where a VisibleString (0x1a) is expected."""
    for suffix in ("in", "hr", "sq", "si", "sd", "ni", "nd", "og", "js"):
        if os.path.exists(f"{d}/prot_v4.p{suffix}"):
            shutil.copy(f"{d}/prot_v4.p{suffix}", f"{d}/prot_v4u8.p{suffix}")
    pin = open(f"{d}/prot_v4u8.pin", "rb").read()
    at = 8
    at += 4 + struct.unpack(">I", pin[at:at + 4])[0]
    at += 4 + struct.unpack(">I", pin[at:at + 4])[0]
    n = struct.unpack(">I", pin[at:at + 4])[0]
    at += 4 + 8 + 4
    offsets = struct.unpack(f">{n + 1}I", pin[at:at + 4 * (n + 1)])
    phr = bytearray(open(f"{d}/prot_v4u8.phr", "rb").read())
    title = phr.index(b"\xa0\x80\x1a", offsets[1], offsets[2]) + 2
    phr[title] = 0x0C
    open(f"{d}/prot_v4u8.phr", "wb").write(bytes(phr))


if sys.argv[1:2] == ["--damaged-only"]:
    out = sys.argv[3] if sys.argv[2:3] == ["--out"] else os.path.dirname(os.path.abspath(__file__))
    damage(out)
    utf8_title(out)
    sys.exit(0)
args = sys.argv[1:]
if len(args) < 2 or args[0] != "--ncbi-bin":
    sys.exit(__doc__)
NCBI = args[1]
D = args[3] if len(args) >= 4 and args[2] == "--out" else os.path.dirname(os.path.abspath(__file__))
rng = random.Random(20261115)
AA = "ACDEFGHIKLMNPQRSTVWY"
CODONS = {}
code = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
bases = "TCAG"
for i, a in enumerate(code):
    c = bases[i // 16] + bases[(i // 4) % 4] + bases[i % 4]
    CODONS.setdefault(a, []).append(c)
L = 300
q = [rng.choice(AA) for _ in range(L)]
qn = [rng.choice(CODONS[a]) for a in q]
def mutate(k, seed):
    r = random.Random(seed)
    p, n = list(q), list(qn)
    for pos in r.sample(range(L), k):
        a = r.choice([x for x in AA if x != p[pos]])
        p[pos] = a
        n[pos] = r.choice(CODONS[a])
    return "".join(p), "".join(n)
# (rank, substitutions, seed): r2a and r2b share the seed (identical sequences)
ranks = [("r0", 0, 1), ("r1", 2, 2), ("r2a", 4, 3), ("r2b", 4, 3)] + [(f"r{i}", 2 * i, 10 + i) for i in range(3, 13)]
seqs = {name: mutate(k, s) for name, k, s in ranks}
# ids per version: lineless (gi only / ti only) for r0, r2a, r6
def ids(name, i, mol, v5):
    # version 5: a gi-only line breaks makeblastdb 2.17.0 (the LMDB id list "not in ascending oid order"),
    # so the lineless sequences have a ti-only line there (`gnl|ti|N`, a general id without an accession)
    if name in ("r0", "r2a", "r6"):
        return f"gnl|ti|{7000 + i}" if v5 else f"gi|{5000 + i}"
    acc = ("AB" if mol == "n" else "BAA") + f"{100000 + i:06d}" + ".1"
    return f"dbj|{acc}|"
order = [n for n, _, _ in ranks]
# OIDs: mix the order so that lineless ones are not all first; r2a before r2b (ties go to the lower OID)
oid_order = ["r3", "r0", "r1", "r2a", "r2b", "r4", "r5", "r6", "r7", "r8", "r9", "r10", "r11", "r12"]
os.makedirs(D, exist_ok=True)
for mol in ("p", "n"):
    for v5 in (False, True):
        fa = f"{D}/{'prot' if mol == 'p' else 'nucl'}_{'v5' if v5 else 'v4'}.fa"
        with open(fa, "w") as f:
            for i, name in enumerate(oid_order):
                s = seqs[name][0 if mol == "p" else 1]
                f.write(f">{ids(name, i, mol, v5)} {name} subject\n{s}\n")
        out = fa[:-3]
        cmd = [f"{NCBI}/makeblastdb", "-in", fa, "-dbtype", "prot" if mol == "p" else "nucl", "-parse_seqids",
               "-out", out, "-title", os.path.basename(out), "-blastdb_version", "5" if v5 else "4"]
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
        os.remove(fa)
with open(f"{D}/q.faa", "w") as f:
    f.write(">qp query\n" + "".join(q) + "\n")
with open(f"{D}/q.fna", "w") as f:
    f.write(">qn query\n" + "".join(qn) + "\n")
# negative list: r12's accession (the OID mask drops it) and an id absent from the databases
r12 = oid_order.index("r12")
with open(f"{D}/neg.txt", "w") as f:
    f.write(f"BAA{100000 + r12:06d}.1\nAB{100000 + r12:06d}.1\nXX999999.1\n")
# a negative list of a gi only: empty for version 5 after the gis are dropped (no id set, no line rule)
with open(f"{D}/neggi.txt", "w") as f:
    f.write("gi|5001\n")
for ext in ("faa", "fna"):
    lines = open(f"{D}/q.{ext}").read().split("\n")
    n = 100 if ext == "faa" else 300
    with open(f"{D}/q2.{ext}", "w") as f:
        f.write(f">{lines[0][1:].split()[0]} query\n{lines[1]}\n>{lines[0][1:].split()[0]}2 query part\n{lines[1][:n]}\n")
damage(D)
utf8_title(D)
print("ok", oid_order)
