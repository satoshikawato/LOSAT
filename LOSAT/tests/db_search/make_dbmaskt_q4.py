#!/usr/bin/env python3
"""Series NP, M4b-3c (2026-11-23): the inputs of the rows db.tmask.*.hard.te.{13,14,15,31}.*, *.tq4.* and their near
rows (PD-LOSAT-NCBI-DEFECTS exception 27) in dbmaskt/ (README.md there), from dbmaskt/edge.fna and dbmaskt/qxe.fna.

A mask that leaves a reading frame an empty first part (within the first two bases of a sequence, within the one or
two bases after its last whole codon, or an empty mask at 0) makes NCBI BLAST+ 2.17.0 crash in every run. The near
input of each such mask is the same database with that mask dropped and its bases stored as N (the Owner's answer
Q6 = B of 2026-11-23: the bases stay masked in every frame; the literal mask removal, with the real bases kept, is
reading A and gives other rows), as in NCBI's own runs where another mask avoids the crash
(ncbi-parity-20261019/M/m4b3c/Q4B_PROGRESS.md). It equals LOSAT's reading only as long as NCBI's preliminary stage
reads no masked residue: makeblastdb stores an N as a random base, which differs from the real one at every masked
position here, and the near input's seed scan and early extensions read that base where LOSAT reads the real one.
Writes:
- edge_n13.fna, edge_n14.fna, edge_n15.fna, edge_n31.fna: edge.fna with the bases of te's masks 13 (e2 [1200, 1202)),
  14 (e3 [0, 1)), 15 (e3 [1, 2)) and 31's [1200, 1201) of e4 stored as N; edge_nq4.fna: e3's bases 1196-1199 as N;
- hedge4_<name>.asnt (algorithm 11 of each database): n13, n14, n15 the last 3 bases of e5 (10 bases, a sequence
  that no search reads: no mask is left on the searched sequences), n31 e4 [0, 3), q4 e3's empty mask at 0 (`from 0,
  to 4294967295`) and [1100, 1200), nq4 e3 [1100, 1196) (the tail 4 bases shorter: the near input of exception 21);
- qneq4.faa: TBLASTN queries, 60 residues of each frame end that holds a masked base (the frames that NCBI searches
  whole and the other strand's frames whose last residue holds the base);
- qxeq4.fna: the TBLASTX query x4t, e4's bases 1001-1200 (qxe.fna has e2's end and the starts of e3 and e4).
With NCBI_BIN, builds te13n, te14n, te15n, te31n, tq4 (edge.fna with hedge4_q4.asnt) and tq4n (version 4) with that
directory's makeblastdb.

usage: make_dbmaskt_q4.py DBMASKT_DIR [NCBI_BIN]
"""
import os
import subprocess
import sys
from pathlib import Path

D = Path(sys.argv[1])
NCBI_BIN = sys.argv[2] if len(sys.argv) > 2 else None
COMP = str.maketrans("ACGT", "TGCA")
CODE = {}
_B = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, a in enumerate(_B):
    for j, b in enumerate(_B):
        for k, c in enumerate(_B):
            CODE[a + b + c] = _AAS[16 * i + 4 * j + k]


def read_fasta(path):
    records, name, seq = [], None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(seq)))
            name, seq = line[1:], []
        else:
            seq.append(line.strip())
    records.append((name, "".join(seq)))
    return records


def fasta(records):
    out = []
    for name, seq in records:
        out.append(f">{name}\n")
        out.extend(seq[i:i + 60] + "\n" for i in range(0, len(seq), 60))
    return "".join(out)


def frame_letters(seq, frame):
    """Frame `frame` (+1..+3 from base frame-1, -1..-3 from the minus strand's base |frame|-1), standard code."""
    strand = seq if frame > 0 else seq.translate(COMP)[::-1]
    strand = strand[abs(frame) - 1:]
    return "".join(CODE[strand[i:i + 3]] for i in range(0, len(strand) - 2, 3))


def mask_file(name, options, masks):
    """masks: [(local id, start, end)], [start, end); the empty [0, 0) as `from 0, to 4294967295`."""
    items = []
    for local, start, end in masks:
        to = 4294967295 if (start, end) == (0, 0) else end - 1
        items.append(f"        {{\n          from {start},\n          to {to},\n          id local id {local}\n        }}")
    text = (f"Blast-db-mask-info ::= {{\n  algo-id 2,\n  algo-program dust,\n  algo-options \"{options}\",\n"
            "  masks {\n    masks {\n      packed-int {\n" + ",\n".join(items) +
            "\n      }\n    },\n    more FALSE\n  }\n}\n")
    (D / f"hedge4_{name}.asnt").write_text(text)


edge = read_fasta(D / "edge.fna")
seqs = {name.split()[0]: seq for name, seq in edge}
assert (len(seqs["e2"]), len(seqs["e3"]), len(seqs["e4"]), len(seqs["e5"])) == (1202, 1200, 1201, 10)

NEAR = {  # file: (sequence, bases stored as N)
    "n13": ("e2", (1200, 1201)),
    "n14": ("e3", (0,)),
    "n15": ("e3", (1,)),
    "n31": ("e4", (1200,)),
    "nq4": ("e3", (1196, 1197, 1198, 1199)),
}
for name, (target, positions) in NEAR.items():
    records = []
    for title, seq in edge:
        if title.split()[0] == target:
            letters = list(seq)
            for p in positions:
                letters[p] = "N"
            seq = "".join(letters)
        records.append((title, seq))
    (D / f"edge_{name}.fna").write_text(fasta(records))

mask_file("n13", "near te 13: e2 [1200, 1202) removed, stored as N", [(5, 7, 10)])
mask_file("n14", "near te 14: e3 [0, 1) removed, stored as N", [(5, 7, 10)])
mask_file("n15", "near te 15: e3 [1, 2) removed, stored as N", [(5, 7, 10)])
mask_file("n31", "near te 31: e4 [1200, 1201) removed, stored as N", [(4, 0, 3)])
mask_file("q4", "e3 empty mask at 0 and tail", [(3, 0, 0), (3, 1100, 1200)])
mask_file("nq4", "near q4: e3 [1100, 1196), 1196-1199 stored as N", [(3, 1100, 1196)])

QUERIES = [  # name, sequence, frame, the frame's start (else its end)
    ("t3p1h", "e3", 1, True),
    ("t3p2h", "e3", 2, True),
    ("t3m1t", "e3", -1, False),
    ("t3m3t", "e3", -3, False),
    ("t2m1h", "e2", -1, True),
    ("t2m2h", "e2", -2, True),
    ("t2p2t", "e2", 2, False),
    ("t2p3t", "e2", 3, False),
    ("t4m1h", "e4", -1, True),
    ("t4p2t", "e4", 2, False),
    ("t4m1t", "e4", -1, False),
]
proteins = []
for name, target, frame, head in QUERIES:
    letters = frame_letters(seqs[target], frame)
    if head:
        title, part = f"{name} {target} frame {frame:+d} residues 0-59", letters[:60]
    else:
        first = len(letters) - 60
        title, part = f"{name} {target} frame {frame:+d} residues {first}-{len(letters) - 1}", letters[first:]
    proteins.append((title, part))
(D / "qneq4.faa").write_text(fasta(proteins))
(D / "qxeq4.fna").write_text(fasta([("x4t e4 1001-1200", seqs["e4"][1001:1201])]))

if NCBI_BIN:
    env = dict(os.environ, BLAST_USAGE_REPORT="false")
    builds = [
        ("te13n", "edge_n13.fna", "n13", "near te 13"),
        ("te14n", "edge_n14.fna", "n14", "near te 14"),
        ("te15n", "edge_n15.fna", "n15", "near te 15"),
        ("te31n", "edge_n31.fna", "n31", "near te 31"),
        ("tq4", "edge.fna", "q4", "empty mask at 0 and a tail"),
        ("tq4n", "edge_nq4.fna", "nq4", "near tq4"),
    ]
    for db, source, mask, title in builds:
        subprocess.run([f"{NCBI_BIN}/makeblastdb", "-in", source, "-dbtype", "nucl", "-out", db, "-title", title,
                        "-blastdb_version", "4", "-mask_data", f"hedge4_{mask}.asnt"], cwd=D, env=env, check=True,
                       stdout=subprocess.DEVNULL)
    for stale in D.glob("*.njs"):
        stale.unlink()
