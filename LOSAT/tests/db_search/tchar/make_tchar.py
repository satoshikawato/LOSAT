#!/usr/bin/env python3
"""Line RG-KB item 9(b): databases whose stored definition-line titles hold bytes outside ' '..'~'.

Writes, from LOSAT/ (`python3 tests/db_search/tchar/make_tchar.py NCBI_BIN`, NCBI_BIN the directory of
NCBI BLAST+ 2.17.0's makeblastdb), the version 4 databases of `tests/db_search/tchar/` (README.md of
`tests/db_search`, section "Titles with bytes outside the visible ASCII range"):

- `tc_p` / `tc_n`: SUBJECTS (seven proteins, or their back-translations), five of whose titles hold
  bytes outside ' '..'~' (one, two or none per title; 0x01, 0x1B, 0x7F, 0x80, 0xC3 0xA9, 0xE9, 0xFF;
  one title of 133 bytes, whose length takes two bytes);
- `tc_b_p` / `tc_b_n`: DBB (three, two with such a byte), the second volume of the alias databases
  `tc_alias_p.pal` / `tc_alias_n.nal` (`DBLIST tc_p tc_b_p`);
- `tc_big_p` / `tc_big_n`: BIG (300; the OIDs 0 and 256, and 1 and 257, which share a slot of NCBI's
  cache of 256 header records, have such a byte and hit the query; the others are random and clean);
- `tc_x_p` / `tc_x_n`: DBX (a title of 300 bytes, one with a NUL, one of 1,100 bytes);
- the queries `tc_q1.faa` (one protein), `tc_q2.faa` (two), `tc_q1.fna` and `tc_q2.fna` (their
  back-translations).

makeblastdb stores each byte outside ' '..'~' of a FASTA title as '#', so the databases are made from
titles with '#' in those places and the stored titles are then rewritten to the raw bytes (the same
length; `.phr` holds each title as a VisibleString, tag 0x1A, after its length).

Python 3 standard library only; makeblastdb is NCBI's (the comparison oracle builds the inputs, under the
caller's `flock` on the oracle lock).
"""
import os
import subprocess
import sys

HERE = os.path.join("tests", "db_search", "tchar")
CODON = {"A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT",
         "I": "ATT", "K": "AAA", "L": "CTG", "M": "ATG", "N": "AAT", "P": "CCG", "Q": "CAG",
         "R": "CGT", "S": "TCT", "T": "ACC", "V": "GTT", "W": "TGG", "Y": "TAT"}
PROT = "MDIVWSTFSARIQFAQYVPPKIIRFTWHNHPSLEGWNGRQFYDIGENTGMFRLGFPKFIG"
AA = "ACDEFGHIKLMNPQRSTVWY"


def lcg(n, seed):
    """n random letters of AA (a fixed linear congruential generator)."""
    out, x = [], seed
    for _ in range(n):
        x = (x * 1103515245 + 12345) % (1 << 31)
        out.append(AA[(x >> 16) % 20])
    return "".join(out)


def mut(s, n):
    """s with n substitutions at its last odd positions (the start stays: megablast's 28-mer)."""
    s = list(s)
    for k in range(n):
        i = len(s) - 1 - 2 * k
        s[i] = "A" if s[i] == "W" else "W"
    return "".join(s)


def nuc(p):
    return "".join(CODON[a] for a in p)


# (title as makeblastdb stores it, '#' for each raw byte; the raw bytes of the '#'s in order; protein)
SUBJECTS = [
    (b"a1 alpha#one", [b"\xe9"], PROT),
    (b"a2 clean two", [], mut(PROT, 3)),
    (b"a3 b#t#a two hsps", [b"\x01", b"\x7f"], mut(PROT, 5) + lcg(40, 7) + mut(PROT, 8)),
    (b"a4 gamma##", [b"\xc3", b"\xa9"], mut(PROT, 10)),
    (b"a5 nohit#", [b"\xff"], lcg(60, 99)),
    (b"#a6 first", [b"\x1b"], mut(PROT, 14)),
    (b"a7 long " + b"x" * 120 + b" end#", [b"\x80"], mut(PROT, 16)),
]
DBB = [(b"b1 beta#one", [b"\xe9"], mut(PROT, 1)), (b"b2 beta clean", [], mut(PROT, 7)),
       (b"b3 beta#three", [b"\x7f"], mut(PROT, 9))]
BIG = []
for i in range(300):
    if i in (0, 256):
        BIG.append((b"g%d bad#a" % i, [b"\xe9"], mut(PROT, 2 if i == 0 else 6)))
    elif i in (1, 257):
        BIG.append((b"g%d bad#b" % i, [b"\x01"], mut(PROT, 4 if i == 1 else 8)))
    else:
        BIG.append((b"g%d clean" % i, [], lcg(60, 1000 + i)))
DBX = [(b"x1 long " + b"y" * 290 + b" en#", [b"\xe9"], PROT),
       (b"x2 nul#here", [b"\x00"], mut(PROT, 3)),
       (b"x3 huge " + b"z" * 1090 + b" #", [b"\x01"], mut(PROT, 5))]
QUERIES = [(b"q1 query one", PROT), (b"q2 query two", mut(PROT, 2))]
DATABASES = {"tc": SUBJECTS, "tc_b": DBB, "tc_big": BIG, "tc_x": DBX}


def fasta(records):
    return b"".join(b">" + title + b"\n" + letters.encode() + b"\n" for title, letters in records)


def raw_title(stored, raws):
    out, k = bytearray(), 0
    for c in stored:
        if c == 0x23 and k < len(raws):
            out += raws[k]
            k += 1
        else:
            out.append(c)
    assert k == len(raws) and len(out) == len(stored), stored
    return bytes(out)


def ber_length(n):
    if n < 128:
        return bytes([n])
    if n < 256:
        return bytes([0x81, n])
    return bytes([0x82, n >> 8, n & 0xFF])


def patch(path, subjects):
    data = open(path, "rb").read()
    for stored, raws, _ in subjects:
        if raws:
            pattern = b"\x1a" + ber_length(len(stored)) + stored
            assert data.count(pattern) == 1, (path, stored[:20])
            data = data.replace(pattern, b"\x1a" + ber_length(len(stored)) + raw_title(stored, raws))
    open(path, "wb").write(data)


def main():
    ncbi_bin = sys.argv[1]
    os.makedirs(HERE, exist_ok=True)
    for suffix, encode in ((".faa", lambda p: p), (".fna", nuc)):
        open(os.path.join(HERE, "tc_q1" + suffix), "wb").write(fasta([(QUERIES[0][0], encode(QUERIES[0][1]))]))
        open(os.path.join(HERE, "tc_q2" + suffix), "wb").write(fasta([(t, encode(s)) for t, s in QUERIES]))
    for stem, subjects in DATABASES.items():
        for mol, dbtype, encode in (("p", "prot", lambda p: p), ("n", "nucl", nuc)):
            name = f"{stem}_{mol}"
            source = os.path.join(HERE, name + ".in")
            open(source, "wb").write(fasta([(t, encode(s)) for t, _, s in subjects]))
            subprocess.run([os.path.join(ncbi_bin, "makeblastdb"), "-in", source, "-dbtype", dbtype,
                            "-blastdb_version", "4", "-title", name, "-out", os.path.join(HERE, name)],
                           check=True, stdout=subprocess.DEVNULL)
            os.remove(source)
            for ext in ("pjs", "njs"):
                if os.path.exists(os.path.join(HERE, f"{name}.{ext}")):
                    os.remove(os.path.join(HERE, f"{name}.{ext}"))
            patch(os.path.join(HERE, f"{name}.{mol}hr"), subjects)
    for mol in "pn":
        open(os.path.join(HERE, f"tc_alias_{mol}.{mol}al"), "w").write(
            f"TITLE tc alias\nDBLIST tc_{mol} tc_b_{mol}\n")


if __name__ == "__main__":
    main()
