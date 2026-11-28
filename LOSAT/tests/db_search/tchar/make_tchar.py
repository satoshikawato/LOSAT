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
  cache of 256 header records, have such a byte and hit the queries; OID 2 has one and OID 258, of
  its slot, is clean, each hit by a query of its own; OID 3 has one; the others are random and clean);
- `tc_big2_p` / `tc_big2_n`: BIG2 (300; OID 215 has such a byte: its OID in the alias databases
  `tc_bigal_p.pal` / `tc_bigal_n.nal`, `DBLIST tc_big tc_big2`, is 515, of the slot of BIG's OID 3
  were the slots taken by the alias's OIDs, not the volume's);
- `tc_x_p` / `tc_x_n`: DBX (a title of 300 bytes, one with a NUL, one of 1,100 bytes);
- `tc_dmg_p` / `tc_dmg_n`: DMG (OID 0 has such a byte, OID 1's title is stored with the tag of a
  GeneralString, 0x1B, which NCBI cannot read: a query that hits both fails after OID 0's posts);
- `tc_ps_p` / `tc_ps_n`: PS (`-parse_seqids`, `lcl|` ids; the title of one has such a byte);
- `tc_u8_p` / `tc_u8_n`: U8 (OID 0 a clean title stored as a UTF8String, tag 0x0C, for which NCBI
  warns once; OID 1 such a byte; OID 2 a UTF8String with such a byte);
- `tc_tax_p` / `tc_tax_n`: the subjects of `tc` with `-taxid 9606` (the loads' taxonomy lookups);
- `tc_fam_p` / `tc_fam_n`: FAM (14 titles with such a byte in three families, one query each:
  `tc_qabc`, `tc_qbac`, `tc_qacb`, `tc_qaba`, the batches' order of NCBI's cache of 10 released sequences);
- `tc_du_p` / `tc_du_n`: DU (a UTF8String title, such a byte and a title NCBI cannot read, one query
  `tc_qdu` hits them in that order);
- the queries `tc_q260` (BIG's OIDs 0 to 259 in one query), `tc_qslotu` (BIG's OID 258 with OIDs
  4 to 13, then OIDs 2 and 258) and `tc_qslotu1` (its first query);
- the queries `tc_q1.faa` (one protein), `tc_q2.faa` (two), `tc_qslot.faa` (eight: OID 0's and OID
  256's letters, which read their slot in both orders, OID 2's, OID 258's and OID 2's again, OID 3's,
  BIG2's OID 215's and OID 3's again), `tc_qu8.faa` (U8's OID 0, then OID 1), `tc_qkeep.faa`
  (three: BIG's OID 2's letters, then the letters of its clean OIDs 10 to 20 in one query, eleven
  subjects whose release pushes OID 2 out of NCBI's cache of 10 released blobs, then OID 2's again),
  and their back-translations `tc_q1.fna`, `tc_q2.fna`, `tc_qslot.fna`, `tc_qu8.fna`,
  `tc_qkeep.fna`.

`make_tchar.py NCBI_BIN STEM ...` writes only the databases of those stems (`tc_big`, ...) and the
queries; `write_queries()` writes the queries alone.

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
# The letters of the slot cases: OID 2 (C) and OID 258 (D) of BIG, OID 3 (H) of BIG, OID 215 (I) of BIG2,
# U8's clean OID 0 (E); random, so that each query hits its own.
C, D, H, I, E = (lcg(80, seed) for seed in (302, 303, 301, 304, 307))
BIG = []
for i in range(300):
    if i in (0, 256):
        BIG.append((b"g%d bad#a" % i, [b"\xe9"], mut(PROT, 2 if i == 0 else 6)))
    elif i in (1, 257):
        BIG.append((b"g%d bad#b" % i, [b"\x01"], mut(PROT, 4 if i == 1 else 8)))
    elif i == 2:
        BIG.append((b"g2 cee#", [b"\xff"], C))
    elif i == 258:
        BIG.append((b"g258 dee clean", [], D))
    elif i == 3:
        BIG.append((b"g3 eta#", [b"\x7f"], H))
    else:
        BIG.append((b"g%d clean" % i, [], lcg(60, 1000 + i)))
BIG2 = [(b"h215 iota#", [b"\xe9"], I) if i == 215 else (b"h%d clean" % i, [], lcg(60, 5000 + i))
        for i in range(300)]
DBX = [(b"x1 long " + b"y" * 290 + b" en#", [b"\xe9"], PROT),
       (b"x2 nul#here", [b"\x00"], mut(PROT, 3)),
       (b"x3 huge " + b"z" * 1090 + b" #", [b"\x01"], mut(PROT, 5))]
DMG = [(b"d0 fam#one", [b"\xe9"], PROT), (b"d1 damaged title", [], mut(PROT, 10)),
       (b"d2 clean", [], lcg(80, 305))]
# The FASTA ids of PS, which makeblastdb parses: the stored titles are the rest of the line.
PS = [(b"fam#one", [b"\xe9"], PROT), (b"clean two", [], lcg(80, 306))]
PS_IDS = [b"lcl|ps1", b"lcl|ps2"]
U8 = [(b"u0 utf8 clean", [], E), (b"u1 fam#one", [b"\xe9"], PROT),
      (b"u2 utf8 fam#bee", [b"\xe9"], mut(PROT, 10))]
QUERIES = [(b"q1 query one", PROT), (b"q2 query two", mut(PROT, 2))]
SLOT_QUERIES = [(b"qs%d" % (k + 1), letters)
                for k, letters in enumerate([mut(PROT, 2), mut(PROT, 6), C, D, C, H, I, H])]
U8_QUERIES = [(b"qu1 utf8", E), (b"qu2 fam", PROT)]
KEEP_QUERIES = [(b"qk1 cee", C), (b"qk2 eleven clean", "".join(BIG[i][2] for i in range(10, 21))),
                (b"qk3 cee again", C)]
# FAM: 14 subjects in three families (A: OIDs 0-5, B: 6-10, C: 11-13), each title with such a byte; one
# query per family, so that BATCH_SIZE gives each its batch (the blob cache across batches).
FAM_A, FAM_B, FAM_C = (lcg(60, seed) for seed in (501, 502, 503))
FAM = [(b"f%02d fam#" % i, [b"\xe9"], mut(FAM_A if i < 6 else FAM_B if i < 11 else FAM_C, 1 + i % 6))
       for i in range(14)]
FAM_QUERIES = {"A": (b"qa famA", FAM_A), "B": (b"qb famB", FAM_B), "C": (b"qc famC", FAM_C)}
# DU: a UTF8String title (v0, tag 0x0C), such a byte (v1) and a title NCBI cannot read (v2, tag 0x1B), which
# one query hits in that order (v0's 120 letters score first).
EU = lcg(120, 308)
DU = [(b"v0 utf8 clean long", [], EU), (b"v1 fam#one", [b"\xe9"], PROT), (b"v2 damaged title", [], mut(PROT, 10))]
DATABASES = {"tc": SUBJECTS, "tc_b": DBB, "tc_big": BIG, "tc_big2": BIG2, "tc_x": DBX, "tc_dmg": DMG,
             "tc_ps": PS, "tc_u8": U8, "tc_tax": SUBJECTS, "tc_fam": FAM, "tc_du": DU}
# More queries: one of the letters of BIG's OIDs 0 to 259 (a batch of 260 Seq-align ids, NCBI's prefetch in
# pieces of 100, 100 and 60); BIG's OID 258 with OIDs 4 to 13 (11 loads, more than NCBI's cache of 10 keeps),
# then OID 2 and OID 258 (OID 2 posts in the header cache slot of OID 258), and its first query alone (the
# rows that LOSAT writes before it refuses the second batch); the families in the orders A B C, B A C, A C B and
# A B A; DU's v0 and v1 letters.
MORE_QUERIES = {
    "tc_q260": [(b"q260 big", "".join(BIG[i][2] for i in range(260)))],
    "tc_qslotu": [(b"qsu1 dee and ten", D + "".join(BIG[i][2] for i in range(4, 14))), (b"qsu2 cee dee", C + D)],
    "tc_qslotu1": [(b"qsu1 dee and ten", D + "".join(BIG[i][2] for i in range(4, 14)))],
    "tc_qabc": [FAM_QUERIES[c] for c in "ABC"],
    "tc_qbac": [FAM_QUERIES[c] for c in "BAC"],
    "tc_qacb": [FAM_QUERIES[c] for c in "ACB"],
    "tc_qaba": [FAM_QUERIES[c] for c in "ABA"],
    "tc_qdu": [(b"qdu utf8 then fam", EU + PROT)],
}


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


def retag(path, titles, old, new):
    """Rewrite the tag `old` of each stored title of `titles` (raw bytes) in `path` to `new`."""
    data = open(path, "rb").read()
    for title in titles:
        pattern = old + ber_length(len(title)) + title
        assert data.count(pattern) == 1, (path, title)
        data = data.replace(pattern, new + pattern[1:])
    open(path, "wb").write(data)


def write_queries():
    os.makedirs(HERE, exist_ok=True)
    for suffix, encode in ((".faa", lambda p: p), (".fna", nuc)):
        open(os.path.join(HERE, "tc_q1" + suffix), "wb").write(fasta([(QUERIES[0][0], encode(QUERIES[0][1]))]))
        open(os.path.join(HERE, "tc_q2" + suffix), "wb").write(fasta([(t, encode(s)) for t, s in QUERIES]))
        open(os.path.join(HERE, "tc_qslot" + suffix), "wb").write(fasta([(t, encode(s)) for t, s in SLOT_QUERIES]))
        open(os.path.join(HERE, "tc_qu8" + suffix), "wb").write(fasta([(t, encode(s)) for t, s in U8_QUERIES]))
        open(os.path.join(HERE, "tc_qkeep" + suffix), "wb").write(fasta([(t, encode(s)) for t, s in KEEP_QUERIES]))
        for name, records in MORE_QUERIES.items():
            open(os.path.join(HERE, name + suffix), "wb").write(fasta([(t, encode(s)) for t, s in records]))


def main():
    ncbi_bin = sys.argv[1]
    stems = sys.argv[2:] or list(DATABASES)
    write_queries()
    for stem in stems:
        subjects = DATABASES[stem]
        for mol, dbtype, encode in (("p", "prot", lambda p: p), ("n", "nucl", nuc)):
            name = f"{stem}_{mol}"
            source = os.path.join(HERE, name + ".in")
            titles = [b"%s %s" % (i, t) for i, (t, _, _) in zip(PS_IDS, subjects)] if stem == "tc_ps" else \
                [t for t, _, _ in subjects]
            open(source, "wb").write(fasta([(t, encode(s)) for t, (_, _, s) in zip(titles, subjects)]))
            subprocess.run([os.path.join(ncbi_bin, "makeblastdb"), "-in", source, "-dbtype", dbtype,
                            "-blastdb_version", "4", "-title", name, "-out", os.path.join(HERE, name)]
                           + (["-parse_seqids"] if stem == "tc_ps" else [])
                           + (["-taxid", "9606"] if stem == "tc_tax" else []),
                           check=True, stdout=subprocess.DEVNULL)
            os.remove(source)
            for ext in ("pjs", "njs"):
                if os.path.exists(os.path.join(HERE, f"{name}.{ext}")):
                    os.remove(os.path.join(HERE, f"{name}.{ext}"))
            header = os.path.join(HERE, f"{name}.{mol}hr")
            patch(header, subjects)
            if stem == "tc_dmg":
                retag(header, [b"d1 damaged title"], b"\x1a", b"\x1b")
            if stem == "tc_u8":
                retag(header, [b"u0 utf8 clean", raw_title(*U8[2][:2])], b"\x1a", b"\x0c")
            if stem == "tc_du":
                retag(header, [b"v0 utf8 clean long"], b"\x1a", b"\x0c")
                retag(header, [b"v2 damaged title"], b"\x1a", b"\x1b")
    for mol in "pn":
        open(os.path.join(HERE, f"tc_alias_{mol}.{mol}al"), "w").write(
            f"TITLE tc alias\nDBLIST tc_{mol} tc_b_{mol}\n")
        open(os.path.join(HERE, f"tc_bigal_{mol}.{mol}al"), "w").write(
            f"TITLE tc big alias\nDBLIST tc_big_{mol} tc_big2_{mol}\n")


if __name__ == "__main__":
    main()
