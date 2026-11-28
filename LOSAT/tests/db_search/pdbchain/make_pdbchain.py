#!/usr/bin/env python3
"""Line RG item K4j: a protein database whose `pdb` Seq-ids hold a chain INTEGER outside 0-255.

Writes, from LOSAT/ (`python3 tests/db_search/pdbchain/make_pdbchain.py NCBI_BIN [NAME ...]`, NCBI_BIN the
directory of NCBI BLAST+ 2.17.0's makeblastdb, NAME one of DATABASES, all by default), the version 4 database
`tests/db_search/pdbchain/pdbchain` (`.phr .pin .psq .pjs .pog .psd .psi`) of CHAINS: one short protein per
record (`pdup.faa` with one letter changed), with an EMPTY title and the Seq-ids `gb|AB60000<n>.1|` and
`pdb { mol "1AB<n>", chain <N> }` with NO chain-id; and (line RG-KB item 12) `tests/db_search/pdbchain/pdbonly`,
the same records with the pdb Seq-id alone, and `tests/db_search/pdbchain/pdbnorm`, the pdb Seq-ids of NORMALIZED
alone (each record's taxid from `-taxid_map`, TAXIDS): ids whose handle `Normalize` builds anew (a lower-case
molecule, a release date that normalizes) or whose chain is in 128..255, and ids whose handle keeps them (a date
that does not normalize: month 16 or a season; a normalized chain that equals the stored one).

Why: NCBI's `CDeflineGenerator::x_SetTitleFromPDB` (c++/src/objmgr/util/create_defline.cpp:1842-1843)
makes the title of a `pdb` Seq-id without a chain-id from `isprint` of `m_PDBChain` cast to
`unsigned char` and a one-character string of `m_PDBChain` cast to `char`: the low byte of the chain
INTEGER. So chain 321 (0x141) and -191 (0xFF41) give `Chain A`, 383 (0x17F) is not printable and
gives no PDB title, 65 gives `Chain A`, and an absent chain (the ASN.1 DEFAULT 32, which is not
"set") gives no PDB title either.

makeblastdb cannot write such a chain: it writes `chain` = the letter and `chain-id` = the letter
(`pdb|1AB1|A`), so the database is made by makeblastdb with `-parse_seqids -blastdb_version 4` from the
FASTA ids `gb|AB60000<n>.1|pdb|1AB<n>|<letter>` and its stored deflines (`.phr`, BER with indefinite
lengths) are then rewritten: the `chain` INTEGER of the pdb Seq-id set to the wanted value, the
`chain-id` removed, and the header offsets of the index (`.pin`) written again (the approach of
`../t5/make_t5.py` `patch_first_taxid`).

In `pdbchain` the gb Seq-id comes first, so that the report loads the sequence by it. In `pdbonly` a record
whose only Seq-id is such a pdb one is found by NCBI's search but its sequence is not loaded for the report: the
data loader looks up the Seq-id of the id's handle, whose chain is the normalized one (the chain's low byte read back as
a signed char: 65 for 321 and -191, 127 for 383), and serves the sequence only where a stored id compares equal
to it, which `CPDB_seq_id::Match` (chain against chain) does not for 321, -191 and 383: outfmt 0 then prints
`Unknown` and `no longer exists in database...alignment skipped`, XML `Unknown`, and the tabular formats the rows
of a failed lookup (a null pointer exception at the subject accession).

In `pdbnorm` the handle's chain is `chain-id[0]` read as a signed char for an id that `Normalize` builds anew
(c++/src/objects/seq/seq_id_tree.cpp `CSeq_id_PDB_Info::Normalize`: a lower-case molecule letter, or a `rel` Date-std
of year 1-4095, month 1-15, day 1-31, hour 0-30, minute and second 0-62 and no season), and the chain's low byte
read as a signed char for a chain without a chain-id: chain 66 with chain-id "A" then looks up chain 65, which the
stored id does not match (no Bioseq), as chain 200 without a chain-id looks up -56. Chain 66 with chain-id "A" and
a date of month 16 or a season, and chain 65 with chain-id "A" and a date of month 13, load.

In `pdbtwo` (review 2 of item 12, N3) each record has two definition lines: the first holds `pdb|1AB<n>|B` rewritten to
chain 321 without a chain-id (the Seq-align id, whose handle looks up chain 65 = `A`), the second `pdb|1AB<n>|A` (chain
65, chain-id `A`) with the title `t1`; the first line of `1AB1` has no title, that of `1AB2` the title `t0` (TWO_LINES).
NCBI's data loader serves the sequence by the second line, so the reports' Bioseq is made of it: the description line,
the organism report and the heading's generated title of an untitled line show `t1`, the heading lists both lines.
Each line has its own taxid (`-taxid_map`: 100001 for `B`, 100002 for `A`). Its keys are not patched: the first line's
`b` key leads to no lookup that the reports make.

`q_norm0.faa` is the query of the letters of `pdbnorm`'s first record (`1ab6`, which has no Bioseq): its alignment comes
first, before that of `1AB8`, whose label NCBI's BLAST XML cannot write.

`q_pdb1ab1.faa` is the query `pdb|1AB1|A` of the duplicated record's letters (decision KB-12-1 of line RG-KB item 12): read
with `-parse_deflines`, its Bioseq holds the id of the handle of `1AB1` chain 321 (chain 65), so NCBI's scope finds the
query for that sequence in `pdbonly` and `pdbtwo`.

`pdbn` (review 3 of item 12, F2) is the nucleotide database of `pdbonly`'s first two records back-translated (one codon
per letter, BACK), with the Seq-ids `pdb|1AB1|A` and `pdb|1AB2|A` as makeblastdb writes them (no rewrite), and
`q_pdb1ab1.fna` the query `pdb|1AB1|A` of the duplicated record's letters back-translated: NCBI's scope finds the query
(`q_pdb1ab1.faa` in TBLASTN, `q_pdb1ab1.fna` in BLASTN and TBLASTX) for `1AB1`.

`pdbnthrow` (line DBID item rf-2, review 1 L1) is a nucleotide database of two records: `lcl|nohit0` (letters
without a hit, NOHIT_SEED), then `pdbn`'s first record with one definition line of the Seq-ids `pdb|1AB1|B`,
rewritten to chain 321 without a chain-id (the Seq-align id, whose handle looks up chain 65 = `A`), and
`pdb|1AB1|A`, titled `t1`: NCBI's load of it throws (no definition line holds the handle's id as the load reads
the ids), but `q_pdb1ab1.faa` read with `-parse_deflines` holds that id, so NCBI's scope finds the query and
nothing is loaded. The record that the query hits is the second (rank 0 among the reported, ordinal 1). Its keys
are not patched, as `pdbtwo`'s.

The sequence index files (`.psd`, `.psi`) keep the ids as makeblastdb wrote them except that the letter
that ends a pdb key is the low byte of the stored chain (`GetChain()` cast to `char`, as
`CPDB_seq_id` forms the chain of a Seq-id without a chain-id), lower case, the keys the search looks
each hit's sequence up by.
The key lines have the same length, so `.psi` (the first key of each page) does not change.

Python 3 standard library only; makeblastdb and blastdbcmd are NCBI's (the comparison oracle builds the
inputs, under the caller's `flock` on the oracle lock).
"""
import os
import random
import struct
import subprocess
import sys

HERE = os.path.join("tests", "db_search", "pdbchain")
PDUP = os.path.join("tests", "db_search", "ids", "pdup.faa")

# (record id, FASTA chain letter, chain INTEGER stored; None: the chain left out, the DEFAULT 32).
CHAINS = [
    ("1AB1", "A", 321),
    ("1AB2", "B", -191),
    ("1AB3", "C", 383),
    ("1AB4", "D", 65),
    ("1AB5", "E", None),
]


def letters(path):
    return "".join(line.strip() for line in open(path) if not line.startswith(">"))


def tlv_end(data, at):
    """The end of the BER value at `at` (indefinite constructed values walked)."""
    at += 1
    first = data[at]
    at += 1
    if first == 0x80:
        while data[at:at + 2] != b"\x00\x00":
            at = tlv_end(data, at)
        return at + 2
    if first < 0x80:
        return at + first
    count = first - 0x80
    return at + count + int.from_bytes(data[at:at + count], "big")


def ber_integer(value):
    """The minimal two's complement content of an INTEGER."""
    size = 1
    while not -(1 << (8 * size - 1)) <= value < (1 << (8 * size - 1)):
        size += 1
    return value.to_bytes(size, "big", signed=True)


def with_chain(record, chain):
    """`record` (a `Blast-def-line-set` of one line, as makeblastdb writes it for `gb|ACC|pdb|ID|X`: title
    [0], seqid [1] holding `gb { accession }` and `pdb { mol [0], chain [1], chain-id [3] }`, taxid [2])
    with the chain INTEGER of the pdb Seq-id `chain` (None: no `chain` element, the DEFAULT 32) and no
    chain-id."""
    assert record[:4] == b"\x30\x80\x30\x80", record[:4].hex()
    at = 4
    while record[at] != 0xA1:
        at = tlv_end(record, at)
    assert record[at:at + 4] == b"\xa1\x80\x30\x80", record[at:at + 4].hex()
    at += 4  # the Seq-ids of the set
    while record[at] != 0xAE:  # pdb [14]
        assert record[at:at + 2] == b"\xa4\x80", record[at:at + 2].hex()  # gb [4]
        at = tlv_end(record, at)
    assert record[at:at + 4] == b"\xae\x80\x30\x80", record[at:at + 4].hex()
    seq = at + 2
    end = tlv_end(record, seq)
    inner = seq + 2
    mol_end = tlv_end(record, inner)
    assert record[inner] == 0xA0, record[inner:mol_end].hex()
    chain_end = tlv_end(record, mol_end)
    assert record[mol_end] == 0xA1, record[mol_end:chain_end].hex()
    id_end = tlv_end(record, chain_end)
    assert record[chain_end] == 0xA3 and id_end == end - 2, record[chain_end:end].hex()
    new = record[inner:mol_end]
    if chain is not None:
        content = ber_integer(chain)
        new += b"\xa1\x80\x02" + bytes([len(content)]) + content + b"\x00\x00"
    return record[:inner] + new + record[end - 2:]


# The records of `pdbnorm`: (molecule, FASTA chain letter, chain INTEGER, release date, chain-id). The date is
# None or the Date-std members as {name: value} (year, month, day, season).
NORMALIZED = [
    ("1ab6", "A", 66, None, "A"),
    ("1AB7", "A", 66, {"year": 2000, "month": 13}, "A"),
    ("1AB8", "C", 200, None, None),
    ("1AB9", "A", 66, {"year": 2000, "month": 16}, "A"),
    ("1ABA", "A", 66, {"year": 2000, "season": "spring"}, "A"),
    ("1ABB", "A", 65, {"year": 2000, "month": 13}, "A"),
]

# The taxid of each record of `pdbnorm` (`-taxid_map`), so that the organism report (outfmt 18) shows the loaded
# ones and leaves out those without a Bioseq.
TAXIDS = [9606, 10090, 9606, 10090, 9606, 10090]

# The records of `pdbtwo`: (molecule, the title of the first definition line or None).
TWO_LINES = [("1AB1", None), ("1AB2", "t0")]

# The records of `pdbdual` (protein) and `pdbndual` (its back-translation): (FASTA line, the k of
# `record_letters`, None for ALT's letters, and the shapes `(mol, chain)` that each pdb Seq-id of the line is
# rewritten to in its order, without a chain-id).
DUAL = [
    ("gb|AB600001.1| first", 0, []),
    ("pdb|1AB1|A|pdb|1AB1|B dual", 1, [("1AB1", 321), ("1AB1", 65)]),
    ("gb|AB600003.1| after", 2, []),
    ("gb|AB600004.1| alt", None, []),
]

# The databases: name -> whether each record has the gb Seq-id before the pdb one.
DATABASES = {
    "pdbchain": True,
    "pdbonly": False,
    "pdbnorm": False,
    "pdbtwo": False,
    "pdbn": False,
    "pdbnonly": False,
    "pdbnthrow": False,
    "pdbdual": False,
    "pdbndual": False,
}

# The seed of the letters of `pdbnthrow`'s first record (`random.Random`, A C G T, three per query letter).
NOHIT_SEED = 20261128

# One codon per letter: the back-translation of `pdbn` and `q_pdb1ab1.fna`.
BACK = {"A": "GCT", "R": "CGT", "N": "AAC", "D": "GAC", "C": "TGC", "Q": "CAG", "E": "GAA", "G": "GGC", "H": "CAC",
        "I": "ATC", "L": "CTG", "K": "AAA", "M": "ATG", "F": "TTC", "P": "CCG", "S": "TCT", "T": "ACC", "W": "TGG",
        "Y": "TAC", "V": "GTG"}

CODON = {
    a + b + c: aa
    for (a, b, c), aa in zip(
        ((a, b, c) for a in "TCAG" for b in "TCAG" for c in "TCAG"),
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    )
}


def back(protein):
    return "".join(BACK[c] for c in protein)


def revcomp(dna):
    return dna[::-1].translate(str.maketrans("ACGT", "TGCA"))


def translate_no_stop(dna):
    return "".join(CODON[dna[k:k + 3]] for k in range(0, len(dna) - 2, 3)).replace("*", "")


def ctx(tag, content):
    """A context-specific constructed value `[tag]` of indefinite length, as makeblastdb writes them."""
    return bytes([0xA0 + tag, 0x80]) + content + b"\x00\x00"


def visible(text):
    return b"\x1a" + bytes([len(text)]) + text.encode("latin-1")


def integer(value):
    content = ber_integer(value)
    return b"\x02" + bytes([len(content)]) + content


DATE_STD = ["year", "month", "day", "season", "hour", "minute", "second"]


def with_pdb(record, mol, chain, rel, chain_id):
    """`record` (a `Blast-def-line-set` of one line whose only Seq-id is `pdb`) with the pdb Seq-id
    `{ mol, chain, rel (Date std), chain-id }` (a None member left out)."""
    assert record[:4] == b"\x30\x80\x30\x80", record[:4].hex()
    at = 4
    while record[at] != 0xA1:
        at = tlv_end(record, at)
    assert record[at:at + 4] == b"\xa1\x80\x30\x80", record[at:at + 4].hex()
    at += 4
    assert record[at:at + 4] == b"\xae\x80\x30\x80", record[at:at + 4].hex()
    seq = at + 2
    end = tlv_end(record, seq)
    new = ctx(0, visible(mol))
    if chain is not None:
        new += ctx(1, integer(chain))
    if rel is not None:
        members = b"".join(ctx(tag, visible(rel[name]) if name == "season" else integer(rel[name]))
                           for tag, name in enumerate(DATE_STD) if name in rel)
        new += ctx(2, ctx(1, b"\x30\x80" + members + b"\x00\x00"))
    if chain_id is not None:
        new += ctx(3, visible(chain_id))
    return record[:seq + 2] + new + record[end - 2:]


def with_first_pdb(record, mol, chain):
    """`record` (a `Blast-def-line-set`) with its first pdb Seq-id `{ mol, chain }` (no chain-id)."""
    at = record.index(b"\xae\x80\x30\x80")
    end = tlv_end(record, at)
    new = ctx(0, visible(mol)) + ctx(1, integer(chain))
    return record[:at] + b"\xae\x80\x30\x80" + new + b"\x00\x00\x00\x00" + record[end:]


def with_pdbs(record, shapes):
    """`record` (a `Blast-def-line-set`) with its pdb Seq-ids, in their order, `{ mol, chain }` of `shapes`
    (no chain-id)."""
    at = 0
    for mol, chain in shapes:
        at = record.index(b"\xae\x80\x30\x80", at)
        end = tlv_end(record, at)
        new = b"\xae\x80\x30\x80" + ctx(0, visible(mol)) + ctx(1, integer(chain)) + b"\x00\x00\x00\x00"
        record = record[:at] + new + record[end:]
        at += len(new)
    return record


def alt_letters(n=220):
    """Letters unrelated to `pdup.faa` (a linear congruential walk over the 20 amino acids), which hit no
    other record of `pdbdual`."""
    x, out = 12345, []
    for _ in range(n):
        x = (1103515245 * x + 12345) % (1 << 31)
        out.append("ACDEFGHIKLMNPQRSTVWY"[(x >> 16) % 20])
    return "".join(out)


def patch(base, db, nucl=False):
    """Rewrites the stored deflines of the version 4 database `base` (one record each of CHAINS, or
    of NORMALIZED for `pdbnorm`, of TWO_LINES for `pdbtwo`, of CHAINS[:1] for `pdbnonly`, the second
    record of two for `pdbnthrow`, of DUAL for `pdbdual` and `pdbndual`) and the header offsets of its
    index file."""
    rows = {"pdbnorm": NORMALIZED, "pdbtwo": TWO_LINES, "pdbnonly": CHAINS[:1],
            "pdbnthrow": [None, None], "pdbdual": DUAL, "pdbndual": DUAL}.get(db, CHAINS)
    idx_ext, hdr_ext = ("nin", "nhr") if nucl else ("pin", "phr")
    index = bytearray(open(f"{base}.{idx_ext}", "rb").read())
    version, = struct.unpack(">I", index[0:4])
    assert version == 4, version
    at = 8
    for _ in range(2):  # title, date
        at += 4 + struct.unpack(">I", index[at:at + 4])[0]
    count, = struct.unpack(">I", index[at:at + 4])
    assert count == len(rows), count
    at += 4 + 8 + 4  # count, total length, longest sequence
    offsets = list(struct.unpack(f">{count + 1}I", index[at:at + 4 * (count + 1)]))
    headers = open(f"{base}.{hdr_ext}", "rb").read()
    records = [headers[offsets[k]:offsets[k + 1]] for k in range(count)]
    if db == "pdbnorm":
        records = [with_pdb(record, mol, chain, rel, chain_id)
                   for record, (mol, _, chain, rel, chain_id) in zip(records, rows)]
    elif db == "pdbtwo":
        records = [with_first_pdb(record, mol, 321) for record, (mol, _) in zip(records, rows)]
    elif db == "pdbnthrow":
        records = [records[0], with_first_pdb(records[1], "1AB1", 321)]
    elif db in ("pdbdual", "pdbndual"):
        records = [with_pdbs(record, shapes) for record, (_, _, shapes) in zip(records, rows)]
    else:
        records = [with_chain(record, chain) for record, (_, _, chain) in zip(records, rows)]
    offsets = [0]
    for record in records:
        offsets.append(offsets[-1] + len(record))
    index[at:at + 4 * (count + 1)] = struct.pack(f">{count + 1}I", *offsets)
    open(f"{base}.{hdr_ext}", "wb").write(b"".join(records))
    open(f"{base}.{idx_ext}", "wb").write(bytes(index))


def patch_keys(base, db, nucl=False):
    """Rewrites the chain letter that ends the keys of the sequence index `base.psd`/`base.nsd`
    (`1ab1 a`, `1ab1|a`, `pdb|1ab1|a`; text lines `key\\x02oid\\n`) to the low byte of the stored
    chain, for the records without a chain-id (`pdbnorm`'s with one keep the keys of their
    chain-id, the FASTA letter)."""
    sd_ext = "nsd" if nucl else "psd"
    lines = open(f"{base}.{sd_ext}", "rb").read().split(b"\n")
    assert lines[-1] == b"", lines[-1]
    out = []
    for line in lines[:-1]:
        key, oid = line.rsplit(b"\x02", 1)
        if db == "pdbnorm":
            _, letter, chain, _, chain_id = NORMALIZED[int(oid)]
            if chain_id is not None:
                out.append(line)
                continue
            letter = letter.lower().encode()
        else:
            letter = CHAINS[int(oid)][1].lower().encode()
            chain = CHAINS[int(oid)][2]
        low = bytes([32 if chain is None else chain & 0xFF]).lower()  # keys are lower case
        if key[-2:] in (b" " + letter, b"|" + letter):
            key = key[:-1] + low
        out.append(key + b"\x02" + oid)
    open(f"{base}.{sd_ext}", "wb").write(b"\n".join(out) + b"\n")


def record_letters(seq, k):
    """The letters of the database records' k-th record: `seq` with one letter changed."""
    changed = list(seq)
    at = 10 * (k + 1)
    changed[at] = "Y" if changed[at] == "W" else "W"
    return "".join(changed)


def main(ncbi_bin, names):
    work = os.path.join(HERE, "work")
    os.makedirs(work, exist_ok=True)
    env = {"PATH": "/usr/bin:/bin", "HOME": os.path.abspath(work), "BLAST_USAGE_REPORT": "false"}
    seq = letters(PDUP)
    with open(os.path.join(HERE, "q_norm0.faa"), "w") as f:
        f.write(">qn0 the letters of pdbnorm's first record\n" + record_letters(seq, 0) + "\n")
    with open(os.path.join(HERE, "q_pdb1ab1.faa"), "w") as f:
        f.write(">pdb|1AB1|A query title\n" + seq + "\n")
    with open(os.path.join(HERE, "q_pdb1ab1.fna"), "w") as f:
        f.write(">pdb|1AB1|A query title\n" + back(seq) + "\n")
    with open(os.path.join(HERE, "q_pdbn_minus.faa"), "w") as f:
        f.write(">qminus minus frame of pdbnonly\n" + translate_no_stop(revcomp(back(record_letters(seq, 0)))) + "\n")
    with open(os.path.join(HERE, "q_pdb1ab1_tail.fna"), "w") as f:
        f.write(">pdb|1AB1|A short tail query\n" + back(seq)[360:660] + "\n")
    with open(os.path.join(HERE, "q_pdb1ab1_straddle.fna"), "w") as f:
        f.write(">pdb|1AB1|A straddle query\n" + back(seq)[150:451] + "\n")
    alt = alt_letters()
    with open(os.path.join(HERE, "q_dual2.faa"), "w") as f:
        f.write(">qalt the letters of pdbdual's alt record\n" + alt + "\n>qdual the duplicated record's letters\n"
                + seq + "\n")
    with open(os.path.join(HERE, "q_dual.fna"), "w") as f:
        f.write(">qndual the duplicated record's letters back-translated\n" + back(seq) + "\n")
    with open(os.path.join(HERE, "q_dual2.fna"), "w") as f:
        f.write(">qnalt pdbndual's alt record\n" + back(alt) + "\n>qndual the duplicated record's letters "
                "back-translated\n" + back(seq) + "\n")
    for db in names:
        fasta = os.path.join(work, f"{db}.faa")
        if db in ("pdbdual", "pdbndual"):
            nucl = db == "pdbndual"
            with open(fasta, "w") as f:
                for line, k, _ in DUAL:
                    letters_k = alt if k is None else record_letters(seq, k)
                    f.write(f">{line}\n" + (back(letters_k) if nucl else letters_k) + "\n")
            base = os.path.join(HERE, db)
            subprocess.run([os.path.join(ncbi_bin, "makeblastdb"), "-in", fasta, "-dbtype",
                            "nucl" if nucl else "prot", "-parse_seqids", "-blastdb_version", "4", "-title", db,
                            "-out", base], check=True, env=env, stdout=subprocess.DEVNULL)
            patch(base, db, nucl=nucl)
            continue
        if db == "pdbn":
            with open(fasta, "w") as f:
                for k, (name, _, _) in enumerate(CHAINS[:2]):
                    f.write(f">pdb|{name}|A\n" + back(record_letters(seq, k)) + "\n")
            subprocess.run([os.path.join(ncbi_bin, "makeblastdb"), "-in", fasta, "-dbtype", "nucl", "-parse_seqids",
                            "-blastdb_version", "4", "-title", db, "-out", os.path.join(HERE, db)],
                           check=True, env=env, stdout=subprocess.DEVNULL)
            continue
        if db == "pdbnthrow":
            rng = random.Random(NOHIT_SEED)
            with open(fasta, "w") as f:
                f.write(">lcl|nohit0 no hit\n" + "".join(rng.choice("ACGT") for _ in range(3 * len(seq))) + "\n")
                f.write(">pdb|1AB1|B|pdb|1AB1|A t1\n" + back(record_letters(seq, 0)) + "\n")
            base = os.path.join(HERE, db)
            subprocess.run([os.path.join(ncbi_bin, "makeblastdb"), "-in", fasta, "-dbtype", "nucl", "-parse_seqids",
                            "-blastdb_version", "4", "-title", db, "-out", base],
                           check=True, env=env, stdout=subprocess.DEVNULL)
            patch(base, db, nucl=True)
            continue
        if db == "pdbnonly":
            name, letter, _ = CHAINS[0]
            with open(fasta, "w") as f:
                f.write(f">pdb|{name}|{letter}\n" + back(record_letters(seq, 0)) + "\n")
            base = os.path.join(HERE, db)
            subprocess.run([os.path.join(ncbi_bin, "makeblastdb"), "-in", fasta, "-dbtype", "nucl", "-parse_seqids",
                            "-blastdb_version", "4", "-title", db, "-out", base],
                           check=True, env=env, stdout=subprocess.DEVNULL)
            patch(base, db, nucl=True)
            patch_keys(base, db, nucl=True)
            continue
        if db == "pdbtwo":
            lines = [f"pdb|{mol}|B{'' if title is None else ' ' + title}\x01pdb|{mol}|A t1"
                     for mol, title in TWO_LINES]
        else:
            ids = [(mol.upper(), letter) for mol, letter, *_ in NORMALIZED] if db == "pdbnorm" else \
                [(name, letter) for name, letter, _ in CHAINS]
            lines = [(f"gb|AB60000{k + 1}.1|" if DATABASES[db] else "") + f"pdb|{name}|{letter}"
                     for k, (name, letter) in enumerate(ids)]
        with open(fasta, "w") as f:
            for k, line in enumerate(lines):
                f.write(f">{line}\n" + record_letters(seq, k) + "\n")
        base = os.path.join(HERE, db)
        taxids = []
        if db in ("pdbnorm", "pdbtwo"):
            taxid_map = os.path.join(work, f"{db}.taxids")
            with open(taxid_map, "w") as f:
                if db == "pdbnorm":
                    for (name, letter), taxid in zip(ids, TAXIDS):
                        f.write(f"pdb|{name}|{letter} {taxid}\n")
                else:
                    for mol, _ in TWO_LINES:
                        f.write(f"pdb|{mol}|B 100001\npdb|{mol}|A 100002\n")
            taxids = ["-taxid_map", taxid_map]
        subprocess.run([os.path.join(ncbi_bin, "makeblastdb"), "-in", fasta, "-dbtype", "prot", "-parse_seqids",
                        "-blastdb_version", "4", "-title", db, "-out", base] + taxids,
                       check=True, env=env, stdout=subprocess.DEVNULL)
        patch(base, db)
        if db != "pdbtwo":
            patch_keys(base, db)
    for name in os.listdir(work):
        os.remove(os.path.join(work, name))
    os.rmdir(work)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:] or list(DATABASES))
