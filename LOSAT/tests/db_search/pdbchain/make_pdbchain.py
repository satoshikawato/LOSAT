#!/usr/bin/env python3
"""Line RG item K4j: a protein database whose `pdb` Seq-ids hold a chain INTEGER outside 0-255.

Writes, from LOSAT/ (`python3 tests/db_search/pdbchain/make_pdbchain.py NCBI_BIN [NAME ...]`, NCBI_BIN the
directory of NCBI BLAST+ 2.17.0's makeblastdb, NAME one of DATABASES, all by default), the version 4 database
`tests/db_search/pdbchain/pdbchain` (`.phr .pin .psq .pjs .pog .psd .psi`) of CHAINS: one short protein per
record (`pdup.faa` with one letter changed), with an EMPTY title and the Seq-ids `gb|AB60000<n>.1|` and
`pdb { mol "1AB<n>", chain <N> }` with NO chain-id; and (line RG-KB item 12) `tests/db_search/pdbchain/pdbonly`,
the same records with the pdb Seq-id alone.

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

The sequence index files (`.psd`, `.psi`) keep the ids as makeblastdb wrote them except that the letter
that ends a pdb key is the low byte of the stored chain (`GetChain()` cast to `char`, as
`CPDB_seq_id` forms the chain of a Seq-id without a chain-id), lower case, the keys the search looks
each hit's sequence up by.
The key lines have the same length, so `.psi` (the first key of each page) does not change.

Python 3 standard library only; makeblastdb and blastdbcmd are NCBI's (the comparison oracle builds the
inputs, under the caller's `flock` on the oracle lock).
"""
import os
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


# The databases: name -> whether each record has the gb Seq-id before the pdb one.
DATABASES = {"pdbchain": True, "pdbonly": False}


def patch(base):
    """Rewrites the stored deflines of the version 4 protein database `base` (one record each of CHAINS)
    and the header offsets of its index file."""
    index = bytearray(open(f"{base}.pin", "rb").read())
    version, = struct.unpack(">I", index[0:4])
    assert version == 4, version
    at = 8
    for _ in range(2):  # title, date
        at += 4 + struct.unpack(">I", index[at:at + 4])[0]
    count, = struct.unpack(">I", index[at:at + 4])
    assert count == len(CHAINS), count
    at += 4 + 8 + 4  # count, total length, longest sequence
    offsets = list(struct.unpack(f">{count + 1}I", index[at:at + 4 * (count + 1)]))
    headers = open(f"{base}.phr", "rb").read()
    records = [headers[offsets[k]:offsets[k + 1]] for k in range(count)]
    records = [with_chain(record, chain) for record, (_, _, chain) in zip(records, CHAINS)]
    offsets = [0]
    for record in records:
        offsets.append(offsets[-1] + len(record))
    index[at:at + 4 * (count + 1)] = struct.pack(f">{count + 1}I", *offsets)
    open(f"{base}.phr", "wb").write(b"".join(records))
    open(f"{base}.pin", "wb").write(bytes(index))


def patch_keys(base):
    """Rewrites the chain letter that ends the keys of the sequence index `base.psd` (`1ab1 a`, `1ab1|a`,
    `pdb|1ab1|a`; text lines `key\\x02oid\\n`) to the low byte of the stored chain."""
    lines = open(f"{base}.psd", "rb").read().split(b"\n")
    assert lines[-1] == b"", lines[-1]
    out = []
    for line in lines[:-1]:
        key, oid = line.rsplit(b"\x02", 1)
        letter = CHAINS[int(oid)][1].lower().encode()
        chain = CHAINS[int(oid)][2]
        low = bytes([32 if chain is None else chain & 0xFF]).lower()  # keys are lower case
        if key[-2:] in (b" " + letter, b"|" + letter):
            key = key[:-1] + low
        out.append(key + b"\x02" + oid)
    open(f"{base}.psd", "wb").write(b"\n".join(out) + b"\n")


def main(ncbi_bin, names):
    work = os.path.join(HERE, "work")
    os.makedirs(work, exist_ok=True)
    env = {"PATH": "/usr/bin:/bin", "HOME": os.path.abspath(work), "BLAST_USAGE_REPORT": "false"}
    seq = letters(PDUP)
    for db in names:
        fasta = os.path.join(work, f"{db}.faa")
        with open(fasta, "w") as f:
            for k, (name, letter, _) in enumerate(CHAINS):
                changed = list(seq)
                at = 10 * (k + 1)
                changed[at] = "Y" if changed[at] == "W" else "W"
                gb = f"gb|AB60000{k + 1}.1|" if DATABASES[db] else ""
                f.write(f">{gb}pdb|{name}|{letter}\n" + "".join(changed) + "\n")
        base = os.path.join(HERE, db)
        subprocess.run([os.path.join(ncbi_bin, "makeblastdb"), "-in", fasta, "-dbtype", "prot", "-parse_seqids",
                        "-blastdb_version", "4", "-title", db, "-out", base],
                       check=True, env=env, stdout=subprocess.DEVNULL)
        patch(base)
        patch_keys(base)
    for name in os.listdir(work):
        os.remove(os.path.join(work, name))
    os.rmdir(work)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:] or list(DATABASES))
