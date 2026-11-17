#!/usr/bin/env python3
"""Line RG item T5b: the texts of a damaged taxonomy database in dustmasker, segmasker and blastdbcmd.

Writes, from LOSAT/ (`python3 tests/db_search/t5b/make_t5b.py NCBI_BIN`, NCBI_BIN the directory of
NCBI BLAST+ 2.17.0's makeblastdb), in `tests/db_search/t5b/`: the version 5 databases `nm` and `pm`
(volume 1) and `nmb` and `pmb` (volume 2), with `-parse_seqids` and the taxids of the records below
by `-taxid_map` (the letters of `../ids/ndup.fna` and `../ids/pdup.faa`), the alias files `nmv`
and `pmv` over both volumes, and `batch.txt` (the ids of the sequences with taxid 99990 alone).

The taxonomy databases are those of `tests/seqdb/taxdb/<variant>` (`make_taxdb.py`): 9606 is named,
99990 is their last entry (in `offset_last_past_end` its offset is past the end of `taxdb.btd`:
NCBI's offset error at each lookup that asks the taxonomy database for it), 150 takes the same slot
of NCBI's 256-slot taxonomy cache of a volume as 99990 (`99990 & 255`) and is not named, nor is
77777. Volume 1: two titled sequences without a taxid (no lookup), titled 9606, untitled 9606,
untitled 99990, titled 150, titled 99990, one of two lines 99990 ^A 150 (99990 hits the slot that
the sequence before filled; blastdbcmd's second `GetBioseq` of it misses after 150), untitled
77777; volume 2: untitled 99990 (a new volume's cache), titled 9606.

`python3 tests/db_search/t5b/make_t5b.py --taxdb` (no NCBI program) writes only the two taxonomy
databases of line RG item T5b's fix 2 in `tests/db_search/t5b/taxdb/`, the `taxdb.btd` of
`tests/seqdb/taxdb/good` (266 bytes, one 4,096-byte page of NCBI's mapping) with another index:
`kingdom_past_page`, where 10090's offset is 5,000, so that the record of 9606 runs from 68 to 5,000,
past the mapping after its third tab (its names inside it, its kingdom not); `names_past_page`,
where 9606's offset is 4,096 and 10090's 5,000, so that the record of 9606 lies wholly past the
mapping (NCBI looks for its tabs in memory past it). Fix 3 adds `zero_kingdom_past_page`, an entry
of taxid 0 at 147 before 562 at 5,000, so that taxid 0's record runs past the mapping after its
third tab (blastdbcmd `-tax_info` writes `N/A` for taxid 0 whatever its lookup reads), and
`dup_last_btd_empty`, an empty `taxdb.btd` (NCBI maps it as a null pointer) with good's index and
its last taxid 99990 again at 266, so that the search for 99990 ends on the entry before the last
one and reads a record (220..266): no lookup writes the offset error and the first that reads a
record faults (the searches' prefetch of a subject with taxid 9606 or 99990).

Python 3 standard library only; makeblastdb is NCBI's (the comparison oracle builds the inputs).
"""
import os
import shutil
import struct
import subprocess
import sys

CTRL_A = "\x01"
# The records of each volume: [(Seq-id, title, taxid or None), ...]. `{k}` is 1 for the nucleotide
# databases and 2 for the protein ones.
VOLUMES = {
    "m": [
        [("gb|AB4{k}0000.1|", "no taxid title 0", None)],
        [("gb|AB4{k}0001.1|", "no taxid title 1", None)],
        [("gb|AB4{k}0002.1|", "titled human", 9606)],
        [("gb|AB4{k}0003.1|", "", 9606)],
        [("gb|AB4{k}0004.1|", "", 99990)],
        [("gb|AB4{k}0005.1|", "t150", 150)],
        [("gb|AB4{k}0006.1|", "t99990 again", 99990)],
        [("gb|AB4{k}0007.1|", "first", 99990), ("gb|AB4{k}0008.1|", "second", 150)],
        [("gb|AB4{k}0009.1|", "", 77777)],
    ],
    "mb": [
        [("gb|AB4{k}0010.1|", "", 99990)],
        [("gb|AB4{k}0011.1|", "titled human b", 9606)],
    ],
}


def letters(path):
    with open(path) as f:
        return "".join(line.strip() for line in f if not line.startswith(">"))


# The index of `tests/seqdb/taxdb/good` (taxid, offset into its taxdb.btd), and for each variant the
# offsets it changes or adds (the index sorted by taxid), the entries it appends after them, and
# whether its taxdb.btd is good's (else empty).
GOOD_INDEX = [(562, 0), (7777, 41), (9606, 68), (10090, 105), (12908, 147), (32630, 196), (99990, 220)]
TAXDB_VARIANTS = {
    "kingdom_past_page": ({10090: 5000}, [], True),
    "names_past_page": ({9606: 4096, 10090: 5000}, [], True),
    "zero_kingdom_past_page": ({0: 147, 562: 5000}, [], True),
    "dup_last_btd_empty": ({}, [(99990, 266)], False),
}


def taxdb_variants():
    good = os.path.join("tests", "seqdb", "taxdb", "good")
    for name, (offsets, appended, good_btd) in TAXDB_VARIANTS.items():
        out = os.path.join("tests", "db_search", "t5b", "taxdb", name)
        os.makedirs(out, exist_ok=True)
        entries = dict(GOOD_INDEX)
        entries.update(offsets)
        entries = sorted(entries.items()) + appended
        index = struct.pack(">II", 0x8739, len(entries)) + bytes(16)
        for taxid, offset in entries:
            index += struct.pack(">II", taxid, offset)
        with open(os.path.join(out, "taxdb.bti"), "wb") as f:
            f.write(index)
        if good_btd:
            shutil.copyfile(os.path.join(good, "taxdb.btd"), os.path.join(out, "taxdb.btd"))
        else:
            open(os.path.join(out, "taxdb.btd"), "wb").close()


def main():
    if sys.argv[1:] == ["--taxdb"]:
        taxdb_variants()
        return
    ncbi_bin = sys.argv[1]
    here = os.path.join("tests", "db_search", "t5b")
    seqs = {"n": letters(os.path.join("tests", "db_search", "ids", "ndup.fna")),
            "p": letters(os.path.join("tests", "db_search", "ids", "pdup.faa"))}
    work = os.path.join(here, "work")
    os.makedirs(work, exist_ok=True)
    env = {"PATH": "/usr/bin:/bin", "HOME": os.path.abspath(work), "BLAST_USAGE_REPORT": "false"}
    for kind, dbtype, k in (("n", "nucl", 1), ("p", "prot", 2)):
        for volume, records in VOLUMES.items():
            name = kind + volume
            fasta = os.path.join(work, f"{name}.fa")
            taxmap = os.path.join(work, f"{name}.map")
            with open(fasta, "w") as f, open(taxmap, "w") as m:
                for lines in records:
                    parts = []
                    for seqid, title, taxid in lines:
                        seqid = seqid.format(k=k)
                        parts.append(seqid + (" " + title if title else ""))
                        if taxid is not None:
                            m.write(f"{seqid.split('|')[1]}\t{taxid}\n")
                    f.write(">" + CTRL_A.join(parts) + "\n" + seqs[kind] + "\n")
            cmd = [os.path.join(ncbi_bin, "makeblastdb"), "-in", fasta, "-dbtype", dbtype,
                   "-parse_seqids", "-blastdb_version", "5", "-taxid_map", taxmap,
                   "-title", name, "-out", os.path.join(here, name)]
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, env=env)
        ext = "nal" if kind == "n" else "pal"
        with open(os.path.join(here, f"{kind}mv.{ext}"), "w") as f:
            f.write(f"TITLE {kind}mv\nDBLIST {kind}m {kind}mb \n")
    with open(os.path.join(here, "batch.txt"), "w") as f:
        f.write("AB410004.1\nAB410006.1\nAB410010.1\n")
    for name in os.listdir(work):
        os.remove(os.path.join(work, name))
    os.rmdir(work)


if __name__ == "__main__":
    main()
