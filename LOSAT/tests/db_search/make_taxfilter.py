#!/usr/bin/env python3
"""Builds tests/db_search/taxfilter/: the databases and the taxonomy file of the rows db.tf.* (DI-8a:
-taxids, -taxidlist, -negative_taxids, -negative_taxidlist, -no_taxid_expansion).

Usage: make_taxfilter.py MAKEBLASTDB TAXONOMY4BLAST_SQLITE   (writes taxfilter/{tfp,tfn}.*, taxfilter/taxonomy4blast.sqlite3
and taxfilter/root/taxonomy4blast.sqlite3)

- `tfp`: version 5 protein database of `tfp.faa` (23 globins of mammals, birds, reptiles, amphibians and fishes;
  four records carry a second definition line of another species), taxids from `tfp.map`.
- `tfn`: version 5 nucleotide database of `tfn.fna` (16 globin cDNAs of human, zebrafish and fly; human HBA1
  carries a second definition line of a zebrafish transcript, zebrafish hbaa1 one of a fly transcript), taxids
  from `tfn.map`, in two volumes (`-max_file_sz 4KB`) under one LMDB.
- `taxonomy4blast.sqlite3`: the rows of NCBI's file for the taxids of both maps and all their ancestors, in a new
  file with the source's own schema text (NCBI reads the table `TaxidInfo` and its index only, so the
  descendants of a taxid within the databases are the same as with the whole file).
- `root/taxonomy4blast.sqlite3`: the same without the rows of 7227 (fly) and 8237 (tuna), whose sequences are then not
  descendants of the root: `-taxids 1` (DECISIONS 10 of DI-8: taxid 1 is its own parent, NCBI's expansion does not end)
  keeps what `-taxids 131567` (the only child of the root here) keeps, and less than the whole database.

`.pjs` and `.njs` (the metadata with the build's paths and dates; no search reads them) are dropped.
"""
import os
import sqlite3
import subprocess
import sys

HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "taxfilter")


def makeblastdb(binary, *args):
    subprocess.run([binary, *args], check=True, cwd=HERE, stdout=subprocess.DEVNULL)


def taxids(*maps):
    found = set()
    for name in maps:
        with open(os.path.join(HERE, name)) as f:
            found.update(int(line.split()[1]) for line in f if line.strip())
    return found


def trim_taxonomy(source, out, leave_out=()):
    src = sqlite3.connect(f"file:{source}?mode=ro", uri=True)
    keep = {}
    todo = sorted(taxids("tfp.map", "tfn.map"))
    while todo:
        taxid = todo.pop()
        if taxid in keep:
            continue
        row = src.execute("SELECT parent FROM TaxidInfo WHERE taxid = ?", (taxid,)).fetchone()
        if row is None:
            raise SystemExit(f"taxid {taxid} is not in {source}")
        keep[taxid] = row[0]
        todo.append(row[0])
    for taxid in leave_out:
        del keep[taxid]
    schema = [sql for (sql,) in src.execute(
        "SELECT sql FROM sqlite_master WHERE sql IS NOT NULL ORDER BY type DESC, name")]
    if os.path.exists(out):
        os.remove(out)
    dst = sqlite3.connect(out)
    for sql in schema:
        dst.execute(sql)
    dst.executemany("INSERT INTO TaxidInfo (taxid, parent) VALUES (?, ?)", sorted(keep.items()))
    dst.commit()
    dst.execute("VACUUM")
    dst.close()


def main(binary, source):
    makeblastdb(binary, "-in", "tfp.faa", "-dbtype", "prot", "-out", "tfp", "-parse_seqids", "-blastdb_version", "5",
                "-taxid_map", "tfp.map", "-title", "tfp")
    makeblastdb(binary, "-in", "tfn.fna", "-dbtype", "nucl", "-out", "tfn", "-parse_seqids", "-blastdb_version", "5",
                "-taxid_map", "tfn.map", "-max_file_sz", "4KB", "-title", "tfn")
    for name in ("tfp.pjs", "tfn.njs"):
        os.remove(os.path.join(HERE, name))
    trim_taxonomy(source, os.path.join(HERE, "taxonomy4blast.sqlite3"))
    os.makedirs(os.path.join(HERE, "root"), exist_ok=True)
    trim_taxonomy(source, os.path.join(HERE, "root", "taxonomy4blast.sqlite3"), leave_out=(7227, 8237))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
