#!/usr/bin/env python3
"""Builds tests/db_search/taxfilter/: the databases and the taxonomy file of the rows db.tf.* (DI-8a:
-taxids, -taxidlist, -negative_taxids, -negative_taxidlist, -no_taxid_expansion).

Usage: make_taxfilter.py MAKEBLASTDB TAXONOMY4BLAST_SQLITE   (writes taxfilter/{tfp,tfn,tfu_n,tfu_p}.*,
taxfilter/taxonomy4blast.sqlite3 and taxfilter/root/taxonomy4blast.sqlite3), or
make_taxfilter.py --unparsed-only MAKEBLASTDB   (taxfilter/tfu_n.*, tfu_p.*, tfu_nr.nal, tfu_pr.pal)

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

- `tfu_n` and `tfu_p` (line DBID ordtax): version 5 databases without parsed Seq-ids, one volume each, of four records of
  `tfn.fna` and `tfp.faa` under short titles (`tfu.fna`, `tfu.faa`, written here), taxids 7955, 9606, 7227 (protein 8237)
  and 9606 in OID order: `-taxids 9606` keeps OIDs 1 and 3, which the restricted search numbers 0 and 1, and NCBI names
  them `gnl|BL_ORD_ID|1` and `gnl|BL_ORD_ID|3`. makeblastdb gives one taxid to a database without parsed Seq-ids
  (`-taxid`; `-taxid_map` needs them), so each record is a database of its own first, and `-input_type blastdb` joins
  them. `tfu_nr.nal` and `tfu_pr.pal` read OIDs 1 to 3 of them (`FIRST_OID 2`, `LAST_OID 4`, which count from 1).
  `make_taxfilter.py --unparsed-only MAKEBLASTDB` builds only these.

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


# The records of `tfu_n` and `tfu_p` in OID order: (FASTA file, first word of the record, short title, taxid).
UNPARSED = {
    "n": [("tfn.fna", "ENSDART00000113098.5", "zebrafish hbaa1 cDNA", 7955),
          ("tfn.fna", "ENST00000251595.11", "human HBA2 cDNA", 9606),
          ("tfn.fna", "FBtr0083245", "fly glob1 cDNA", 7227),
          ("tfn.fna", "ENST00000320868.9", "human HBA1 cDNA", 9606)],
    "p": [("tfp.faa", "sp|Q90487.3|", "zebrafish hemoglobin alpha", 7955),
          ("tfp.faa", "sp|P02008.2|", "human hemoglobin zeta", 9606),
          ("tfp.faa", "sp|P11748.1|", "tuna hemoglobin alpha", 8237),
          ("tfp.faa", "sp|P09105.2|", "human hemoglobin theta-1", 9606)],
}


def fasta_records(name):
    records, current = {}, None
    with open(os.path.join(HERE, name)) as f:
        for line in f:
            if line.startswith(">"):
                current = line[1:].split()[0]
                records[current] = []
            else:
                records[current].append(line.strip())
    return {key: "".join(lines) for key, lines in records.items()}


def build_unparsed(binary):
    for kind, dbtype, extension, alias in (("n", "nucl", "fna", "nal"), ("p", "prot", "faa", "pal")):
        parts, fasta = [], []
        for i, (source, key, title, taxid) in enumerate(UNPARSED[kind]):
            letters = fasta_records(source)[key]
            record = f">{title}\n" + "".join(letters[at:at + 60] + "\n" for at in range(0, len(letters), 60))
            fasta.append(record)
            part = f"tfu_{kind}{i}"
            with open(os.path.join(HERE, part + ".tmp"), "w") as f:
                f.write(record)
            makeblastdb(binary, "-in", part + ".tmp", "-dbtype", dbtype, "-out", part, "-blastdb_version", "5",
                        "-taxid", str(taxid), "-title", part)
            parts.append(part)
        with open(os.path.join(HERE, f"tfu.{extension}"), "w") as f:
            f.write("".join(fasta))
        makeblastdb(binary, "-in", " ".join(parts), "-input_type", "blastdb", "-dbtype", dbtype, "-out", f"tfu_{kind}",
                    "-blastdb_version", "5", "-title", f"tfu_{kind}")
        for name in os.listdir(HERE):
            if any(name.startswith(part + ".") for part in parts):
                os.remove(os.path.join(HERE, name))
        os.remove(os.path.join(HERE, f"tfu_{kind}.{dbtype[0]}js"))
        with open(os.path.join(HERE, f"tfu_{kind}r.{alias}"), "w") as f:
            f.write(f"#\n# Alias file of line DBID, ordtax\n#\nTITLE tfu_{kind}r\nDBLIST tfu_{kind}\nFIRST_OID 2\nLAST_OID 4\n")


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
    if sys.argv[1] == "--unparsed-only":
        build_unparsed(sys.argv[2])
    else:
        main(sys.argv[1], sys.argv[2])
        build_unparsed(sys.argv[1])
