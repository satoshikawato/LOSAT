#!/usr/bin/env python3
"""Line RG item T5: database sequences without a title whose taxid a taxonomy database names.

Writes, from LOSAT/ (`python3 tests/db_search/t5/make_t5.py NCBI_BIN`, NCBI_BIN the directory of
NCBI BLAST+ 2.17.0's makeblastdb):

- `tests/seqdb/taxdb/t5.tsv` and the taxonomy database `tests/seqdb/taxdb/t5/` of it
  (`tests/seqdb/taxdb/make_taxdb.py write`): one name of each form that NCBI's
  `CDeflineGenerator` changes when it makes a title from it (TAXA below).
- in `tests/db_search/t5/`: for each taxid of TAXA and 77777 (which the taxonomy database does
  not hold), the version 4 databases `n_<taxid>` and `p_<taxid>` of one sequence without a title
  (`>` alone; the letters of `../ids/ndup.fna` and `../ids/pdup.faa`), made with `-taxid <taxid>`
  and without `-parse_seqids`, and for the taxids of PARSED `np_<taxid>` and `pp_<taxid>` with
  `-parse_seqids` (the id `gb|AB0000<n>.1|` alone); alias files over them (ALIASES).

Python 3 standard library only; makeblastdb is NCBI's (the comparison oracle builds the inputs).
"""
import os
import subprocess
import sys

TAXA = [
    # taxid, scientific name, common name, blast name, kingdom
    (100001, b"uncultured bacterium", b"", b"bacteria", b"Bacteria"),
    (100002, b"Bacillus sp.", b"", b"firmicutes", b"Bacteria"),
    (100003, b"  Spaced   name  ", b"sp", b"x", b"Bacteria"),
    (100004, b"A &amp; B &lt;x&gt; &#65;", b"", b"x", b"Bacteria"),
    (100005, b"TPA: tpa organism", b"", b"x", b"Bacteria"),
    (100006, b"MAG uncultured thing", b"", b"x", b"Bacteria"),
    # A byte that NCBI's HtmlDecode cannot place in an encoding (it throws).
    (100009, b"Bad \x81 name", b"", b"x", b"Bacteria"),
    # An empty scientific name: no name.
    (100010, b"", b"only common", b"x", b"Bacteria"),
]
UNNAMED = 77777
PARSED = [100001, 100005, 100010]
NAMED = [100001, 100002, 100003, 100004, 100005, 100006]
ALIASES = {
    # outfmt 0 and 18: every form, the unknown encoding and the taxids without a name.
    "nall": ("n", NAMED + [100009, 100010, UNNAMED]),
    "pall": ("p", NAMED + [100009, 100010, UNNAMED]),
    # blastdbcmd, dustmasker and segmasker: no title NCBI cannot decode (NCBI stops there).
    "nclean": ("n", NAMED + [100010, UNNAMED]),
    "pclean": ("p", NAMED + [100010, UNNAMED]),
    # The tabular subject ids of a nucleotide: an empty title without a name has no id.
    "ntab": ("n", NAMED),
    "npall": ("np", PARSED),
    "ppall": ("pp", PARSED),
}


def letters(path):
    with open(path) as f:
        return "".join(line.strip() for line in f if not line.startswith(">"))


def main():
    ncbi_bin = sys.argv[1]
    here = os.path.join("tests", "db_search", "t5")
    taxdir = os.path.join("tests", "seqdb", "taxdb")
    tsv = os.path.join(taxdir, "t5.tsv")
    with open(tsv, "wb") as f:
        for row in TAXA:
            f.write(b"\t".join([str(row[0]).encode()] + list(row[1:])) + b"\n")
    subprocess.run([sys.executable, os.path.join(taxdir, "make_taxdb.py"), "write", tsv,
                    os.path.join(taxdir, "t5")], check=True)
    seqs = {"n": letters(os.path.join("tests", "db_search", "ids", "ndup.fna")),
            "p": letters(os.path.join("tests", "db_search", "ids", "pdup.faa"))}
    work = os.path.join(here, "work")
    os.makedirs(work, exist_ok=True)
    taxids = [row[0] for row in TAXA] + [UNNAMED]
    for kind, dbtype in (("n", "nucl"), ("p", "prot")):
        for taxid in taxids:
            fasta = os.path.join(work, f"{kind}.fa")
            with open(fasta, "w") as f:
                f.write(">\n" + seqs[kind] + "\n")
            makeblastdb(ncbi_bin, fasta, dbtype, taxid, os.path.join(here, f"{kind}_{taxid}"), False)
        for number, taxid in enumerate(PARSED, 1):
            fasta = os.path.join(work, f"{kind}p.fa")
            with open(fasta, "w") as f:
                f.write(f">gb|AB0000{number:02d}.1|\n" + seqs[kind] + "\n")
            makeblastdb(ncbi_bin, fasta, dbtype, taxid, os.path.join(here, f"{kind}p_{taxid}"), True)
    for name in os.listdir(work):
        os.remove(os.path.join(work, name))
    os.rmdir(work)
    for name, (prefix, members) in ALIASES.items():
        ext = "nal" if prefix.startswith("n") else "pal"
        with open(os.path.join(here, f"{name}.{ext}"), "w") as f:
            f.write(f"TITLE {name}\nDBLIST {' '.join(f'{prefix}_{t}' for t in members)} \n")


def makeblastdb(ncbi_bin, fasta, dbtype, taxid, out, parse):
    cmd = [os.path.join(ncbi_bin, "makeblastdb"), "-in", fasta, "-dbtype", dbtype,
           "-blastdb_version", "4", "-taxid", str(taxid), "-title", os.path.basename(out),
           "-out", out]
    if parse:
        cmd.append("-parse_seqids")
    env = {"PATH": "/usr/bin:/bin", "HOME": os.path.abspath(os.path.dirname(out)),
           "BLAST_USAGE_REPORT": "false"}
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, env=env)
    # makeblastdb's version 4 JSON metadata file is not read by a search.
    for ext in (".njs", ".pjs"):
        if os.path.exists(out + ext):
            os.remove(out + ext)


if __name__ == "__main__":
    main()
