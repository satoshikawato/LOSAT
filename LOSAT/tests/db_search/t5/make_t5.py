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

With `--mix-only` (`python3 tests/db_search/t5/make_t5.py NCBI_BIN --mix-only`) it writes only the
version 5 databases `nmix` and `pmix` (fix 1 of the review of T5: the cases the rows above do not
reach) and leaves everything else as it is: one sequence for each record of MIX (the letters of
`ndup.fna` / `pdup.faa`), with `-parse_seqids` and the taxids of MIX by `-taxid_map`: several
definition lines with other taxids, the target's line without a taxid (blastdbcmd `-target_only`
takes the first source from the other lines), titled and untitled lines, and one line of each Seq-id
kind whose title NCBI makes from the source (tpg, emb, gi, pdb with and without a chain, lcl, gnl;
for the protein also sp, ref NP_ and WP_, and a title ending in a ", partial [name]" suffix that
`x_AdjustProteinTitleSuffix` rewrites from the first source). Version 5, for `-taxids` and
`-negative_taxids` with `-no_taxid_expansion` (NCBI filters no version 4 database by taxid). The
record with the suffix goes to the version 4 database `psuffix` of its own (blastdbcmd only: the
search programs refuse its outfmt 0 heading, whose generated title NCBI rewrites the same way,
`ParsedSubject::title_suffix_adjusted`).

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

CTRL_A = "\x01"
# The records of `nmix` / `pmix` / `psuffix`: ([(Seq-ids, title, taxid or None), ...], databases
# ("np": nmix and pmix, "p": pmix, "f": psuffix). `{k}` is 1 for the nucleotide database and 2 for
# the protein ones.
MIX = [
    # Two lines without titles, taxids 100001 then 100002.
    ([("gb|AB30{k}001.1|", "", 100001), ("gb|AB30{k}002.1|", "", 100002)], "np"),
    # The first line without a taxid: its source is the second line's (also with -target_only).
    ([("gb|AB31{k}001.1|", "", None), ("gb|AB31{k}002.1|", "", 100002)], "np"),
    # The first line untitled (100002), the second titled (100001).
    ([("gb|AB32{k}001.1|", "", 100002), ("gb|AB32{k}002.1|", "second title", 100001)], "np"),
    # The first line titled (100001), the second untitled (100002).
    ([("gb|AB33{k}001.1|", "first title", 100001), ("gb|AB33{k}002.1|", "", 100002)], "np"),
    ([("tpg|BK3{k}0041.1|", "", 100001)], "np"),
    ([("emb|CA3{k}0051.1|", "", 100005)], "np"),
    ([("gi|300006{k}|gb|AB34{k}061.1|", "", 100002)], "np"),
    ([("pdb|1XY{k}|", "", 100001)], "np"),
    ([("pdb|1XZ{k}|A", "", 100002)], "np"),
    ([("lcl|Subject_1", "", 100001)], "np"),
    ([("gnl|mydb|seq{k}1", "", 100004)], "np"),
    ([("sp|P99901.1|NAME_BACSP", "", 100002)], "p"),
    ([("ref|NP_300001.1|", "", 100001)], "p"),
    ([("ref|WP_300000001.1|", "", 100002)], "p"),
    # A title with the suffix of a source the first line does not have (taxid 0 / 100002).
    ([("gb|AB36{k}001.1|", "alpha, partial [Bacillus sp.]", None), ("gb|AB36{k}002.1|", "", 100002)], "f"),
]


def map_key(seqid):
    """The key of `-taxid_map` for the Seq-id `seqid`: its accession.version, a PDB molecule
    with `_` and its chain, the local or general tag."""
    parts = seqid.split("|")
    kind = parts[0]
    if kind == "gi":
        return map_key("|".join(parts[2:]))
    if kind == "pdb":
        return parts[1] + ("_" + parts[2] if len(parts) > 2 and parts[2] else "")
    if kind == "lcl":
        return parts[1]
    if kind == "gnl":
        return parts[2]
    return parts[1]


def mix_main(ncbi_bin):
    here = os.path.join("tests", "db_search", "t5")
    seqs = {"n": letters(os.path.join("tests", "db_search", "ids", "ndup.fna")),
            "p": letters(os.path.join("tests", "db_search", "ids", "pdup.faa"))}
    work = os.path.join(here, "work")
    os.makedirs(work, exist_ok=True)
    for kind, name, dbtype, k, version in (("n", "nmix", "nucl", 1, "5"), ("p", "pmix", "prot", 2, "5"),
                                           ("f", "psuffix", "prot", 2, "4")):
        fasta = os.path.join(work, f"{name}.fa")
        taxmap = os.path.join(work, f"{name}.map")
        with open(fasta, "w") as f, open(taxmap, "w") as m:
            for lines, molecules in MIX:
                if kind not in molecules:
                    continue
                parts = []
                for seqid, title, taxid in lines:
                    seqid = seqid.format(k=k)
                    parts.append(seqid + (" " + title if title else ""))
                    if taxid is not None:
                        m.write(f"{map_key(seqid)}\t{taxid}\n")
                f.write(">" + CTRL_A.join(parts) + "\n" + seqs[dbtype[0]] + "\n")
        out = os.path.join(here, name)
        cmd = [os.path.join(ncbi_bin, "makeblastdb"), "-in", fasta, "-dbtype", dbtype,
               "-parse_seqids", "-blastdb_version", version, "-taxid_map", taxmap,
               "-title", name, "-out", out]
        env = {"PATH": "/usr/bin:/bin", "HOME": os.path.abspath(work),
               "BLAST_USAGE_REPORT": "false"}
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, env=env)
        # makeblastdb's version 4 JSON metadata file is not read.
        if version == "4" and os.path.exists(out + ".pjs"):
            os.remove(out + ".pjs")
    for name in os.listdir(work):
        os.remove(os.path.join(work, name))
    os.rmdir(work)


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
    if "--mix-only" in sys.argv[2:]:
        mix_main(sys.argv[1])
    else:
        main()
