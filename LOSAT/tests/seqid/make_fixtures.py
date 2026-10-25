#!/usr/bin/env python3
"""Build the fixtures of LOSAT/tests/seqid (tests of LOSAT/src/seqid).

NCBI BLAST+ 2.17.0 is the comparison oracle only: makeblastdb writes the databases (version 4,
-parse_seqids), blastdbcmd writes what NCBI reads back from each definition line, and blastp
and blastn write the reports whose subject id and title fields the tests compare
(expected/*.tsv). The Rust tests read the committed files and never run an NCBI program.

usage: make_fixtures.py --ncbi-bin DIR --work DIR --db-18s BASE [--out LOSAT/tests/seqid]

Inputs (fasta/): prot_ids.fa and nucl_ids.fa (every Seq-id form makeblastdb writes, several
definition lines joined by Ctrl-A, GIs, a definition line without a title, lcl|Subject_ ids)
with their taxid maps, prot_extra.fa (a GI in the second definition line only, 12 long titles,
an empty title), variants.fa (one record per further form: tpe, tpd, gpp, gibbsq, gibbmt, pdb
chains, a GI above 2^31, a 320-byte title, integer tags, patents and a pre-grant patent, nat,
prf with two fields, lcl|Subject_ ids; the gibbsq, PIG and nat records hold a sequence that no
query hits: NCBI's formatters cannot fetch their Bioseqs, and blastp's tabular output stops
with an exception), qprot.fa and qnucl.fa (the queries). BASE is NCBI's
distributed database 18S_fungal_sequences (version 5, -parse_seqids, memberships): a sample of
its header records is kept as expected/18S_sample.ber.

--work must not exist, or be a directory that this script made (it holds WORK_MARKER); the
script deletes it first.

Outputs:
  db/NAME.{pin,phr,psq} / .{nin,nhr,nsq}: the databases.
  expected/NAME.deflines.tsv: blastdbcmd -entry all, one line per definition line: ordinal id,
    %i (FASTA form of the WorstRank-best id), %a (its bare id), %g (first GI or N/A), %T
    (taxid), %e (sum of the memberships), %t (title).
  expected/NAME.asn.tsv: blastdbcmd -outfmt %d (the Blast-def-line-set in ASN.1 text), one
    line per Seq-id: ordinal id, definition line, position, the Seq-id (white space collapsed).
  expected/18S_sample.ber: records of the 18S database (u32 ordinal id, u32 length, the
    record; big-endian); expected/18S_sample.deflines.tsv as above for those ordinal ids.
  expected/CASE.fmt6.tsv: outfmt "6 qaccver sseqid sallseqid sacc saccver sallacc sgi sallgi
    stitle salltitles staxid staxids", the subject columns.
  expected/CASE.fmt5.tsv: outfmt 5, Hit_id, Hit_def and Hit_accession of each hit.
  expected/CASE.fmt8.tsv: outfmt 8, the subject Seq-id of each hit (white space collapsed).
  expected/CASE.fmt12.tsv: outfmt 12, the subject Seq-id of each hit (compact JSON).
  expected/CASE.fmt0.tsv: outfmt 0 with and without -show_gis: the description lines and the
    definition-line block of each alignment (kind, text).
Fields use backslash escapes (\\, \t, \n).
"""
import argparse
import json
import re
import shutil
import struct
import subprocess
from pathlib import Path

NCBI_ENV_DROP = ("BLASTDB",)
# One HSP per subject: each report row belongs to one sequence.
SEARCH = ["-num_threads", "1", "-max_target_seqs", "50", "-max_hsps", "1"]
WORK_MARKER = ".losat-seqid-fixtures-work"
FIELDS6 = ("qaccver sseqid sallseqid sacc saccver sallacc sgi sallgi stitle salltitles "
           "staxid staxids")
CASES = [
    # case, program, query, database
    ("blastp_prot_ids", "blastp", "qprot.fa", "prot_ids"),
    ("blastn_nucl_ids", "blastn", "qnucl.fa", "nucl_ids"),
    ("blastp_prot_extra", "blastp", "qprot.fa", "prot_extra"),
    ("blastp_variants", "blastp", "qprot.fa", "variants"),
]
# 18S records kept: the first 20 (ordinal id 9 has two definition lines), the first three
# whose GI takes five bytes, one more with several definition lines, the one with memberships
# { 12 } and the one with 23 definition lines.
SAMPLE_18S = list(range(20)) + [33, 75, 91, 131, 865, 1479]


def run(args, cwd, out=None):
    env = {k: v for k, v in __import__("os").environ.items()
           if not k.startswith("LOSAT_") and k not in NCBI_ENV_DROP}
    result = subprocess.run([str(a) for a in args], cwd=cwd, env=env, check=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if out is not None:
        Path(out).write_bytes(result.stdout)
    return result.stdout.decode()


def escape(text):
    return text.replace("\\", "\\\\").replace("\t", "\\t").replace("\n", "\\n")


def write_tsv(path, rows):
    with open(path, "w") as f:
        for row in rows:
            f.write("\t".join(escape(field) for field in row) + "\n")


def deflines(bin_dir, base, dbtype, cwd, oids=None):
    # '%' separates the fields: no title of these databases holds one.
    text = run([bin_dir / "blastdbcmd", "-db", base, "-dbtype", dbtype, "-entry", "all",
                "-outfmt", "%o%%%i%%%a%%%g%%%T%%%e%%%t"], cwd)
    rows = [line.split("%", 6) for line in text.splitlines()]
    assert all(len(row) == 7 for row in rows)
    if oids is not None:
        rows = [row for row in rows if int(row[0]) in oids]
    return rows


def header_records(base, oids):
    """The header records of `oids` of a version 4 or 5 volume BASE.nin / BASE.nhr."""
    index = Path(str(base) + ".nin").read_bytes()
    version, = struct.unpack(">I", index[0:4])
    pos = 8 + (4 if version == 5 else 0)
    for _ in range(3 if version == 5 else 2):  # title, (LMDB name,) date
        length, = struct.unpack(">I", index[pos:pos + 4])
        pos += 4 + length
    num_oids, = struct.unpack(">I", index[pos:pos + 4])
    pos += 4 + 8 + 4
    offsets = struct.unpack(">%dI" % (num_oids + 1), index[pos:pos + 4 * (num_oids + 1)])
    header = Path(str(base) + ".nhr").read_bytes()
    return [(oid, header[offsets[oid]:offsets[oid + 1]]) for oid in oids]


def fmt5(text):
    unescape = lambda s: (s.replace("&gt;", ">").replace("&lt;", "<").replace("&quot;", '"')
                          .replace("&apos;", "'").replace("&amp;", "&"))
    rows = []
    for hit in re.findall(r"<Hit>(.*?)</Hit>", text, re.S):
        row = [re.search(r"<%s>(.*?)</%s>" % (tag, tag), hit, re.S) for tag in
               ("Hit_id", "Hit_def", "Hit_accession")]
        rows.append([unescape(m.group(1)) if m else "" for m in row])
    return rows


def fmt8(text):
    """The second element of each `ids { query, subject }` of the Seq-aligns, once per hit."""
    rows = []
    for match in re.finditer(r"ids \{", text):
        depth, start, parts = 1, match.end(), []
        i = start
        while depth:
            c = text[i]
            if c == '"':
                i = text.index('"', i + 1)
            elif c in "{":
                depth += 1
            elif c == "}":
                depth -= 1
            elif c == "," and depth == 1:
                parts.append(text[start:i])
                start = i + 1
            i += 1
        parts.append(text[start:i - 1])
        subject = " ".join(parts[1].split())
        if not rows or rows[-1][0] != subject:
            rows.append([subject])
    return rows


def fmt12(text):
    """The second element of each `ids` pair of the Seq-aligns, once per hit (one JSON document
    per query)."""
    decoder, pos, rows = json.JSONDecoder(), 0, []

    def walk(value):
        if isinstance(value, dict):
            for key, member in value.items():
                if key == "ids" and isinstance(member, list) and len(member) == 2:
                    subject = json.dumps(member[1], separators=(",", ":"))
                    if not rows or rows[-1][0] != subject:
                        rows.append([subject])
                else:
                    walk(member)
        elif isinstance(value, list):
            for member in value:
                walk(member)

    while True:
        while pos < len(text) and text[pos].isspace():
            pos += 1
        if pos == len(text):
            return rows
        document, pos = decoder.raw_decode(text, pos)
        walk(document)


def split_top(text):
    """The elements of ASN.1 text `text` separated by commas outside braces and strings."""
    parts, depth, start, i = [], 0, 0, 0
    while i < len(text):
        c = text[i]
        if c == '"':
            i = text.index('"', i + 1)
            while text[i + 1:i + 2] == '"':  # a doubled quote inside the string
                i = text.index('"', i + 2)
        elif c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
        elif c == "," and depth == 0:
            parts.append(text[start:i])
            start = i + 1
        i += 1
    parts.append(text[start:])
    return [part.strip() for part in parts if part.strip()]


def inside_braces(text):
    text = text.strip()
    assert text[0] == "{" and text[-1] == "}", text
    return text[1:-1]


def asn_seqids(text):
    """The Seq-ids of blastdbcmd -outfmt %d: ordinal id, definition line, position, text."""
    rows = []
    for oid, document in enumerate(text.split("Blast-def-line-set ::= ")[1:]):
        for line, defline in enumerate(split_top(inside_braces(document))):
            for member in split_top(inside_braces(defline)):
                if member.startswith("seqid "):
                    ids = split_top(inside_braces(member[len("seqid "):]))
                    for at, seqid in enumerate(ids):
                        rows.append([str(oid), str(line), str(at), " ".join(seqid.split())])
    return rows


def fmt0(text, kind):
    rows = []
    lines = text.split("\n")
    in_descriptions = False
    for i, line in enumerate(lines):
        if line.startswith("Sequences producing significant alignments:"):
            in_descriptions = True
            continue
        if in_descriptions:
            if line.startswith(">") or line.startswith("Lambda"):
                in_descriptions = False
            elif line.strip():
                rows.append([kind + "_description", line])
                continue
        if line.startswith(">"):
            block = []
            for follow in lines[i:]:
                if follow.startswith("Length="):
                    break
                block.append(follow)
            rows.append([kind + "_block", "\n".join(block)])
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ncbi-bin", required=True, type=Path)
    parser.add_argument("--work", required=True, type=Path)
    parser.add_argument("--db-18s", required=True, type=Path)
    parser.add_argument("--out", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    bin_dir, work, out = args.ncbi_bin.resolve(), args.work.resolve(), args.out.resolve()
    if work.exists():
        if any(work.iterdir()) and not (work / WORK_MARKER).exists():
            parser.error(f"--work {work} is not empty and was not made by this script")
        shutil.rmtree(work)
    work.mkdir(parents=True)
    (work / WORK_MARKER).touch()
    for f in (out / "fasta").iterdir():
        shutil.copyfile(f, work / f.name)

    def makedb(name, fasta, dbtype, *extra):
        run([bin_dir / "makeblastdb", "-in", fasta, "-dbtype", dbtype, "-out", name,
             "-parse_seqids", "-blastdb_version", "4", "-title", name, *extra], work)

    makedb("prot_ids", "prot_ids.fa", "prot", "-taxid_map", "taxmap_prot.txt")
    makedb("nucl_ids", "nucl_ids.fa", "nucl", "-taxid_map", "taxmap_nucl.txt")
    makedb("prot_extra", "prot_extra.fa", "prot")
    makedb("variants", "variants.fa", "prot")
    (out / "db").mkdir(exist_ok=True)
    (out / "expected").mkdir(exist_ok=True)
    for name, dbtype, letter in (("prot_ids", "prot", "p"), ("nucl_ids", "nucl", "n"),
                                 ("prot_extra", "prot", "p"), ("variants", "prot", "p")):
        for ext in ("in", "hr", "sq"):
            shutil.copyfile(work / f"{name}.{letter}{ext}", out / "db" / f"{name}.{letter}{ext}")
        write_tsv(out / "expected" / f"{name}.deflines.tsv", deflines(bin_dir, name, dbtype, work))
        dump = run([bin_dir / "blastdbcmd", "-db", name, "-dbtype", dbtype, "-entry", "all",
                    "-outfmt", "%d"], work)
        write_tsv(out / "expected" / f"{name}.asn.tsv", asn_seqids(dump))

    db18s = args.db_18s.resolve()
    with open(out / "expected" / "18S_sample.ber", "wb") as f:
        for oid, record in header_records(db18s, SAMPLE_18S):
            f.write(struct.pack(">II", oid, len(record)) + record)
    write_tsv(out / "expected" / "18S_sample.deflines.tsv",
              deflines(bin_dir, db18s, "nucl", work, set(SAMPLE_18S)))

    for case, program, query, db in CASES:
        search = [bin_dir / program, "-query", query, "-db", db, *SEARCH]
        tabular = run(search + ["-outfmt", "6 " + FIELDS6], work)
        write_tsv(out / "expected" / f"{case}.fmt6.tsv",
                  [line.split("\t")[1:] for line in tabular.splitlines()])
        write_tsv(out / "expected" / f"{case}.fmt5.tsv", fmt5(run(search + ["-outfmt", "5"], work)))
        write_tsv(out / "expected" / f"{case}.fmt8.tsv", fmt8(run(search + ["-outfmt", "8"], work)))
        write_tsv(out / "expected" / f"{case}.fmt12.tsv",
                  fmt12(run(search + ["-outfmt", "12"], work)))
        write_tsv(out / "expected" / f"{case}.fmt0.tsv",
                  fmt0(run(search + ["-outfmt", "0"], work), "plain")
                  + fmt0(run(search + ["-outfmt", "0", "-show_gis"], work), "show_gis"))


if __name__ == "__main__":
    main()
