#!/usr/bin/env python3
"""Build the BLAST database fixtures of LOSAT/tests/seqdb (tests of LOSAT/src/seqdb).

NCBI BLAST+ is used here as the comparison oracle only: makeblastdb and
blastdb_aliastool write the databases, and blastdbcmd writes what NCBI reads
back from them (expected/*.tsv). The Rust tests read the committed files and
never run an NCBI program.

usage: make_fixtures.py --ncbi-bin DIR --work DIR [--out LOSAT/tests/seqdb]

Committed per database: the files LOSAT reads (.pin .psq .phr / .nin .nsq .nhr,
.pal / .nal) and the LMDB file of a version 5 database (.pdb / .ndb), which must
exist. The other files of a version 5 database (.pot/.ptf/.pto, .pos, .pog, .pjs)
and the ISAM files of -parse_seqids stay in the work directory.
Rerunning changes the dates stored in the index files and in expected/*.info.tsv.
"""
import argparse
import random
import re
import shutil
import subprocess
from pathlib import Path

# The LMDB file of a version 5 database (.pdb / .ndb) must exist for NCBI and for LOSAT.
KEEP = {".pin", ".psq", ".phr", ".nin", ".nsq", ".nhr", ".pal", ".nal", ".pdb", ".ndb"}


# Date strings of a first and a second volume (at most 27 bytes each).
DATE_PAIRS = [
    ("Oct 12, 2026  3:22 PM", "Oct 12, 2026  3:23 PM"),
    ("Oct 12, 2026  3:23 PM", "Oct 12, 2026  3:22 PM"),
    ("Oct 12, 2026  11:59 AM", "Oct 12, 2026  12:01 PM"),
    ("Oct 12, 2026  12:30 AM", "Oct 12, 2026  1:00 AM"),
    ("Oct 12, 2026  1:00 PM", "Oct 12, 2026  12:30 PM"),
    ("Oct 12, 2026  12:30 PM", "Oct 12, 2026  12:30 AM"),
    ("Dec 31, 2025  11:59 PM", "Jan 1, 2026  12:00 AM"),
    ("Mar 9, 2026  9:05 AM", "Feb 28, 2026  11:00 PM"),
    ("oct 12, 2026  3:22 pm", "OCT 12, 2026  3:23 PM"),
    ("Oct 12, 2026  3:22 PM", "Oct 12, 2026 3:23 PM"),
    ("Oct 12, 2026  3:22 PM", "Oct 12, 2026   03:23 PM"),
    ("Oct 12, 2026  3:22 PM", "Oct 12, 2026  3:23 PM "),
    ("Oct 12, 2026  3:22 PM", "Oct 12, 2026  3:23"),
    ("Oct 12, 2026  3:22", "Oct 12, 2026  3:21 PM"),
    ("Oct 12, 2026  3:22 PM", "Oct 12, 2026  13:23 PM"),
    ("Oct 12, 2026  3:22 PM", "Oct 12, 2026  15:23"),
    ("Oct 12, 2026  3:22 PM", "Feb 30, 2027  3:22 PM"),
    ("Oct 12, 2026  3:22 PM", "2027-01-01"),
    ("2027-01-01", "Oct 12, 2026  3:22 PM"),
    ("Oct 12, 2026  3:22 PM", "Oct , 2027  3:22 PM"),
    ("Oct , 2027  3:22 PM", "Oct 12, 2026  3:22 PM"),
    ("Oct 12, 2026  3:22 PM", "Oct 12, 2027"),
    ("Oct 12, 2026  3:22 PM", "Oct 12 2027  3:22 PM"),
    ("Oct 12, 2026  3:22 PM", "October 12, 2027  3:22 PM"),
    ("Oct 12, 2026  3:22 PM", "Oct 12, 27  3:22 PM"),
    ("Oct 12, 2026  3:22 PM", "Oct 12, 1500  3:22 PM"),
    ("", "Oct 12, 2026  3:22 PM"),
    ("Oct 12, 2026  3:22 PM", ""),
    ("not a date", "Oct 12, 2026  3:22 PM"),
]


def write_fasta(path, records, width=60):
    with open(path, "w", newline="\n") as f:
        for title, seq in records:
            f.write(f">{title}\n")
            for i in range(0, len(seq), width):
                f.write(seq[i:i + width] + "\n")


def protein_records():
    rng = random.Random(20261012)
    aa = "ACDEFGHIKLMNPQRSTVWY"

    def rand(n):
        return "".join(rng.choice(aa) for _ in range(n))

    long_title = "p05 a title longer than 127 bytes " + " ".join(f"word{i:02d}" for i in range(20))
    longer_title = "p06 a title longer than 255 bytes " + " ".join(f"token{i:03d}" for i in range(40))
    return [
        ("p01 first protein [Test organism]", rand(61)),
        ("p02", rand(1)),
        ("p03 every letter of NCBIstdaa but the gap", "ABCDEFGHIKLMNPQRSTVWXYZUOJ" + rand(7)),
        ("p04 lower case and ambiguity codes", rand(20).lower() + "BXZ" + rand(20) + "xbz"),
        (long_title, rand(130)),
        (longer_title, rand(88)),
        ("p07  two  spaces  and  punctuation: (a) [b] {c} 100% a=b;c|d", rand(45)),
        ("p08 selenocysteine pyrrolysine", "MU" + rand(30) + "O" + rand(12) + "J"),
        ("p09 stop inside", rand(25) + "*" + rand(25)),
        ("p10 last protein", rand(300)),
    ]


def nucleotide_records():
    rng = random.Random(1012)

    def rand(n):
        return "".join(rng.choice("ACGT") for _ in range(n))

    s5 = rand(30) + "N" * 5 + rand(20) + "R" + rand(3) + "Y" + rand(14)
    s6 = rand(17) + "RYKMSWBDHVN" + rand(9) + "N" * 40 + rand(22)
    s9 = rand(50) + "N" * 16 + rand(3) + "N" * 17 + rand(31)
    long_title = "n11 a title longer than 127 bytes " + " ".join(f"word{i:02d}" for i in range(20))
    return [
        ("n01 length 40, no ambiguity", rand(40)),
        ("n02 length 41", rand(41)),
        ("n03 length 42", rand(42)),
        ("n04 length 43", rand(43)),
        ("n05 short ambiguity runs (old format)", s5),
        ("n06 every ambiguity code and a run of 40 (new format)", s6),
        ("n07 one base", "G"),
        ("n08 all N", "N" * 23),
        ("n09 runs of 16 and 17", s9),
        ("n10 lower case", rand(33).lower() + "nnnn" + rand(10).lower()),
        (long_title, rand(259)),
        ("n12 ambiguity at both ends", "N" + rand(62) + "K"),
    ]


def run(cmd, cwd):
    """Run an NCBI program; its standard output as bytes (blastdbcmd %s is not always UTF-8)."""
    result = subprocess.run(cmd, cwd=cwd, capture_output=True)
    if result.returncode != 0:
        raise SystemExit(f"{' '.join(map(str, cmd))}\n{result.stdout!r}{result.stderr!r}")
    return result.stdout


def patch_date(index_path, new_date):
    """Replace the date string of a version 4 index file, keeping its padded length."""
    data = bytearray(index_path.read_bytes())
    title_len = int.from_bytes(data[8:12], "big")
    at = 12 + title_len
    date_len = int.from_bytes(data[at:at + 4], "big")
    encoded = new_date.encode("ascii")
    assert len(encoded) <= date_len
    data[at + 4:at + 4 + date_len] = encoded + b"\0" * (date_len - len(encoded))
    index_path.write_bytes(bytes(data))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ncbi-bin", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--out", default=str(Path(__file__).resolve().parent))
    args = ap.parse_args()
    bin_dir = Path(args.ncbi_bin)
    work = Path(args.work)
    out = Path(args.out)
    if work.exists():
        raise SystemExit(f"{work} exists; give a new work directory")
    work.mkdir(parents=True)

    write_fasta(work / "prot.faa", protein_records())
    write_fasta(work / "nucl.fna", nucleotide_records())

    def makedb(name, fasta, dbtype, version, *extra):
        run([bin_dir / "makeblastdb", "-in", fasta, "-dbtype", dbtype, "-out", name,
             "-blastdb_version", str(version), *extra], work)

    makedb("prot_v4", "prot.faa", "prot", 4, "-title", "LOSAT seqdb fixture, protein, version 4")
    makedb("prot_v5", "prot.faa", "prot", 5)
    makedb("nucl_v4", "nucl.fna", "nucl", 4, "-title", "LOSAT seqdb fixture, nucleotide, version 4")
    makedb("nucl_v5", "nucl.fna", "nucl", 5)
    makedb("prot_mv", "prot.faa", "prot", 4, "-max_file_sz", "300B", "-title", "protein in volumes")
    makedb("nucl_mv", "nucl.fna", "nucl", 5, "-max_file_sz", "300B", "-title", "nucleotide in volumes")
    makedb("prot_ids_v4", "prot.faa", "prot", 4, "-parse_seqids")
    makedb("prot_ids_v5", "prot.faa", "prot", 5, "-parse_seqids")

    # A copy of prot_v4 with an earlier date: volumes of one database with
    # different dates.
    for ext in (".pin", ".psq", ".phr"):
        shutil.copyfile(work / f"prot_v4{ext}", work / f"prot_old{ext}")
    patch_date(work / "prot_old.pin", "Jan 5, 2020  11:07 PM")

    # Volume order. NCBI sorts volumes by file name without a final .pal/.nal/
    # .pin/.nin/.pdb/.ndb, so "ord.pdb" (compared as "ord") comes before
    # "ord-extra" although '-' sorts before '.'. The two volumes hold different
    # proteins, and the second one is named with a directory.
    records = protein_records()
    write_fasta(work / "ord_a.faa", records[:4])
    write_fasta(work / "ord_b.faa", records[4:7])
    makedb("ord_a", "ord_a.faa", "prot", 4, "-title", "first four proteins")
    makedb("ord_b", "ord_b.faa", "prot", 4, "-title", "next three proteins")
    (work / "sub").mkdir()
    for ext in (".pin", ".psq", ".phr"):
        (work / f"ord_a{ext}").rename(work / f"ord.pdb{ext}")
        (work / f"ord_b{ext}").rename(work / "sub" / f"ord-extra{ext}")
    (work / "alias_order.pal").write_bytes(b"DBLIST sub/ord-extra ord.pdb\n")

    # A volume without sequences (NCBI allows it), first in the volume order.
    index = (work / "prot_v4.pin").read_bytes()
    title_len = int.from_bytes(index[8:12], "big")
    date_end = 12 + title_len + 4 + int.from_bytes(index[12 + title_len:16 + title_len], "big")
    (work / "empty.pin").write_bytes(index[:date_end] + bytes(4) + bytes(8) + bytes(4)
                                     + (0).to_bytes(4, "big") + (1).to_bytes(4, "big"))
    (work / "empty.psq").write_bytes(b"\0")
    (work / "empty.phr").write_bytes(b"")
    (work / "alias_emptyvol.pal").write_bytes(b"TITLE with an empty volume\nDBLIST prot_v4 empty\n")

    run([bin_dir / "blastdb_aliastool", "-dblist", "prot_v4 prot_old", "-dbtype", "prot",
         "-out", "alias_nseq", "-title", "alias with NSEQ and LENGTH"], work)
    (work / "alias_nested.pal").write_bytes(
        b"#\n# no TITLE: the titles below are joined\n#\nDBLIST prot_v4 \"prot_mv\" prot_old\n")
    (work / "alias_titled.pal").write_bytes(
        b"# comment\r\n  TITLE\ta title\twith tabs  \r\n\r\nFUTURE_KEY ignored\r\n"
        b"DBLIST   prot_old\tprot_v4\r\nTITLE\tthe second TITLE line wins \r\n")
    # NSEQ and LENGTH written by hand: the sums over prot_v4 and prot_old are 20 and 1602
    # (alias_nseq, above). One too many, one too few, values that are not plain digits.
    for name, nseq, length in [("alias_nseq_high", b"21", b"1602"), ("alias_length_low", b"20", b"1601"),
                               ("alias_nseq_plus", b"+20", b"1602"), ("alias_nseq_exp", b"2e1", b"1602"),
                               ("alias_length_space", b"20", b"1 602")]:
        (work / f"{name}.pal").write_bytes(
            b"TITLE hand-written sums\nDBLIST prot_v4 prot_old\nNSEQ " + nseq + b"\nLENGTH " + length + b"\n")
    # Nested: alias_nseq (right values, 20 and 1602) and prot_mv (10 and 801) under an outer
    # node with right values (30, 2403), then with wrong ones.
    (work / "alias_sums_nested.pal").write_bytes(
        b"TITLE nested sums\nDBLIST alias_nseq prot_mv\nNSEQ 30\nLENGTH 2403\n")
    (work / "alias_sums_nested_bad.pal").write_bytes(
        b"TITLE nested sums, outer wrong\nDBLIST alias_nseq prot_mv\nNSEQ 30\nLENGTH 2402\n")
    # The inner node is wrong, the outer one has the key too (NCBI never reads the inner value).
    (work / "alias_sums_inner_bad.pal").write_bytes(
        b"TITLE nested sums, inner wrong\nDBLIST alias_nseq_high prot_mv\nNSEQ 31\nLENGTH 2403\n")
    (work / "alias_gilist.pal").write_bytes(b"TITLE filtered\nDBLIST prot_v4\nGILIST some.gil\n")
    (work / "alias_missing.pal").write_bytes(b"TITLE missing volume\nDBLIST prot_v4 no_such_volume\n")
    (work / "alias_loop_a.pal").write_bytes(b"DBLIST alias_loop_b\n")
    (work / "alias_loop_b.pal").write_bytes(b"DBLIST alias_loop_a\n")
    (work / "alias_twice.pal").write_bytes(b"DBLIST prot_v4 alias_titled\n")
    (work / "alias_empty.pal").write_bytes(b"TITLE no DBLIST\n")
    (work / "alias_mixed.pal").write_bytes(b"TITLE versions 4 and 5\nDBLIST prot_v4 prot_v5\n")

    expected = out / "expected"
    db = out / "db"
    for stale in (expected, db):
        if stale.exists():
            shutil.rmtree(stale)
    expected.mkdir(parents=True)
    db.mkdir(parents=True)
    readable = [("prot_v4", "prot"), ("prot_v5", "prot"), ("nucl_v4", "nucl"), ("nucl_v5", "nucl"),
                ("prot_mv", "prot"), ("nucl_mv", "nucl"), ("alias_nested", "prot"),
                ("alias_titled", "prot"), ("prot_old", "prot"), ("alias_order", "prot"),
                ("alias_emptyvol", "prot"), ("alias_nseq", "prot"),
                ("alias_sums_nested", "prot")]
    for name, dbtype in readable:
        info = run([bin_dir / "blastdbcmd", "-db", name, "-dbtype", dbtype, "-info"], work).decode("ascii")
        m = re.search(r"Database: (.*?)\n\t([\d,]+) sequences; ([\d,]+) total (?:residues|bases)\n\n"
                      r"Date: (.*?)\tLongest sequence: ([\d,]+) (?:residues|bases)\n\n"
                      r"BLASTDB Version: (\d+)\n\nVolumes:\n(.*)", info, re.S)
        if not m:
            raise SystemExit(f"unexpected blastdbcmd -info output for {name}:\n{info}")
        volumes = [Path(line.strip()).name for line in m.group(7).splitlines() if line.strip()]
        rows = [("title", m.group(1)), ("sequences", m.group(2).replace(",", "")),
                ("letters", m.group(3).replace(",", "")), ("date", m.group(4)),
                ("longest", m.group(5).replace(",", "")), ("version", m.group(6)),
                ("volumes", " ".join(volumes))]
        (expected / f"{name}.info.tsv").write_text("".join(f"{k}\t{v}\n" for k, v in rows))
        entries = run([bin_dir / "blastdbcmd", "-db", name, "-dbtype", dbtype, "-entry", "all",
                       "-outfmt", "%o\t%l\t%i\t%s\t%t"], work)
        # blastdbcmd %s writes the stop residue (NCBIstdaa 25, '*') as the byte 0xFF.
        (expected / f"{name}.entries.tsv").write_bytes(entries.replace(b"\xff", b"*"))

    # Which date NCBI reports for two volumes with different date strings
    # (CSeqDBImpl::GetDate). The volumes stay in the work directory; the rows
    # are the expectation of the test of LOSAT/src/seqdb/date.rs.
    rows = []
    for n, (first, second) in enumerate(DATE_PAIRS):
        pair = work / f"dates{n:02d}"
        pair.mkdir()
        for vol, date in (("d_a", first), ("d_b", second)):
            for ext in (".pin", ".psq", ".phr"):
                shutil.copyfile(work / f"prot_v4{ext}", pair / f"{vol}{ext}")
            patch_date(pair / f"{vol}.pin", date)
        (pair / "d.pal").write_bytes(b"TITLE two dates\nDBLIST d_a d_b\n")
        result = subprocess.run([bin_dir / "blastdbcmd", "-db", "d", "-dbtype", "prot", "-info"],
                                cwd=pair, capture_output=True)
        m = re.search(rb"Date: (.*?)\tLongest sequence", result.stdout)
        if result.returncode == 0 and m:
            reported = m.group(1).decode("ascii")
        else:
            reported = "ERROR"
        rows.append(f"{first}\t{second}\t{reported}\n")
    (expected / "date_pairs.tsv").write_text("".join(rows))

    (db / "sub").mkdir()
    for path in sorted(work.iterdir()) + sorted((work / "sub").iterdir()):
        if path.is_file() and path.suffix in KEEP:
            shutil.copyfile(path, db / path.relative_to(work))
    shutil.copyfile(work / "prot.faa", out / "prot.faa")
    shutil.copyfile(work / "nucl.fna", out / "nucl.fna")
    total = sum(p.stat().st_size for p in out.rglob("*") if p.is_file())
    print(f"fixtures written to {out}: {total} bytes in {sum(1 for p in out.rglob('*') if p.is_file())} files")


if __name__ == "__main__":
    main()
