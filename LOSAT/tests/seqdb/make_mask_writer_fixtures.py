#!/usr/bin/env python3
"""Build the fixtures of the mask-column writer (LOSAT/tests/seqdb/mask_writer).

`make_mask_fixtures.py` builds the databases that the mask-column reader is tested with. This
script builds what the writer (`writedb::mask_data`) is tested with: small makeblastdb runs, each
with its own mask files, that exercise one rule of NCBI's writer: how Seq-locs of one sequence are
combined, how algorithm ids are assigned, what happens when no sequence matches a mask, when a
range lies past the end of its sequence, when a volume has no masked sequence.

NCBI BLAST+ 2.17.0 is the comparison oracle only: makeblastdb writes the databases (and the
error text), dustmasker / segmasker write `dust.asnt` / `seg.asnt`, and blastdbcmd writes what NCBI
reads back from each database (`expected/<case>.masks.tsv`, `expected/<case>.algorithms.txt`).
The Rust tests (`src/writedb/mask_data/tests.rs`) read the committed files and run no NCBI program;
they feed the same mask files and FASTA files through the writer and compare every file byte for
byte (the dates of the index and column files are read from NCBI's files).

usage: make_mask_writer_fixtures.py --ncbi-bin DIR --work DIR [--out LOSAT/tests/seqdb/mask_writer]

No database here uses -parse_seqids: the masks name the sequences by their number in the FASTA file
(`local id N`); a second -in file starts counting again at 1.
Rerunning changes the dates stored in the index and column files.

Output files:
  <name>.asnt, in_*.fna        the inputs (names as given to makeblastdb, which runs in the work dir)
  cases.tsv                    one row per case: name, dbtype, -in files (space separated), -mask_data
                               files (comma separated), -mask_id (or -), -mask_desc (or -),
                               -max_file_sz (or -), -title
  db/<case>.*                  NCBI's database files of a successful case (.njs/.pjs not kept)
  expected/<case>.rc           makeblastdb's exit status
  expected/<case>.stderr       makeblastdb's standard error
  expected/<case>.masks.tsv    blastdbcmd -entry all -outfmt "%o\t%m<algo>" for each algorithm
  expected/<case>.algorithms.txt   blastdbcmd -info's table of algorithms
"""
import argparse
import re
import shutil
import subprocess
from pathlib import Path

NUCL_LENGTHS = [640, 350, 40, 390, 1, 810, 320, 285]  # mask_nucl.fna
PROT_LENGTHS = [164, 70, 147, 12, 123]  # mask_prot.faa
TITLE = "mask writer fixture"
# The files LOSAT writes (makeblastdb also writes .njs / .pjs, which LOSAT does not).
DB_EXTENSIONS = {"nin", "nsq", "nhr", "naa", "nab", "nac", "nal",
                 "pin", "psq", "phr", "paa", "pab", "pac", "pal"}


def interval(a, b, lid, strand=None, indent=10):
    pad = " " * indent
    st = f"{pad}  strand {strand},\n" if strand else ""
    return f"{pad}{{\n{pad}  from {a},\n{pad}  to {b},\n{st}{pad}  id local id {lid}\n{pad}}}"


def packed(lid, ranges):
    """A packed-int Seq-loc; a range is (from, to) or (from, to, strand)."""
    ints = ",\n".join(interval(r[0], r[1], lid, r[2] if len(r) > 2 else None) for r in ranges)
    return "      packed-int {\n" + ints + "\n      }"


def single(lid, a, b, strand=None):
    """An int Seq-loc."""
    st = f"        strand {strand},\n" if strand else ""
    return f"      int {{\n        from {a},\n        to {b},\n{st}        id local id {lid}\n      }}"


def mask_object(program, options, locs):
    """One Blast-db-mask-info in ASN.1 text; `locs` are strings made by packed() / single()."""
    body = ",\n".join(locs)
    return (f"Blast-db-mask-info ::= {{\n  algo-id 0,\n  algo-program {program},\n"
            f"  algo-options \"{options}\",\n  masks {{\n    masks {{\n{body}\n    }},\n"
            f"    more FALSE\n  }}\n}}\n")


def fasta_records(path):
    """The records of a FASTA file, each from its definition line (without '>') to its last line."""
    return Path(path).read_text().split(">")[1:]


def write_fasta(path, records):
    Path(path).write_text("".join(">" + r for r in records))


def run(cmd, cwd, check=True):
    result = subprocess.run([str(c) for c in cmd], cwd=cwd, capture_output=True)
    if check and result.returncode != 0:
        raise SystemExit(f"{' '.join(map(str, cmd))}\n{result.stdout!r}{result.stderr!r}")
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ncbi-bin", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--out", default=str(Path(__file__).resolve().parent / "mask_writer"))
    args = ap.parse_args()
    bin_dir = Path(args.ncbi_bin)
    work = Path(args.work)
    out = Path(args.out)
    seqdb = Path(__file__).resolve().parent
    if work.exists():
        raise SystemExit(f"{work} exists; give a new work directory")
    work.mkdir(parents=True)
    (out / "db").mkdir(parents=True, exist_ok=True)
    (out / "expected").mkdir(parents=True, exist_ok=True)

    # Inputs: the sequences of the mask fixtures, the first three and the last five nucleotide
    # sequences as two files, and four 40-base sequences for the case where a column file fills up.
    shutil.copyfile(seqdb / "mask_nucl.fna", work / "mask_nucl.fna")
    shutil.copyfile(seqdb / "mask_prot.faa", work / "mask_prot.faa")
    records = fasta_records(work / "mask_nucl.fna")
    write_fasta(work / "in_a.fna", records[:3])
    write_fasta(work / "in_b.fna", records[3:])
    write_fasta(work / "in_col.fna", [f"c{i} forty bases\n" + "ACGT" * 10 + "\n" for i in range(1, 5)])
    # 105 sequences of 20 bases: with -max_file_sz 100B each is a volume, and more than 100 volumes
    # are renamed with three digits at the end.
    write_fasta(work / "in_many.fna", [f"q{i} short\n" + "ACGT" * 5 + "\n" for i in range(1, 106)])
    run([bin_dir / "dustmasker", "-in", "mask_nucl.fna", "-outfmt", "maskinfo_asn1_text", "-out",
         "dust.asnt"], work)
    run([bin_dir / "segmasker", "-in", "mask_prot.faa", "-outfmt", "maskinfo_asn1_text", "-out",
         "seg.asnt"], work)
    # The hand-written masks of make_mask_fixtures.py (same ranges), as wm.asnt / wm_opts.asnt.
    for name in ("wm", "wm_opts", "wm_colon"):
        shutil.copyfile(seqdb / f"mask_{name}.asnt", work / f"{name}.asnt")

    def write(name, *objects):
        (work / name).write_text("".join(objects))

    wm = "windowmasker"
    # Several Seq-locs for one sequence: overlapping, abutting, contained, unsorted (sequence 1);
    # an int among packed-ints, one inside another (2); single-base abutting ranges (4); a
    # range list in the middle of another (6).
    write("merge.asnt", mask_object(wm, "", [
        packed(1, [(300, 309), (0, 9)]), packed(1, [(5, 20), (310, 320), (100, 100)]),
        packed(2, [(50, 60)]), packed(2, [(0, 100)]), single(2, 20, 30), single(2, 101, 101),
        single(2, 300, 310), single(2, 305, 305),
        single(4, 0, 0), single(4, 1, 1), single(4, 3, 3), single(4, 10, 12), single(4, 13, 20),
        packed(6, [(800, 809), (700, 709)]), single(6, 705, 720),
    ]))
    # A Seq-loc alone is kept as given: unsorted, overlapping, even a repeated range.
    write("alone.asnt", mask_object(wm, "", [
        packed(1, [(300, 309), (0, 9), (5, 20), (8, 8)]), single(2, 7, 7),
        packed(4, [(100, 150), (100, 150)]),
    ]))
    # The same sequence in two objects of one file: the second adds to the first, and its options
    # are not read (the algorithm is registered by the first object).
    write("objs.asnt",
          mask_object(wm, "", [packed(1, [(300, 309), (0, 9)]), single(2, 3, 4)]),
          mask_object(wm, "-ignored 2", [packed(1, [(5, 20), (310, 320)]), single(3, 0, 39)]))
    # Strands: the plus-like and the reverse ranges are merged apart (plus first, then reverse by
    # decreasing end); one Seq-loc alone stays as given.
    write("strand.asnt", mask_object(wm, "", [
        packed(1, [(300, 309, "minus"), (0, 9, "plus")]),
        packed(1, [(5, 20, "minus"), (310, 320, "plus"), (7, 8, "minus")]),
        single(2, 0, 9, "both"), single(2, 10, 19, "plus"), single(2, 20, 29, "unknown"),
        single(2, 40, 49, "both-rev"), single(2, 30, 39, "both-rev"),
        packed(4, [(300, 309, "minus"), (0, 9, "minus")]),
    ]))
    # Masks of sequences the input does not have.
    write("nomatch.asnt", mask_object(wm, "", [packed(99, [(0, 9)]), single(100, 0, 0)]))
    # A range one past the end of sequence 2 (350 bases), after a good mask of sequence 1.
    write("oob_end.asnt", mask_object(wm, "", [packed(1, [(0, 9)]), packed(2, [(349, 350)])]))
    # A bad range in sequence 3 (40 bases) after good masks of 1 and 2; and the last base.
    write("oob_later.asnt", mask_object(wm, "", [
        packed(1, [(0, 9)]), packed(2, [(340, 349)]), packed(3, [(39, 40)])]))
    write("last_base.asnt", mask_object(wm, "", [
        packed(1, [(0, 9)]), packed(2, [(340, 349)]), packed(3, [(39, 39)])]))
    write("oob_reversed.asnt", mask_object(wm, "", [packed(1, [(0, 9)]), single(2, 10, 5)]))
    write("empty_range.asnt", mask_object(wm, "", [single(1, 6, 5), packed(2, [(6, 5), (20, 25)])]))
    # Algorithm ids.
    # Eleven dust files with different options, each masking a sequence of at least 285 bases.
    long_sequences = [1, 2, 4, 6, 7, 8]
    for i in range(11):
        write(f"dd{i}.asnt", mask_object("dust", f"-x {i}" if i else "",
                                         [packed(long_sequences[i % 6], [(i, i + 5)])]))
    write("seg_d.asnt", mask_object("seg", "", [packed(2, [(1, 4)])]))
    write("seg_o.asnt", mask_object("seg", "-window 5", [packed(2, [(2, 5)]), packed(4, [(0, 11)])]))
    write("wm_d.asnt", mask_object(wm, "", [packed(3, [(2, 9)])]))
    write("wm_d2.asnt", mask_object(wm, "", [packed(5, [(0, 0)])]))
    write("wm_o1.asnt", mask_object(wm, "-x 1", [packed(6, [(10, 29)])]))
    write("wm_o2.asnt", mask_object(wm, "-x 2", [packed(7, [(100, 199)])]))
    write("rep_d.asnt", mask_object("repeat", "", [packed(8, [(0, 99)])]))
    write("rep_o.asnt", mask_object("repeat", "-x 1", [packed(1, [(600, 639)])]))
    write("oth_d.asnt", mask_object("other", "", [packed(2, [(0, 49)])]))
    write("oth_o.asnt", mask_object("other", "-x:1", [packed(4, [(0, 49)])]))
    write("p55.asnt", mask_object("55", "", [packed(1, [(0, 9)])]))
    write("pnotset.asnt", mask_object("not-set", "", [packed(1, [(0, 9)])]))
    write("pmax.asnt", mask_object("max", "", [packed(1, [(0, 9)])]))
    # A file without masks: the algorithm is registered, no sequence gets mask data.
    write("none.asnt", mask_object(wm, "-e", []))
    # Two input files: the masks name sequences 1 and 2 of each.
    write("two_in.asnt", mask_object(wm, "", [packed(1, [(10, 19)]), packed(2, [(0, 0)])]))
    # One volume without a masked sequence (-max_file_sz 400B makes volumes of 4 / 4 sequences).
    write("vol_first.asnt", mask_object(wm, "", [packed(1, [(0, 9)]), packed(2, [(5, 5)])]))
    write("vol_last.asnt", mask_object(wm, "", [packed(8, [(0, 9)])]))
    write("vol_mid.asnt", mask_object(wm, "", [packed(3, [(0, 9)]), packed(6, [(1, 1)])]))
    # A column file that fills up before the sequence files do: 20 ranges of 8 bytes in a blob.
    write("col_roll.asnt", mask_object(wm, "", [
        packed(i, [(2 * j, 2 * j) for j in range(20)]) for i in range(1, 5)]))

    write("many.asnt", mask_object(wm, "", [
        packed(1, [(0, 3)]), packed(50, [(2, 9)]), packed(101, [(1, 1)]), packed(105, [(0, 19)])]))

    cases = [
        # name, dbtype, inputs, masks, ids, descs, max file size
        ("merge", "nucl", "mask_nucl.fna", "merge.asnt", "-", "-", "-"),
        ("merge_small", "nucl", "mask_nucl.fna", "merge.asnt", "-", "-", "150B"),
        ("alone", "nucl", "mask_nucl.fna", "alone.asnt", "-", "-", "-"),
        ("objs", "nucl", "mask_nucl.fna", "objs.asnt", "-", "-", "-"),
        ("strand", "nucl", "mask_nucl.fna", "strand.asnt", "-", "-", "-"),
        ("nomatch", "nucl", "mask_nucl.fna", "nomatch.asnt", "-", "-", "-"),
        ("oob_end", "nucl", "mask_nucl.fna", "oob_end.asnt", "-", "-", "-"),
        ("oob_later", "nucl", "mask_nucl.fna", "oob_later.asnt", "-", "-", "-"),
        ("oob_multi", "nucl", "mask_nucl.fna", "oob_end.asnt", "-", "-", "400B"),
        ("last_base", "nucl", "mask_nucl.fna", "last_base.asnt", "-", "-", "-"),
        ("oob_reversed", "nucl", "mask_nucl.fna", "oob_reversed.asnt", "-", "-", "-"),
        ("empty_range", "nucl", "mask_nucl.fna", "empty_range.asnt", "-", "-", "-"),
        ("ids_dust", "nucl", "mask_nucl.fna", "dd0.asnt,dd1.asnt,dd2.asnt", "-", "-", "-"),
        ("ids_dust_nine", "nucl", "mask_nucl.fna",
         ",".join(f"dd{i}.asnt" for i in range(10)), "-", "-", "-"),
        ("ids_dust_toomany", "nucl", "mask_nucl.fna",
         ",".join(f"dd{i}.asnt" for i in range(11)), "-", "-", "-"),
        ("ids_mix", "nucl", "mask_nucl.fna",
         "wm_o1.asnt,dd0.asnt,seg_d.asnt,rep_d.asnt,rep_o.asnt,oth_d.asnt,oth_o.asnt,wm_d.asnt,"
         "wm_o2.asnt", "-", "-", "-"),
        ("ids_user", "nucl", "mask_nucl.fna", "wm_d.asnt,dd0.asnt,seg_d.asnt",
         "alpha,beta:2,gamma", "first,second: desc,third", "-"),
        ("dup", "nucl", "mask_nucl.fna", "wm_d.asnt,wm_d.asnt", "-", "-", "-"),
        ("dup_id", "nucl", "mask_nucl.fna", "wm_d.asnt,dd0.asnt", "a,a", "-", "-"),
        ("dup_default", "nucl", "mask_nucl.fna", "wm_d.asnt,wm_d2.asnt", "-", "-", "-"),
        ("invalid55", "nucl", "mask_nucl.fna", "p55.asnt", "-", "-", "-"),
        ("invalid_notset", "nucl", "mask_nucl.fna", "pnotset.asnt", "-", "-", "-"),
        ("invalid_max", "nucl", "mask_nucl.fna", "pmax.asnt", "-", "-", "-"),
        ("none", "nucl", "mask_nucl.fna", "none.asnt", "-", "-", "-"),
        ("none_then_masks", "nucl", "mask_nucl.fna", "none.asnt,wm_d.asnt", "-", "-", "-"),
        ("two_in", "nucl", "in_a.fna in_b.fna", "two_in.asnt", "-", "-", "-"),
        ("vol_first", "nucl", "mask_nucl.fna", "vol_first.asnt", "-", "-", "400B"),
        ("vol_last", "nucl", "mask_nucl.fna", "vol_last.asnt", "-", "-", "400B"),
        ("vol_mid", "nucl", "mask_nucl.fna", "vol_mid.asnt", "-", "-", "400B"),
        ("col_roll", "nucl", "in_col.fna", "col_roll.asnt", "-", "-", "400B"),
        ("many_volumes", "nucl", "in_many.fna", "many.asnt", "-", "-", "100B"),
        ("prot_two", "prot", "mask_prot.faa", "seg.asnt,seg_o.asnt", "-", "-", "-"),
        ("prot_ids", "prot", "mask_prot.faa", "seg.asnt,seg_o.asnt", "s1,s2", "x,y", "-"),
        ("nucl_dust_wm", "nucl", "mask_nucl.fna", "dust.asnt,wm.asnt", "-", "-", "-"),
        ("plain", "nucl", "mask_nucl.fna", "-", "-", "-", "-"),
    ]

    rows = []
    for name, dbtype, inputs, masks, ids, descs, max_size in cases:
        cmd = [bin_dir / "makeblastdb", "-in", inputs, "-dbtype", dbtype, "-out", name,
               "-blastdb_version", "4", "-title", TITLE]
        if masks != "-":
            cmd += ["-mask_data", masks]
        if ids != "-":
            cmd += ["-mask_id", ids]
        if descs != "-":
            cmd += ["-mask_desc", descs]
        if max_size != "-":
            cmd += ["-max_file_sz", max_size]
        # -in with several files is one space separated argument.
        result = run(cmd, work, check=False)
        (out / "expected" / f"{name}.rc").write_text(f"{result.returncode}\n")
        (out / "expected" / f"{name}.stderr").write_bytes(result.stderr)
        rows.append("\t".join([name, dbtype, inputs, masks, ids, descs, max_size, TITLE]) + "\n")
        files = sorted(p for p in work.iterdir() if p.name.startswith(name + ".") and
                       p.suffix[1:] in DB_EXTENSIONS)
        for path in files:
            shutil.copyfile(path, out / "db" / path.name)
        if result.returncode == 0:
            info = run([bin_dir / "blastdbcmd", "-db", name, "-dbtype", dbtype, "-info"], work).stdout.decode()
            m = re.search(r"\nAvailable filtering algorithms applied to database sequences:\n\n.*?\n(?=\nVolumes:)",
                          info, re.S)
            if m:
                (out / "expected" / f"{name}.algorithms.txt").write_text(m.group(0))
                ids_found = [int(line.split()[0]) for line in m.group(0).splitlines()[4:] if line.strip()]
            else:
                ids_found = []
            tsv = ""
            for algo in ids_found:
                masked = run([bin_dir / "blastdbcmd", "-db", name, "-dbtype", dbtype, "-entry", "all",
                              "-outfmt", f"%o\t%m{algo}"], work).stdout.decode()
                tsv += "".join(f"{algo}\t{line}\n" for line in masked.splitlines())
            (out / "expected" / f"{name}.masks.tsv").write_text(tsv)

    (out / "cases.tsv").write_text("".join(rows))
    for path in sorted(work.iterdir()):
        if path.suffix in (".asnt", ".fna", ".faa") and not path.name.startswith("mask_"):
            shutil.copyfile(path, out / path.name)


if __name__ == "__main__":
    main()
