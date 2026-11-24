#!/usr/bin/env python3
"""Build the fixtures of the GI-mask readers (LOSAT/tests/seqdb/gi_mask_read).

A database that `makeblastdb -gi_mask` wrote keeps its masks in GI-indexed files named by the
`MASKLIST` line of its alias file. This script builds such databases with NCBI BLAST+ 2.17.0's
makeblastdb (the inputs of `make_gi_mask_fixtures.py`, the writer's fixtures) and freezes what NCBI's
readers print on them: blastdbcmd (`-info`, `-outfmt %m`, `-mask_sequence_with`), the search
programs (`-db_soft_mask` / `-db_hard_mask`, and plain searches of every program family), makeblastdb
`-input_type blastdb` and makembindex (`-db_mask`, `-show_filters`). NCBI is the comparison oracle
only: the Rust test (`LOSAT/tests/gi_mask_read.rs`) reads the committed files and runs no NCBI
program.

usage: make_gi_mask_read_fixtures.py --ncbi-bin DIR --work DIR [--out LOSAT/tests/seqdb/gi_mask_read]

The work directory is the root `<R>` of the runs: each database is built in `<R>/<db>` with its
inputs copied in first (and removed after the build), and every read runs in `<R>` or below it, so
that names on the command line are relative as a user types them. `{R}` in a word is the root.

The single-volume databases (`sv5`, `sv4`, `svp`, `big`): NCBI's makeblastdb writes an alias whose
`DBLIST` is the absolute path of the database itself, which no NCBI program can open ("DB alias files
are mutually recursive."). LOSAT's makeblastdb writes the database's file name there instead (the
Owner's decision of 2026-11-20, 「Q6=B」; PD-LOSAT-NCBI-DEFECTS). `losat_single_volume_alias` makes
that one change to NCBI's alias, so that NCBI's readers on the result are the oracle of LOSAT's
readers on LOSAT's single-volume databases. Everything else is NCBI's bytes.

The database `wide` (review S3 of M4c-2, 2026-11-24) has GI masks of 3,000, 5,000 and 3,000 ranges
(`write_wide_inputs`): its data file holds more than 64 KiB, and the third record crosses byte
65,536, where LOSAT's readers of the files change chunk and window. The variants `sv4n*` (review S4
and N12) are copies of `sv4` whose `.nog` (the GI of each sequence, which NCBI reads through its
mapping without a check) is cut or edited (`VARIANTS`).

Output files:
  in/                 the queries `q.fa` (pieces of nucl.fa) and `pq.fa` (two records of prot.fa); the
                      databases' inputs are those of the writer's fixtures (`gi_mask/in/`, the same
                      `write_inputs` and `write_masks`) and of `write_wide_inputs`
  dbs.tsv             one row per database: name, what is done after the build (`-`, `losat-sv`:
                      `losat_single_volume_alias`, `alias`: also an alias file `al` over it made by
                      blastdb_aliastool), setup (as in make_gi_mask_fixtures.py), the words after
                      `makeblastdb`
  variants.tsv        one row per variant: name, the database it copies, the change of its `.nog`
  db/<db>/            every file of the database (the inputs removed); in the text files (alias and
                      metadata) the root is `<R>`
  cases.tsv           one row per read: name, databases it needs (comma separated), the directory
                      it runs in (relative to the root), environment (`-` or `NAME=VALUE`), then the
                      program and its words, one per field; the rows run in order in one root, so
                      a `makeblastdb -input_type blastdb` row's copy is there for the rows after it
  expected/<case>.rc, .stdout, .stderr
                      NCBI's exit status, standard output and error, the root as `<R>`; for the
                      cases of `DEVIATION_5` the standard output is the approved difference 5 of
                      PD-LOSAT-BLASTDBCMD-DEVIATIONS (a mask on a protein in FASTA without a
                      range): NCBI's headers with the letters of NCBI's `-outfmt %s` run of the
                      same words, NCBI's own FASTA kept as `<case>.ncbi_stdout`
  expected/<case>.files
                      for the makembindex rows: each file the run made or changed under the root,
                      `path<TAB>size<TAB>FNV-1a 64 digest` (16 hex digits), in path order
"""
import argparse
import os
import random
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import make_gi_mask_fixtures as writer  # noqa: E402

N = writer.N
NM = writer.NM
P = writer.P
DN = writer.DN
DN2 = writer.DN2
V4 = ["-blastdb_version", "4"]

# name, after the build, setup, words after `makeblastdb`
DBS = [
    ("m01", "alias", DN2, NM + ["-mask_data", "dust_nucl.asnb,dust10_nucl.asnb", "-gi_mask",
                                "-gi_mask_name", "m1,m2"]),
    ("m10", "-", DN2, NM + V4 + ["-mask_data", "dust_nucl.asnb,dust10_nucl.asnb", "-gi_mask",
                                  "-gi_mask_name", "m1,m2"]),
    ("m09", "-", "prot.fa seg_prot.asnb",
     P + ["-max_file_sz", "200", "-mask_data", "seg_prot.asnb", "-gi_mask", "-gi_mask_name", "pm1"]),
    ("m03", "-", DN, NM + V4 + ["-gi_mask"]),
    ("m06", "-", "nucl.fa conv_other.asnb",
     NM + V4 + ["-mask_data", "conv_other.asnb", "-gi_mask", "-gi_mask_name", "oth"]),
    ("m06n", "-", "nucl.fa conv_other.asnb", NM + V4 + ["-mask_data", "conv_other.asnb", "-gi_mask"]),
    ("m07", "-", "gimax.fa dust_gimax.asnb",
     ["-in", "gimax.fa", "-dbtype", "nucl", "-parse_seqids", "-out", "db", "-title", "t"] + V4
     + ["-max_file_sz", "100", "-mask_data", "dust_gimax.asnb", "-gi_mask", "-gi_mask_name", "gx"]),
    ("g01", "-", "dupgi.fa acc_dup.asnt",
     ["-in", "dupgi.fa", "-dbtype", "nucl", "-parse_seqids", "-out", "db", "-title", "t"] + V4
     + ["-max_file_sz", "200", "-mask_data", "acc_dup.asnt", "-gi_mask", "-gi_mask_name", "ad"]),
    ("g02", "-", DN, NM + V4 + ["-mask_data", "dust_nucl.asnb", "-gi_mask", "-gi_mask_name",
                                "{R}/g02/abs_gm"]),
    ("d10", "-", "nucl.fa dust10_nucl.asnb",
     N + V4 + ["-max_file_sz", "10", "-mask_data", "dust10_nucl.asnb", "-gi_mask", "-gi_mask_name",
               "m"]),
    ("d20", "-", "nucl.fa dust10_nucl.asnb",
     N + V4 + ["-max_file_sz", "20", "-mask_data", "dust10_nucl.asnb", "-gi_mask", "-gi_mask_name",
               "m"]),
    ("sp", "-", DN, NM + V4 + ["-mask_data", "dust_nucl.asnb", "-gi_mask", "-gi_mask_name", "a b"]),
    ("sv5", "losat-sv", DN, N + ["-mask_data", "dust_nucl.asnb", "-gi_mask", "-gi_mask_name", "gm1"]),
    ("sv4", "losat-sv", DN, N + V4 + ["-mask_data", "dust_nucl.asnb", "-gi_mask", "-gi_mask_name",
                                      "gm1"]),
    ("svp", "losat-sv", "prot.fa seg_prot.asnb",
     P + ["-mask_data", "seg_prot.asnb", "-gi_mask", "-gi_mask_name", "pm1"]),
    ("big", "losat-sv", "big.fa dust_big.asnb",
     ["-in", "big.fa", "-dbtype", "nucl", "-parse_seqids", "-out", "db", "-title", "t"] + V4
     + ["-mask_data", "dust_big.asnb", "-gi_mask", "-gi_mask_name", "big"]),
    ("wide", "losat-sv", "wide.fa wide.asnt",
     ["-in", "wide.fa", "-dbtype", "nucl", "-parse_seqids", "-out", "db", "-title", "t"] + V4
     + ["-mask_data", "wide.asnt", "-gi_mask", "-gi_mask_name", "wd"]),
]

# name, the database copied, the change of its `db.nog` (`edit_nog`): sv4 has 12 sequences, so
# its `.nog` is a header of 32 bytes and 12 entries of 4 bytes.
VARIANTS = [
    ("sv4n12", "sv4", "cut 12"),      # inside the header: the count reads 0
    ("sv4n15", "sv4", "cut 15"),      # inside the count: its last byte reads 0
    ("sv4n42", "sv4", "cut 42"),      # OIDs 0 and 1 whole, two bytes of OID 2's, then zeros
    ("sv4nx", "sv4", "remove"),       # no .nog: the definition lines' GIs
    ("sv4nf", "sv4", "count3swap"),   # a count of 3, the entries of OIDs 0 and 1 swapped
]

FMT = "%o %a %g %m"

# The cases whose standard output is approved difference 5 of PD-LOSAT-BLASTDBCMD-DEVIATIONS:
# NCBI's FASTA of a protein with a mask and no range keeps the masks only in the first third of
# the sequence; LOSAT masks the whole length, as NCBI's `-outfmt %s` of the same words does.
DEVIATION_5 = {"bdc_m09_msw_0"}
# blastdbcmd's FASTA line length (`-line_length` default).
LINE_LENGTH = 80


def with_letters_of(fasta, letters):
    """NCBI's FASTA `fasta` with the sequence lines of each record replaced by the line of
    `letters` (the `-outfmt %s` run of the same records, one per line), wrapped at
    `LINE_LENGTH`."""
    rows = letters.decode("latin-1").split("\n")
    out = []
    record = -1
    for line in fasta.decode("latin-1").split("\n"):
        if line.startswith(">"):
            record += 1
            out.append(line)
            seq = rows[record]
            out.extend(seq[i:i + LINE_LENGTH] for i in range(0, len(seq), LINE_LENGTH))
        elif record < 0 or line == "":
            out.append(line)
    if record + 1 != len([r for r in rows if r]):
        raise SystemExit("deviation 5: the record counts of the two runs differ")
    return "\n".join(out).encode("latin-1")


def bdc(*words):
    return ["blastdbcmd"] + list(words)


def mbi(*words):
    return ["makembindex"] + list(words)


def bn(*words):
    return ["blastn", "-query", "{R}/q.fa", "-dust", "no"] + list(words)


def cases():
    """name, databases, directory, environment, argv."""
    rows = []

    def add(name, dbs, cwd, env, argv):
        rows.append((name, dbs, cwd, env, argv))

    # blastdbcmd on the multi-volume database of two algorithms (version 5).
    add("bdc_m01_info", "m01", "m01", "-", bdc("-db", "db", "-info"))
    for algo in ("0", "1", "2", "11"):
        add(f"bdc_m01_m{algo}", "m01", "m01", "-", bdc("-db", "db", "-entry", "all", "-outfmt", FMT + algo))
    for algo in ("0", "1", "m1", "m2", "2", "11", "dust", "-1"):
        add(f"bdc_m01_msw_{algo}", "m01", "m01", "-",
            bdc("-db", "db", "-entry", "AB000002.1,AB000008.1,AB000004.1", "-mask_sequence_with", algo))
    # Two algorithms of one database in one run.
    add("bdc_m01_msw1_m0", "m01", "m01", "-",
        bdc("-db", "db", "-entry", "all", "-mask_sequence_with", "1", "-outfmt", "%a %s %m0"))
    add("bdc_m01_msw0_m1", "m01", "m01", "-",
        bdc("-db", "db", "-entry", "all", "-mask_sequence_with", "0", "-outfmt", "%a %s %m1"))
    # Where the GI-mask files are looked for: the current directory and BLASTDB.
    add("bdc_m01_parent_m0", "m01", ".", "-", bdc("-db", "m01/db", "-entry", "all", "-outfmt", FMT + "0"))
    add("bdc_m01_parent_info", "m01", ".", "-", bdc("-db", "m01/db", "-info"))
    add("bdc_m01_parent_msw", "m01", ".", "-",
        bdc("-db", "m01/db", "-entry", "AB000002.1", "-mask_sequence_with", "m1"))
    add("bdc_m01_abs_m0", "m01", ".", "-", bdc("-db", "{R}/m01/db", "-entry", "all", "-outfmt", FMT + "0"))
    add("bdc_m01_blastdb_m0", "m01", ".", "BLASTDB={R}/m01",
        bdc("-db", "db", "-entry", "all", "-outfmt", FMT + "0"))
    add("bdc_m01_blastdb_info", "m01", ".", "BLASTDB={R}/m01", bdc("-db", "db", "-info"))
    # An alias over it, and the database twice.
    add("bdc_m01_al_info", "m01", "m01", "-", bdc("-db", "al", "-info"))
    add("bdc_m01_al_m0", "m01", "m01", "-", bdc("-db", "al", "-entry", "all", "-outfmt", FMT + "0"))
    add("bdc_m01_two_info", "m01", "m01", "-", bdc("-db", "db db", "-info"))
    add("bdc_m01_two_m0", "m01", "m01", "-", bdc("-db", "db db", "-entry", "all", "-outfmt", FMT + "0"))
    # Other databases.
    add("bdc_m10_info", "m10", "m10", "-", bdc("-db", "db", "-info"))
    add("bdc_m10_m1", "m10", "m10", "-", bdc("-db", "db", "-entry", "all", "-outfmt", FMT + "1"))
    add("bdc_m09_info", "m09", "m09", "-", bdc("-db", "db", "-info"))
    add("bdc_m09_m0", "m09", "m09", "-", bdc("-db", "db", "-entry", "all", "-outfmt", FMT + "0"))
    add("bdc_m09_msw_0", "m09", "m09", "-", bdc("-db", "db", "-entry", "all", "-mask_sequence_with", "0"))
    add("bdc_m03_info", "m03", "m03", "-", bdc("-db", "db", "-info"))
    add("bdc_m03_msw_0", "m03", "m03", "-",
        bdc("-db", "db", "-entry", "AB000002.1", "-mask_sequence_with", "0"))
    add("bdc_m03_m0", "m03", "m03", "-", bdc("-db", "db", "-entry", "all", "-outfmt", FMT + "0"))
    add("bdc_m06_info", "m06", "m06", "-", bdc("-db", "db", "-info"))
    add("bdc_m06_m0", "m06", "m06", "-", bdc("-db", "db", "-entry", "all", "-outfmt", FMT + "0"))
    add("bdc_m06_msw_oth", "m06", "m06", "-",
        bdc("-db", "db", "-entry", "AB000002.1", "-mask_sequence_with", "oth"))
    add("bdc_m06n_info", "m06n", "m06n", "-", bdc("-db", "db", "-info"))
    add("bdc_m07_m0", "m07", "m07", "-", bdc("-db", "db", "-entry", "all", "-outfmt", FMT + "0"))
    add("bdc_g01_m0", "g01", "g01", "-", bdc("-db", "db", "-entry", "all", "-outfmt", FMT + "0"))
    add("bdc_g02_parent_m0", "g02", ".", "-", bdc("-db", "g02/db", "-entry", "all", "-outfmt", FMT + "0"))
    add("bdc_d10_m0", "d10", "d10", "-", bdc("-db", "db", "-entry", "all", "-outfmt", FMT + "0"))
    add("bdc_d10_info", "d10", "d10", "-", bdc("-db", "db", "-info"))
    add("bdc_d20_m0", "d20", "d20", "-", bdc("-db", "db", "-entry", "all", "-outfmt", FMT + "0"))
    add("bdc_sp_info", "sp", "sp", "-", bdc("-db", "db", "-info"))
    add("bdc_sp_m1", "sp", "sp", "-", bdc("-db", "db", "-entry", "all", "-outfmt", FMT + "1"))
    add("bdc_sp_msw_b", "sp", "sp", "-", bdc("-db", "db", "-entry", "AB000002.1", "-mask_sequence_with", "b"))
    add("bdc_big_m0", "big", "big", "-", bdc("-db", "db", "-entry", "all", "-outfmt", FMT + "0"))
    for db in ("sv5", "sv4"):
        add(f"bdc_{db}_info", db, db, "-", bdc("-db", "db", "-info"))
        add(f"bdc_{db}_m0", db, db, "-", bdc("-db", "db", "-entry", "all", "-outfmt", FMT + "0"))
        add(f"bdc_{db}_msw_gm1", db, db, "-",
            bdc("-db", "db", "-entry", "AB000002.1,AB000008.1", "-mask_sequence_with", "gm1"))
    add("bdc_svp_m0", "svp", "svp", "-", bdc("-db", "db", "-entry", "all", "-outfmt", FMT + "0"))
    add("bdc_svp_info", "svp", "svp", "-", bdc("-db", "db", "-info"))

    # blastn with the database's masks.
    for kind in ("soft", "hard"):
        for algo in ("0", "1", "m1"):
            for fmt in ("6", "0"):
                add(f"bn_m01_{kind}{algo}_f{fmt}", "m01", "m01", "-",
                    bn("-db", "db", f"-db_{kind}_mask", algo, "-outfmt", fmt))
        for algo in ("2", "11", "nosuch", "-1", "-2"):
            add(f"bn_m01_{kind}{algo}_f6", "m01", "m01", "-",
                bn("-db", "db", f"-db_{kind}_mask", algo, "-outfmt", "6"))
    add("bn_m01_softdust_f6", "m01", "m01", "-", bn("-db", "db", "-db_soft_mask", "dust", "-outfmt", "6"))
    add("bn_m01_plain_f0", "m01", "m01", "-", bn("-db", "db", "-outfmt", "0"))
    add("bn_m01_parent_soft0", "m01", ".", "-", bn("-db", "m01/db", "-db_soft_mask", "0", "-outfmt", "6"))
    add("bn_m01_parent_plain", "m01", ".", "-", bn("-db", "m01/db", "-outfmt", "6"))
    add("bn_m01_blastdb_soft0", "m01", ".", "BLASTDB={R}/m01",
        bn("-db", "db", "-db_soft_mask", "0", "-outfmt", "6"))
    add("bn_m01_al_soft0", "m01", "m01", "-", bn("-db", "al", "-db_soft_mask", "0", "-outfmt", "6"))
    add("bn_m01_al_plain", "m01", "m01", "-", bn("-db", "al", "-outfmt", "6"))
    add("bn_m01_two_hard0_f6", "m01", "m01", "-", bn("-db", "db db", "-db_hard_mask", "0", "-outfmt", "6"))
    add("bn_m01_two_soft0_f0", "m01", "m01", "-", bn("-db", "db db", "-db_soft_mask", "0", "-outfmt", "0"))
    add("bn_m01_two_soft1_f6", "m01", "m01", "-", bn("-db", "db db", "-db_soft_mask", "1", "-outfmt", "6"))
    add("bn_m01_two_softm1_f6", "m01", "m01", "-",
        bn("-db", "db db", "-db_soft_mask", "m1", "-outfmt", "6"))
    add("bn_m10_hard1_f6", "m10", "m10", "-", bn("-db", "db", "-db_hard_mask", "1", "-outfmt", "6"))
    add("bn_m03_soft0_f6", "m03", "m03", "-", bn("-db", "db", "-db_soft_mask", "0", "-outfmt", "6"))
    add("bn_m06_soft0_f0", "m06", "m06", "-", bn("-db", "db", "-db_soft_mask", "0", "-outfmt", "0"))
    add("bn_m07_soft0_f6", "m07", "m07", "-", bn("-db", "db", "-db_soft_mask", "0", "-outfmt", "6"))
    add("bn_d10_soft0_f6", "d10", "d10", "-", bn("-db", "db", "-db_soft_mask", "0", "-outfmt", "6"))
    add("bn_d20_hard0_f6", "d20", "d20", "-", bn("-db", "db", "-db_hard_mask", "0", "-outfmt", "6"))
    add("bn_sp_soft0_f6", "sp", "sp", "-", bn("-db", "db", "-db_soft_mask", "0", "-outfmt", "6"))
    add("bn_g02_parent_soft0", "g02", ".", "-", bn("-db", "g02/db", "-db_soft_mask", "0", "-outfmt", "6"))
    for db in ("sv5", "sv4"):
        add(f"bn_{db}_soft0_f0", db, db, "-", bn("-db", "db", "-db_soft_mask", "0", "-outfmt", "0"))
        add(f"bn_{db}_hard0_f6", db, db, "-", bn("-db", "db", "-db_hard_mask", "gm1", "-outfmt", "6"))

    # Plain searches of a GI-mask database by each program family.
    add("tbn_m01_plain", "m01", "m01", "-",
        ["tblastn", "-query", "{R}/pq.fa", "-db", "db", "-seg", "no", "-outfmt", "6"])
    add("tbx_m01_plain", "m01", "m01", "-",
        ["tblastx", "-query", "{R}/q.fa", "-db", "db", "-seg", "no", "-outfmt", "6"])
    add("bp_m09_plain", "m09", "m09", "-",
        ["blastp", "-query", "{R}/pq.fa", "-db", "db", "-seg", "no", "-outfmt", "6"])
    add("bx_m09_plain", "m09", "m09", "-",
        ["blastx", "-query", "{R}/q.fa", "-db", "db", "-seg", "no", "-outfmt", "6"])
    add("psi_m09_plain", "m09", "m09", "-",
        ["psiblast", "-query", "{R}/pq.fa", "-db", "db", "-seg", "no", "-outfmt", "6"])
    add("bp_svp_plain", "svp", "svp", "-",
        ["blastp", "-query", "{R}/pq.fa", "-db", "db", "-seg", "no", "-outfmt", "6"])

    # The protein programs' database masks (series NP, M4b) read the GI masks through the same
    # entry (`seqdb::mask`): blastp and blastx on the protein databases, tblastn and tblastx on the
    # nucleotide ones, soft and hard, by id and by name.
    prot = {"blastp": "{R}/pq.fa", "blastx": "{R}/q.fa", "tblastn": "{R}/pq.fa", "tblastx": "{R}/q.fa"}
    short = {"blastp": "bp", "blastx": "bx", "tblastn": "tbn", "tblastx": "tbx"}
    masked = [("blastp", "m09", "pm1"), ("blastp", "svp", "pm1"), ("blastx", "m09", "pm1"),
              ("tblastn", "m01", "m1"), ("tblastn", "sv5", "gm1"), ("tblastx", "m01", "m2"),
              ("tblastx", "sv5", "gm1")]
    for program, db, name in masked:
        for kind in ("soft", "hard"):
            for algo, fmt in (("0", "6"), (name, "0")):
                add(f"{short[program]}_{db}_{kind}{algo}_f{fmt}", db, db, "-",
                    [program, "-query", prot[program], "-db", "db", "-seg", "no",
                     f"-db_{kind}_mask", algo, "-outfmt", fmt])

    # `makeblastdb -input_type blastdb` (series RG, D3c) reads the input's masks through the same
    # entry: a GI-mask database's GI masks become the copy's mask columns, with the ids and names
    # that NCBI's raw database source registers (`dust11`, `dust12`, `seg21`). These rows run in
    # table order in one root: each copy is read back by the rows after it.
    for db, dbtype, ids in (("m01", "nucl", ("11", "12")), ("sv5", "nucl", ("11",)),
                            ("svp", "prot", ("21",))):
        add(f"raw_{db}", db, db, "-",
            ["makeblastdb", "-in", "db", "-input_type", "blastdb", "-dbtype", dbtype, "-out", "copy"])
        add(f"raw_{db}_info", db, db, "-", bdc("-db", "copy", "-dbtype", dbtype, "-info"))
        # One `%m` a row: with two, NCBI takes the first one's algorithm for both and leaves the
        # second one's digits in the text, an approved difference of LOSAT's (item 4 of
        # PD-LOSAT-BLASTDBCMD-DEVIATIONS).
        for i in ids:
            add(f"raw_{db}_m{i}", db, db, "-",
                bdc("-db", "copy", "-dbtype", dbtype, "-entry", "all", "-outfmt", f"%a %m{i}"))
    add("raw_m01_hard12", "m01", "m01", "-", bn("-db", "copy", "-db_hard_mask", "12", "-outfmt", "6"))
    # The GI-mask files are looked for from the current directory, as for every reader.
    add("raw_m01_parent", "m01", ".", "-",
        ["makeblastdb", "-in", "m01/db", "-input_type", "blastdb", "-dbtype", "nucl", "-out", "m01/copy_up"])
    # Two names: MASKLIST is passed over, the copy has no masks.
    add("raw_m01_two", "m01", "m01", "-",
        ["makeblastdb", "-in", "db db", "-input_type", "blastdb", "-dbtype", "nucl", "-title", "t2",
         "-out", "copy2"])
    add("raw_m01_two_info", "m01", "m01", "-", bdc("-db", "copy2", "-info"))

    # Records of a data file past 64 KiB (review S3), and the `.nog` cut or edited (S4, N12).
    add("bdc_wide_m0", "wide", "wide", "-", bdc("-db", "db", "-entry", "all", "-outfmt", FMT + "0"))
    add("bdc_wide_info", "wide", "wide", "-", bdc("-db", "db", "-info"))
    add("bdc_wide_msw0", "wide", "wide", "-",
        bdc("-db", "db", "-entry", "WD000003.1", "-mask_sequence_with", "0"))
    for name, _, _ in VARIANTS:
        add(f"bdc_{name}_m0", name, name, "-", bdc("-db", "db", "-entry", "all", "-outfmt", FMT + "0"))
    add("bn_sv4n42_soft0_f6", "sv4n42", "sv4n42", "-", bn("-db", "db", "-db_soft_mask", "0", "-outfmt", "6"))
    add("bn_sv4nf_hard0_f6", "sv4nf", "sv4nf", "-", bn("-db", "db", "-db_hard_mask", "0", "-outfmt", "6"))

    # makembindex (review S2): a new style index opens each volume of `m01` alone (no MASKLIST
    # there) and the one volume of `sv5` with it; an old style index opens the database; `d10`'s
    # data volume 00 was never written. These rows write index files, listed in `<case>.files`.
    nmer = ["-nmer", "8"]
    old = ["-old_style_index", "true", "-iformat", "blastdb", "-output", "o"]
    add("mbi_m01_id0", "m01", "m01", "-", mbi("-input", "db", "-db_mask", "0", *nmer))
    add("mbi_m01_show", "m01", "m01", "-", mbi("-input", "db", "-show_filters"))
    add("mbi_m01_old_m1", "m01", "m01", "-", mbi(*old, "-input", "db", "-db_mask", "m1", *nmer))
    add("mbi_m01_old_show", "m01", "m01", "-", mbi(*old[:4], "-input", "db", "-show_filters"))
    add("mbi_sv5_id0", "sv5", "sv5", "-", mbi("-input", "db", "-db_mask", "0", *nmer))
    add("mbi_sv5_gm1", "sv5", "sv5", "-", mbi("-input", "db", "-db_mask", "gm1", *nmer))
    add("mbi_sv5_nosuch", "sv5", "sv5", "-", mbi("-input", "db", "-db_mask", "nosuch", *nmer))
    add("mbi_sv5_id7", "sv5", "sv5", "-", mbi("-input", "db", "-db_mask", "7", *nmer))
    add("mbi_sv5_show", "sv5", "sv5", "-", mbi("-input", "db", "-show_filters"))
    add("mbi_sv5_old_gm1", "sv5", "sv5", "-", mbi(*old, "-input", "db", "-db_mask", "gm1", *nmer))
    add("mbi_sv5_parent", "sv5", ".", "-", mbi("-input", "sv5/db", "-db_mask", "0", *nmer))
    add("mbi_d10_old0", "d10", "d10", "-", mbi(*old, "-input", "db", "-db_mask", "0", *nmer))
    return rows


def write_wide_inputs(d):
    """`wide.fa`: three sequences with GIs; `wide.asnt`: masks of 3,000, 5,000 and 3,000 ranges of
    five letters, ten letters apart, so that the GI-mask data file has 88,012 bytes (records of
    24,004, 40,004 and 24,004 bytes, the third from byte 64,008 to 88,012)."""
    rng = random.Random(24)
    counts = (3000, 5000, 3000)
    fasta, ints = [], []
    for k, n in enumerate(counts):
        gi = 8001 + k
        seq = "".join(rng.choice("ACGT") for _ in range(10 * n + 40))
        fasta.append(">gi|%d|gb|WD%06d.1| wide %d\n%s\n" % (gi, k + 1, k + 1, writer.wrap(seq)))
        ints.extend("      int { from %d, to %d, id gi %d }" % (10 * i + 2, 10 * i + 6, gi)
                    for i in range(n))
    (d / "wide.fa").write_text("".join(fasta))
    (d / "wide.asnt").write_text(
        "Blast-db-mask-info ::= {\n  algo-id 2,\n  algo-program dust,\n  algo-options \"wide\",\n"
        "  masks {\n    masks {\n" + ",\n".join(ints) + "\n    },\n    more FALSE\n  }\n}\n")


def edit_nog(path, edit):
    """The change `edit` of a `.nog` (`VARIANTS`)."""
    if edit == "remove":
        path.unlink()
        return
    data = bytearray(path.read_bytes())
    if edit.startswith("cut "):
        data = data[:int(edit[4:])]
    elif edit == "count3swap":
        data[12:16] = (3).to_bytes(4, "big")
        data[32:36], data[36:40] = data[36:40], data[32:36]
    else:
        raise SystemExit(f"unknown .nog change {edit}")
    path.write_bytes(bytes(data))


def fnv1a64(data):
    h = 0xcbf29ce484222325
    for b in data:
        h = ((h ^ b) * 0x100000001b3) & 0xFFFFFFFFFFFFFFFF
    return h


def file_states(root):
    """Each file under `root`: its path relative to it (with `/`) and (size, mtime)."""
    out = {}
    for path in root.rglob("*"):
        if path.is_file():
            st = path.stat()
            out[path.relative_to(root).as_posix()] = (st.st_size, st.st_mtime_ns)
    return out


def made_files(root, before):
    """The `.files` text of the files under `root` made or changed since `before`."""
    lines = []
    for name, state in sorted(file_states(root).items()):
        if before.get(name) != state:
            data = (root / name).read_bytes()
            lines.append(f"{name}\t{len(data)}\t{fnv1a64(data):016x}\n")
    return "".join(lines)


def losat_single_volume_alias(db_dir, extension):
    """The one change from NCBI's files (PD-LOSAT-NCBI-DEFECTS, the Owner's 「Q6=B」 of 2026-11-20):
    the alias of a single-volume database names the volume by its file name (`db`), the form of
    NCBI's reader rule `CSeqDBAliasNode::x_ExpandAliases` (seqdbalias.cpp:764-797), instead of the
    absolute path of the database that NCBI's makeblastdb writes."""
    alias = db_dir / f"db.{extension}al"
    lines = alias.read_bytes().split(b"\n")
    changed = [b"DBLIST db" if line.startswith(b"DBLIST ") else line for line in lines]
    if changed == lines:
        raise SystemExit(f"{alias}: no DBLIST line")
    alias.write_bytes(b"\n".join(changed))


def write_queries(d):
    def records(name):
        out = []
        for rec in (d / name).read_text().split(">")[1:]:
            head, _, body = rec.partition("\n")
            out.append((head, body.replace("\n", "")))
        return out

    nucl = records("nucl.fa")
    prot = records("prot.fa")
    w = writer.wrap
    (d / "q.fa").write_text(">q1\n" + w(nucl[0][1][120:400]) + "\n>q2\n" + w(nucl[1][1]) + "\n>q6\n"
                            + w(nucl[5][1][150:330]) + "\n")
    (d / "pq.fa").write_text(">pq1\n" + w(prot[0][1]) + "\n>pq2\n" + w(prot[1][1]) + "\n")


def build(bin_dir, root, inputs, out, env):
    rows = []
    for name, after, setup_text, words in DBS:
        db_dir = root / name
        db_dir.mkdir(parents=True)
        put = writer.setup(db_dir, setup_text.split(), inputs)
        argv = [w.replace("{R}", str(root)) for w in words]
        result = subprocess.run([str(bin_dir / "makeblastdb")] + argv, cwd=db_dir,
                                capture_output=True, env=env)
        want = 2 if after == "losat-sv" else 0
        if result.returncode != want:
            raise SystemExit(f"{name}: makeblastdb rc {result.returncode}\n{result.stderr!r}")
        for f in put:
            (db_dir / f).unlink()
        extension = "p" if "prot" in words else "n"
        if after == "losat-sv":
            losat_single_volume_alias(db_dir, extension)
        if after == "alias":
            r = subprocess.run([str(bin_dir / "blastdb_aliastool"), "-dblist", "db", "-dbtype", "nucl",
                                "-out", "al", "-title", "al"], cwd=db_dir, capture_output=True, env=env)
            if r.returncode != 0:
                raise SystemExit(f"{name}: blastdb_aliastool rc {r.returncode}")
        rows.append("\t".join([name, after, setup_text] + words))
        keep_database(root, db_dir, out / "db" / name)
    (out / "dbs.tsv").write_text("".join(r + "\n" for r in rows))
    rows = []
    for name, base, edit in VARIANTS:
        shutil.copytree(root / base, root / name)
        edit_nog(root / name / "db.nog", edit)
        rows.append("\t".join([name, base, edit]))
        keep_database(root, root / name, out / "db" / name)
    (out / "variants.tsv").write_text("".join(r + "\n" for r in rows))


def keep_database(root, db_dir, keep):
    """Every file of `db_dir` into `keep`, the root `<R>` in the text files."""
    for path in sorted(db_dir.rglob("*")):
        if path.is_dir():
            continue
        target = keep / path.relative_to(db_dir)
        target.parent.mkdir(parents=True, exist_ok=True)
        if path.suffix in (".nal", ".pal", ".njs", ".pjs"):
            text = path.read_bytes().decode("latin-1").replace(str(root), "<R>")
            target.write_bytes(text.encode("latin-1"))
        else:
            shutil.copyfile(path, target)


def read(bin_dir, root, out, env):
    rows = []
    exp = out / "expected"
    for name, dbs, cwd, run_env, argv in cases():
        words = [w.replace("{R}", str(root)) for w in argv]
        e = dict(env)
        if run_env != "-":
            k, _, v = run_env.partition("=")
            e[k] = v.replace("{R}", str(root))
        before = file_states(root) if words[0] == "makembindex" else None
        result = subprocess.run([str(bin_dir / words[0])] + words[1:], cwd=root / cwd,
                                capture_output=True, env=e)
        rows.append("\t".join([name, dbs, cwd, run_env] + argv))
        (exp / f"{name}.rc").write_text(f"{result.returncode}\n")
        if before is not None:
            (exp / f"{name}.files").write_text(made_files(root, before))
        stdout = result.stdout
        if name in DEVIATION_5:
            letters = subprocess.run([str(bin_dir / words[0])] + words[1:] + ["-outfmt", "%s"],
                                     cwd=root / cwd, capture_output=True, env=e, check=True)
            (exp / f"{name}.ncbi_stdout").write_bytes(stdout)
            stdout = with_letters_of(stdout, letters.stdout)
        for stream, data in (("stdout", stdout), ("stderr", result.stderr)):
            text = data.decode("latin-1").replace(str(root), "<R>")
            (exp / f"{name}.{stream}").write_bytes(text.encode("latin-1"))
    (out / "cases.tsv").write_text("".join(r + "\n" for r in rows))
    return len(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ncbi-bin", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--out", default=str(Path(__file__).resolve().parent / "gi_mask_read"))
    args = ap.parse_args()
    bin_dir = Path(args.ncbi_bin).resolve()
    root = Path(args.work).resolve()
    out = Path(args.out)
    if root.exists():
        raise SystemExit(f"{root} exists; give a new work directory")
    inputs = root / "in"
    inputs.mkdir(parents=True)
    writer.write_inputs(inputs)
    writer.write_masks(bin_dir, inputs)
    write_wide_inputs(inputs)
    write_queries(inputs)
    if out.exists():
        shutil.rmtree(out)
    for sub in ("in", "expected", "db"):
        (out / sub).mkdir(parents=True)
    for f in ("q.fa", "pq.fa"):
        shutil.copyfile(inputs / f, out / "in" / f)
    shutil.copyfile(inputs / "q.fa", root / "q.fa")
    shutil.copyfile(inputs / "pq.fa", root / "pq.fa")
    env = dict(os.environ, BLAST_USAGE_REPORT="false")
    for name in ("BLASTDB", "NCBI_CONFIG__BLAST__BLASTDB"):
        env.pop(name, None)
    build(bin_dir, root, inputs, out, env)
    print(f"{read(bin_dir, root, out, env)} cases")


if __name__ == "__main__":
    main()
