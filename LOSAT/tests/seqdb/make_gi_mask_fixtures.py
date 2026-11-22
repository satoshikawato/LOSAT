#!/usr/bin/env python3
"""Build the fixtures of the GI-mask writer (LOSAT/tests/seqdb/gi_mask).

`makeblastdb -gi_mask` writes the masks of each algorithm into GI-indexed files (`<name>.gmd`,
`.gnd`, `.gmo`, `.gno`, `.gmi`, `.gni`) instead of the mask column, and an alias file that names
them (`MASKLIST`). This script runs NCBI BLAST+ 2.17.0's makeblastdb on the cases of the M4c survey
(line NP, M/m4c/SURVEY.md section 4) and keeps what it wrote. NCBI is the comparison oracle only:
the Rust tests (`LOSAT/tests/makeblastdb_gi_mask.rs`, `LOSAT/src/writedb/gi_mask_tests.rs`) read the
committed files and run no NCBI program.

usage: make_gi_mask_fixtures.py --ncbi-bin DIR --work DIR [--out LOSAT/tests/seqdb/gi_mask]

Every case runs makeblastdb in a directory of its own (a case whose setup starts with `^` runs in
the directory of the case before it), with the input files copied in first, so that the names on
the command line, in the log and in the files are relative as a user types them.

One expectation leaves NCBI's bytes, in one place: `single_volume_exception` (PD-LOSAT-NCBI-DEFECTS
approved exception 23, the Owner's decision of 2026-11-20). NCBI's makeblastdb gives a
single-volume database built with `-gi_mask` an alias file whose DBLIST is the database's own
absolute path; no program can open that database, and makeblastdb itself ends with "BLAST Database
error: Illegal configuration: DB alias files are mutually recursive." (exit status 2, no metadata
file). LOSAT writes the database's file name without a directory there, which NCBI's reader takes
as the volume, and goes on as a build of several volumes does. So after each such run of NCBI this
script rewrites that one line in the run directory (before the next case of a `^` chain runs) and,
where NCBI ended on that error, takes the run as the successful run it then is: the line goes from
the standard error (or the `-logfile` file), the exit status is 0, and the metadata file is the
one NCBI's `blastdbcmd -metadata` writes for the database (the same `CSeqDB::GetDBMetaData` and
JSON writer that makeblastdb's own metadata step runs). Everything else is NCBI's run as it was.

Output files:
  in/                         the inputs: FASTA files (deterministic, from fixed seeds) and the mask
                              files that dustmasker, segmasker and convert2blastmask 2.17.0 wrote
                              for them (`-parse_seqids`, ASN.1 binary or text); `acc_dup.asnt` is
                              written by hand; `ncbi_single.nal` is the alias file that makeblastdb
                              wrote for a single-volume database with `-gi_mask` (`{C}` for its
                              directory), which LOSAT never writes: the cases x10 put it beside
                              a database as an earlier run of NCBI's would have left it;
                              `src4.*` is a version 4 database of nucl.fa with the mask column of
                              dust_nucl.asnb that makeblastdb wrote, the input of the cases s01 and
                              s02 (`-input_type blastdb`; `@SRC4` in the table of this script)
  cases.tsv                   one row per case: name, losat, setup, then the words after
                              `makeblastdb`, one per field. `losat`: `same` (LOSAT does as NCBI),
                              `refuse` (LOSAT refuses the input: two algorithms with one name that
                              would write one data file at the same time, through two streams of
                              NCBI's; see gi_mask.rs), `parser` (an
                              argument rule: LOSAT's parser refuses the words, exit status 2, where
                              NCBI's prints its usage). `setup`: space separated, `-` for none:
                              `F` copies in/F, `F>G` copies in/F to G, `F^N>G` writes the first N
                              FASTA records of in/F to G, `F=>G` copies in/F to G with `{C}`
                              replaced by the case directory (a file of a database, not an input:
                              it is listed and compared like the files of the run), `D/` makes the
                              directory D, `^` (first) keeps the directory of the case before.
                              `{C}` in a word is the absolute path of the case directory.
  expected/<case>.rc          makeblastdb's exit status (after `single_volume_exception`, as are
                              the other expectations)
  expected/<case>.stdout      its standard output, `.stderr` its standard error, with the case
                              directory as <C>, the time of "Building a new DB" as <T> and the
                              seconds of "Adding sequences" as <S>
  expected/<case>.files       the files of the case directory after the run (not the inputs), one
                              per line: path<TAB>size, `*` for a file whose size follows the date
                              or the directory of the run (.nin, .pin, .njs, .pjs, .nal, .pal, the
                              log file)
  db/<case>/<path>            the files kept: those that -gi_mask writes or changes (the GI-mask
                              files, the alias file, the metadata file, the mask column files and
                              the log file), and for the cases of FULL every file but the LMDB
                              files (.ndb, .pdb, .ntf, .ptf) and the index files (.nin, .pin), which
                              hold the date. In the alias, metadata and log files the case
                              directory is <C>. The other files are listed with their sizes only:
                              their bytes are those of a build without -gi_mask, which the fixtures
                              of the mask writer and of the database writer compare.

The "bytes-total" of every metadata file that a run writes is checked to be the sum of the sizes of
the files it lists (the Rust test compares LOSAT's with that sum).

Rerunning changes the dates in the `.gmi`/`.gni` files and the alias files (the tests mask them), and
in the files of `in/src4`.
"""
import argparse
import json
import os
import random
import re
import shutil
import subprocess
import time
from pathlib import Path

DATE_SIZED = (".nin", ".pin", ".njs", ".pjs", ".nal", ".pal")
NOT_KEPT = (".ndb", ".pdb", ".ntf", ".ptf", ".nin", ".pin")
GI_EXTENSIONS = (".gmd", ".gnd", ".gmo", ".gno", ".gmi", ".gni")
TEXT = (".nal", ".pal", ".njs", ".pjs")
KEPT = GI_EXTENSIONS + TEXT + (".naa", ".nab", ".nac", ".paa", ".pab", ".pac")
# The cases whose volume files are kept as well: single-volume version 5 and 4, several volumes
# of version 5 and 4, protein with one and with several volumes.
FULL = ("w02", "w05", "m01", "m10", "w17", "m09")


def wrap(s, w=60):
    return "\n".join(s[i:i + w] for i in range(0, len(s), w))


def write_inputs(d):
    """The FASTA inputs of the survey (M/m4c/scripts/mkinputs.py, seed 4; b2.sh, seed 7)."""
    rng = random.Random(4)

    def rnd(n, alpha="ACGT"):
        return "".join(rng.choice(alpha) for _ in range(n))

    nucl = [
        ("gi|1001|gb|AB000001.1| one", rnd(200) + "A" * 60 + rnd(200)),
        ("gi|1002|gb|AB000002.1| two", rnd(150) + "AT" * 30 + rnd(100) + "CAG" * 20 + rnd(150)),
        ("gi|2000000000|gb|AB000003.1| three", rnd(180) + "A" * 50 + rnd(170)),
        ("gi|3000000000|gb|AB000004.1| four above INT4_MAX", rnd(160) + "T" * 55 + rnd(190)),
        ("gb|AB000005.1| five no gi", rnd(120) + "A" * 60 + rnd(220)),
        ("gi|1006|emb|CAA00006.1| six", rnd(210) + "GT" * 30 + rnd(140)),
        ("gi|1007|gb|AB000007.1| seven no low complexity", rnd(400)),
        ("gi|999|gb|AB000008.1| eight\x01gi|1500|gb|AB000009.1| eight-b",
         rnd(130) + "A" * 70 + rnd(200)),
        ("lcl|seq9 nine local id", rnd(140) + "A" * 60 + rnd(200)),
        ("gi|5|gb|AB000010.1| ten small gi", rnd(100) + "C" * 60 + rnd(240)),
        ("gi|1011|ref|NM_000011.1| eleven", rnd(190) + "TTA" * 25 + rnd(160)),
        ("gi|1012 twelve gi only", rnd(170) + "A" * 60 + rnd(170)),
    ]
    with open(d / "nucl.fa", "w") as f:
        for h, s in nucl:
            f.write(">" + h + "\n" + wrap(s) + "\n")
    # nucl2.fa and the next records of mkinputs.py keep the generator's sequence.
    rnd(150), rnd(150), rnd(150), rnd(150)
    with open(d / "nogi.fa", "w") as f:
        f.write(">gb|AB003001.1| nogi one\n" + wrap(rnd(150) + "A" * 60 + rnd(150)) + "\n")
        f.write(">gb|AB003002.1| nogi two\n" + wrap(rnd(150) + "AT" * 30 + rnd(150)) + "\n")
    with open(d / "dupgi.fa", "w") as f:
        f.write(">gi|4001|gb|AB004001.1| dup gi first\n" + wrap(rnd(150) + "A" * 60 + rnd(150)) + "\n")
        f.write(">gi|4001|gb|AB004002.1| dup gi second\n" + wrap(rnd(100) + "AT" * 30 + rnd(200)) + "\n")
        f.write(">gi|4002|gb|AB004003.1| other\n" + wrap(rnd(100) + "C" * 60 + rnd(200)) + "\n")
    gis = list(range(10000, 10000 + 3 * 1100, 3))
    rng.shuffle(gis)
    with open(d / "big.fa", "w") as f:
        for k, gi in enumerate(gis):
            body = rnd(40) + "A" * 40 + rnd(40)
            if k % 7 == 0:
                body += "AT" * 25 + rnd(30)
            f.write(">gi|%d|gb|BG%06d.1| big %d\n%s\n" % (gi, k, k, wrap(body)))
    # The last three records of big.fa, and the records of nucl.fa without a GI: the inputs of
    # the masks of one algorithm that writes after another of its name (x12) or never (x11).
    big_records = (d / "big.fa").read_text().split(">")[1:]
    (d / "big_late.fa").write_text("".join(">" + r for r in big_records[-3:]))
    (d / "nucl_nogionly.fa").write_text("".join(
        ">" + h + "\n" + wrap(seq) + "\n" for h, seq in nucl if not h.startswith("gi|")
        and "|gi|" not in h))
    aa = "ACDEFGHIKLMNPQRSTVWY"
    prot = [
        ("gi|7001|gb|AAA00001.1| p one", rnd(80, aa) + "Q" * 25 + rnd(90, aa)),
        ("gi|7002|gb|AAA00002.1| p two", rnd(60, aa) + "PS" * 15 + rnd(70, aa) + "E" * 20 + rnd(40, aa)),
        ("gb|AAA00003.1| p three no gi", rnd(70, aa) + "Q" * 25 + rnd(80, aa)),
        ("gi|7004|gb|AAA00004.1| p four no lc", rnd(150, aa)),
        ("gi|7005|ref|NP_000005.1| p five", rnd(50, aa) + "K" * 22 + rnd(100, aa)),
    ]
    with open(d / "prot.fa", "w") as f:
        for h, s in prot:
            f.write(">" + h + "\n" + wrap(s) + "\n")

    r = random.Random(7)

    def rnd7(n):
        return "".join(r.choice("ACGT") for _ in range(n))

    a = rnd7(300) + "A" * 60 + rnd7(100)
    b = rnd7(300) + "A" * 60 + rnd7(100)
    c = rnd7(150) + "A" * 60 + rnd7(150)
    gi_seq, nogi_seq = ">gi|6001|gb|ZZ000001.1| gi seq\n", ">gb|ZZ000002.1| no gi seq\n"
    (d / "long.fa").write_text(gi_seq + wrap(a) + "\n" + nogi_seq + wrap(b) + "\n")
    (d / "short_nogi.fa").write_text(gi_seq + wrap(a) + "\n" + nogi_seq + wrap(b[:200]) + "\n")
    (d / "short_gi.fa").write_text(gi_seq + wrap(a[:200]) + "\n" + nogi_seq + wrap(b) + "\n")
    with open(d / "gimax.fa", "w") as f:
        f.write(">gi|4294967295|gb|ZZ000010.1| gi max\n" + wrap(c) + "\n")
        f.write(">gb|ZZ000011.1| no gi, no mask\n" + wrap(rnd7(360)) + "\n")
        f.write(">gi|10|gb|ZZ000012.1| gi 10\n" + wrap(rnd7(100) + "T" * 60 + rnd7(200)) + "\n")
        f.write(">gi|4294967296|gb|ZZ000013.1| gi above 32 bits\n"
                + wrap(rnd7(100) + "C" * 60 + rnd7(200)) + "\n")
    (d / "nolc.fa").write_text(">gi|6101|gb|ZZ000101.1| no low complexity\n" + wrap(rnd7(300)) + "\n")

    # nucl.fa with its first low-complexity stretch in lower case, for convert2blastmask.
    out = []
    for rec in (d / "nucl.fa").read_text().split(">")[1:]:
        h, _, body = rec.partition("\n")
        body = body.replace("\n", "")
        m = re.search("A{40,}|(AT){20,}|C{40,}|(GT){20,}", body)
        if m:
            body = body[:m.start()] + body[m.start():m.end()].lower() + body[m.end():]
        out.append(">" + h + "\n" + body + "\n")
    (d / "lower_nucl.fa").write_text("".join(out))
    # One file of the sequences with and without GIs, for version 4 (with version 5 and
    # -parse_seqids, sequences after those of nucl.fa end in an abort of makeblastdb that has
    # nothing to do with masks: "Input id list not in ascending oid order" of its LMDB writer).
    (d / "nucl_nogi.fa").write_text((d / "nucl.fa").read_text() + (d / "nogi.fa").read_text())
    # Three masks of the version 4 duplicate GI 4001 (b6.sh), named by accession.
    (d / "acc_dup.asnt").write_text(
        "Blast-db-mask-info ::= {\n  algo-id 2,\n  algo-program dust,\n  algo-options \"acc\",\n"
        "  masks {\n    masks {\n"
        "      int { from 10, to 19, id genbank { accession \"AB004001\", version 1 } },\n"
        "      int { from 30, to 39, id genbank { accession \"AB004002\", version 1 } },\n"
        "      int { from 50, to 59, id genbank { accession \"AB004003\", version 1 } }\n"
        "    },\n    more FALSE\n  }\n}\n")


def run(cmd, cwd, env=None):
    result = subprocess.run([str(c) for c in cmd], cwd=cwd, capture_output=True, env=env)
    if result.returncode != 0:
        raise SystemExit(f"{' '.join(map(str, cmd))}\n{result.stdout!r}{result.stderr!r}")


def write_masks(bin_dir, d):
    dust = [bin_dir / "dustmasker", "-infmt", "fasta", "-parse_seqids"]
    for name in ("nucl", "nogi", "dupgi", "big", "long", "gimax", "nolc"):
        run(dust + ["-in", f"{name}.fa", "-outfmt", "maskinfo_asn1_bin", "-out", f"dust_{name}.asnb"], d)
    run(dust + ["-in", "nucl.fa", "-outfmt", "maskinfo_asn1_bin", "-level", "10", "-out",
                "dust10_nucl.asnb"], d)
    run(dust + ["-in", "nucl.fa", "-outfmt", "maskinfo_asn1_text", "-out", "dust_nucl.asnt"], d)
    run([bin_dir / "dustmasker", "-infmt", "fasta", "-in", "nucl.fa", "-outfmt", "maskinfo_asn1_bin",
         "-out", "dust_nucl_noparse.asnb"], d)
    for name in ("big_late", "nucl_nogionly"):
        run(dust + ["-in", f"{name}.fa", "-outfmt", "maskinfo_asn1_bin", "-level", "10", "-out",
                    f"dust10_{name}.asnb"], d)
        (d / f"{name}.fa").unlink()
    run([bin_dir / "segmasker", "-infmt", "fasta", "-parse_seqids", "-in", "prot.fa", "-outfmt",
         "maskinfo_asn1_bin", "-out", "seg_prot.asnb"], d)
    run([bin_dir / "convert2blastmask", "-in", "lower_nucl.fa", "-masking_algorithm", "mymask",
         "-masking_options", "opt:a", "-outfmt", "maskinfo_asn1_bin", "-parse_seqids", "-out",
         "conv_other.asnb"], d)
    (d / "lower_nucl.fa").unlink()


N = ["-in", "nucl.fa", "-dbtype", "nucl", "-parse_seqids", "-out", "db", "-title", "t"]
BN = ["-in", "big.fa", "-dbtype", "nucl", "-parse_seqids", "-out", "db", "-title", "t"]
# -input_type blastdb from in/src4.* (`write_source_db`; `@SRC4` in a setup lists its files).
S = ["-input_type", "blastdb", "-in", "src4", "-dbtype", "nucl", "-parse_seqids", "-out", "db",
     "-title", "t", "-gi_mask"]
NM = N + ["-max_file_sz", "400"]
P = ["-in", "prot.fa", "-dbtype", "prot", "-parse_seqids", "-out", "db", "-title", "pt"]
DN = "nucl.fa dust_nucl.asnb"
DN2 = "nucl.fa dust_nucl.asnb dust10_nucl.asnb"
GM = ["-gi_mask"]


def big(n):
    src = "big.fa" if n == 1100 else f"big.fa^{n}>big{n}.fa"
    name = "big.fa" if n == 1100 else f"big{n}.fa"
    return (f"big{n}", "same", f"{src} dust_big.asnb",
            ["-in", name, "-dbtype", "nucl", "-parse_seqids", "-out", "db", "-title", "t",
             "-mask_data", "dust_big.asnb", "-gi_mask", "-gi_mask_name", "big"])


BIG = ["-in", "big.fa", "-dbtype", "nucl", "-parse_seqids", "-out", "db", "-title", "t",
       "-mask_data", "dust_big.asnb", "-gi_mask", "-gi_mask_name", "big"]
L60 = "n" * 60

# name, losat, setup, words (survey M4c section 4; the probe names of M/m4c/probe/b1..b6).
CASES = [
    ("w02", "same", DN, N + ["-mask_data", "dust_nucl.asnb"] + GM),
    ("w04", "same", DN, N + ["-mask_data", "dust_nucl.asnb", "-gi_mask", "-gi_mask_name", "gm1"]),
    ("w05", "same", DN, N + ["-blastdb_version", "4", "-mask_data", "dust_nucl.asnb"] + GM),
    ("w06", "same", DN, N + ["-blastdb_version", "4", "-mask_data", "dust_nucl.asnb", "-gi_mask",
                             "-gi_mask_name", "gm1"]),
    ("w08", "same", DN2, N + ["-mask_data", "dust_nucl.asnb,dust10_nucl.asnb", "-gi_mask",
                              "-gi_mask_name", "m1,m2"]),
    ("w09", "same", DN2, N + ["-mask_data", "dust_nucl.asnb,dust10_nucl.asnb"] + GM),
    ("w10", "same", "nucl.fa dust_nucl.asnt",
     N + ["-mask_data", "dust_nucl.asnt", "-gi_mask", "-gi_mask_name", "mt"]),
    ("w11", "same", "nucl.fa masks/ dust_nucl.asnb>masks/dust_nucl.asnb",
     N + ["-mask_data", "masks/dust_nucl.asnb"] + GM),
    ("w12", "same", DN, ["-in", "nucl.fa", "-dbtype", "nucl", "-parse_seqids", "-out", "sub/db",
                         "-title", "t", "-mask_data", "dust_nucl.asnb", "-gi_mask",
                         "-gi_mask_name", "gm1"]),
    ("w13", "same", DN + " sub2/",
     N + ["-mask_data", "dust_nucl.asnb", "-gi_mask", "-gi_mask_name", "sub2/gm1"]),
    ("w14", "same", DN, N + ["-mask_data", "dust_nucl.asnb", "-gi_mask", "-gi_mask_name",
                             "nodir/gm1"]),
    ("w17", "same", "prot.fa seg_prot.asnb",
     P + ["-mask_data", "seg_prot.asnb", "-gi_mask", "-gi_mask_name", "pm1"]),
    ("w18", "same", "prot.fa seg_prot.asnb",
     P + ["-blastdb_version", "4", "-mask_data", "seg_prot.asnb"] + GM),
    ("w21", "same", "nogi.fa dust_nogi.asnb",
     ["-in", "nogi.fa", "-dbtype", "nucl", "-parse_seqids", "-out", "db", "-title", "t",
      "-mask_data", "dust_nogi.asnb", "-gi_mask", "-gi_mask_name", "ng"]),
    ("w22", "same", DN, N + ["-mask_data", "missing.asnb"] + GM),
    ("w23", "same", DN, N + ["-mask_data", "dust_nucl.asnb", "-gi_mask", "-gi_mask_name", "a b"]),
    ("w24", "same", DN + " dust_nucl.asnb>dust_copy.asnb",
     N + ["-mask_data", "dust_nucl.asnb,dust_copy.asnb", "-gi_mask", "-gi_mask_name", "x,x"]),
    ("w25", "refuse", DN2, N + ["-mask_data", "dust_nucl.asnb,dust10_nucl.asnb", "-gi_mask",
                                "-gi_mask_name", "x,x"]),
    ("w25b", "refuse", DN2, N + ["-mask_data", "dust_nucl.asnb,dust10_nucl.asnb", "-gi_mask",
                                 "-gi_mask_name", "x,./x"]),
    ("w26", "same", DN, N + ["-mask_data", "dust_nucl.asnb,dust_nucl.asnb"] + GM),
    ("w27", "same", DN, N + ["-mask_data", "dust_nucl.asnb", "-gi_mask", "-gi_mask_name", ""]),
    ("w28", "same", DN, N + ["-mask_data", "dust_nucl.asnb", "-gi_mask", "-gi_mask_name", "a,"]),
    ("w29", "same", DN, N + ["-mask_data", "dust_nucl.asnb", "-gi_mask", "-gi_mask_name",
                             " gm1 "]),
    ("w30", "same", "dupgi.fa dust_dupgi.asnb",
     ["-in", "dupgi.fa", "-dbtype", "nucl", "-parse_seqids", "-out", "db", "-title", "t",
      "-mask_data", "dust_dupgi.asnb", "-gi_mask", "-gi_mask_name", "dg"]),
    ("w31", "same", DN, N + GM),
    ("w32", "same", DN, N + ["-mask_data", "dust_nucl.asnb", "-gi_mask", "-gi_mask_name", "db"]),
    ("w33", "same", DN, N + ["-mask_data", "dust_nucl.asnb", "-gi_mask", "-gi_mask_name", "gm1",
                             "-logfile", "log.txt"]),
    ("e01", "parser", DN, ["-in", "nucl.fa", "-dbtype", "nucl", "-out", "db", "-mask_data",
                           "dust_nucl.asnb", "-gi_mask"]),
    ("e02", "parser", DN, N + ["-mask_data", "dust_nucl.asnb", "-gi_mask", "-mask_id", "foo"]),
    ("e03", "parser", DN, N + ["-mask_data", "dust_nucl.asnb", "-gi_mask", "-mask_desc", "d"]),
    ("e04", "parser", DN, N + ["-mask_data", "dust_nucl.asnb", "-gi_mask_name", "x"]),
    ("e05", "parser", DN, N + ["-gi_mask", "-gi_mask_name", "x"]),
    ("e06", "same", DN, N + ["-mask_data", "dust_nucl.asnb", "-gi_mask", "-gi_mask_name", "a,b"]),
    ("e07", "same", DN2, N + ["-mask_data", "dust_nucl.asnb,dust10_nucl.asnb", "-gi_mask",
                              "-gi_mask_name", "a"]),
    ("e08", "parser", DN, N + ["-mask_data", "dust_nucl.asnb", "-gi_mask", "true"]),
    ("e09", "parser", DN, ["-in", "nucl.fa", "-dbtype", "nucl", "-out", "db", "-mask_data",
                           "dust_nucl.asnb", "-gi_mask_name", "x"]),
    big(1), big(512), big(513), big(1024), big(1025), big(1100),
    ("bigmax", "same", "big.fa dust_big.asnb", BIG + ["-max_file_sz", "4000"]),
    ("bigmax4", "same", "big.fa dust_big.asnb", BIG + ["-blastdb_version", "4", "-max_file_sz",
                                                        "4000"]),
    ("x01_1", "same", DN, N + ["-mask_data", "dust_nucl.asnb", "-gi_mask", "-gi_mask_name", "gm1"]),
    ("x01_2", "same", "^", N + ["-mask_data", "dust_nucl.asnb", "-gi_mask", "-gi_mask_name", "gm1"]),
    ("x01_3", "same", "^", N + ["-mask_data", "dust_nucl.asnb"]),
    ("x02", "same", "gimax.fa dust_gimax.asnb",
     ["-in", "gimax.fa", "-dbtype", "nucl", "-parse_seqids", "-out", "db", "-title", "t",
      "-mask_data", "dust_gimax.asnb", "-gi_mask", "-gi_mask_name", "gx"]),
    ("x03nogi", "same", "short_nogi.fa dust_long.asnb",
     ["-in", "short_nogi.fa", "-dbtype", "nucl", "-parse_seqids", "-out", "db", "-title", "t",
      "-mask_data", "dust_long.asnb", "-gi_mask", "-gi_mask_name", "m"]),
    ("x03gi", "same", "short_gi.fa dust_long.asnb",
     ["-in", "short_gi.fa", "-dbtype", "nucl", "-parse_seqids", "-out", "db", "-title", "t",
      "-mask_data", "dust_long.asnb", "-gi_mask", "-gi_mask_name", "m"]),
    ("x04a", "same", "nucl_nogi.fa dust_nucl.asnb dust_nogi.asnb",
     ["-in", "nucl_nogi.fa", "-dbtype", "nucl", "-parse_seqids", "-out", "db", "-title", "t",
      "-blastdb_version", "4", "-mask_data", "dust_nucl.asnb,dust_nogi.asnb", "-gi_mask", "-gi_mask_name", "m1,m2"]),
    ("x04b", "same", "nucl_nogi.fa dust_nucl.asnb dust_nogi.asnb",
     ["-in", "nucl_nogi.fa", "-dbtype", "nucl", "-parse_seqids", "-out", "db", "-title", "t",
      "-blastdb_version", "4", "-mask_data", "dust_nogi.asnb,dust_nucl.asnb", "-gi_mask", "-gi_mask_name", "m1,m2"]),
    ("x04c", "same", "nogi.fa dust_nogi.asnb",
     ["-in", "nogi.fa", "-dbtype", "nucl", "-parse_seqids", "-out", "db", "-title", "t",
      "-mask_data", "dust_nogi.asnb", "-gi_mask", "-gi_mask_name", "m1"]),
    ("x05", "same", "nucl.fa conv_other.asnb",
     N + ["-mask_data", "conv_other.asnb", "-gi_mask", "-gi_mask_name", "oth"]),
    ("x05n", "same", "nucl.fa conv_other.asnb", N + ["-mask_data", "conv_other.asnb"] + GM),
    ("x07", "same", DN, N + ["-mask_data", "dust_nucl.asnb", "-gi_mask", "-gi_mask_name",
                             "{C}/abs_gm"]),
    ("x08", "same", "nolc.fa dust_nolc.asnb",
     ["-in", "nolc.fa", "-dbtype", "nucl", "-parse_seqids", "-out", "db", "-title", "t",
      "-mask_data", "dust_nolc.asnb", "-gi_mask", "-gi_mask_name", "e"]),
    ("x09", "same", "nucl.fa dust_nucl_noparse.asnb",
     N + ["-mask_data", "dust_nucl_noparse.asnb", "-gi_mask", "-gi_mask_name", "np"]),
    # A database beside the alias file that an earlier run of NCBI's makeblastdb wrote for a
    # single-volume database with -gi_mask (survey D2): the removal of the old database cannot
    # read it; a build with -gi_mask writes its own alias file over it, a build without -gi_mask
    # leaves it and ends on it as NCBI's makeblastdb does.
    ("x10_1", "same", DN, N + ["-mask_data", "dust_nucl.asnb"]),
    ("x10_2", "same", "^ ncbi_single.nal=>db.nal",
     N + ["-mask_data", "dust_nucl.asnb", "-gi_mask", "-gi_mask_name", "gm1"]),
    ("x10_3", "same", "^ ncbi_single.nal=>db.nal", N + ["-mask_data", "dust_nucl.asnb"]),
    # Two algorithms of one name (x11, x12): one never writes, or each data volume is written by
    # one algorithm after the other's closed it (the files are NCBI's; w25, w25b and m04 write
    # one data file at the same time and are refused).
    ("x11_1", "same", "nucl.fa dust_nucl.asnb dust10_nucl_nogionly.asnb",
     N + ["-mask_data", "dust_nucl.asnb,dust10_nucl_nogionly.asnb", "-gi_mask", "-gi_mask_name",
          "x,x"]),
    ("x11_2", "same", "nucl.fa dust_nucl.asnb dust10_nucl_nogionly.asnb",
     N + ["-mask_data", "dust10_nucl_nogionly.asnb,dust_nucl.asnb", "-gi_mask", "-gi_mask_name",
          "x,x"]),
    ("x12", "same", "big.fa dust_big.asnb dust10_big_late.asnb",
     BN + ["-max_file_sz", "4000", "-mask_data", "dust_big.asnb,dust10_big_late.asnb",
           "-gi_mask", "-gi_mask_name", "x,x"]),
    # -input_type blastdb: the input's algorithms write GI masks under the GIs of every
    # definition line of a sequence; s02 with an algorithm of -mask_data first.
    ("s01", "same", "@SRC4", S),
    ("s02", "same", "@SRC4 dust10_nucl.asnb",
     S + ["-mask_data", "dust10_nucl.asnb", "-gi_mask_name", "m"]),
    # An exclusion is reported before the missing -dbtype (x_PostCheck's first walk).
    ("e10", "parser", DN, ["-in", "nucl.fa", "-taxid", "1", "-taxid_map", "nucl.fa",
                           "-parse_seqids", "-out", "db"]),
    ("e11", "parser", DN, ["-in", "nucl.fa", "-gi_mask", "-mask_id", "foo", "-mask_data",
                           "dust_nucl.asnb", "-parse_seqids", "-out", "db"]),
    ("m01", "same", DN2, NM + ["-mask_data", "dust_nucl.asnb,dust10_nucl.asnb", "-gi_mask",
                               "-gi_mask_name", "m1,m2"]),
    ("m03", "same", DN, NM + GM),
    ("m04", "refuse", DN2, NM + ["-mask_data", "dust_nucl.asnb,dust10_nucl.asnb", "-gi_mask",
                                 "-gi_mask_name", "x,x"]),
    ("m05", "same", "nucl.fa dust_nogi.asnb dust_nucl.asnb",
     NM + ["-mask_data", "dust_nogi.asnb,dust_nucl.asnb", "-gi_mask", "-gi_mask_name", "m1,m2"]),
    ("m06", "same", "nucl.fa conv_other.asnb",
     NM + ["-mask_data", "conv_other.asnb", "-gi_mask", "-gi_mask_name", "oth"]),
    ("m06n", "same", "nucl.fa conv_other.asnb", NM + ["-mask_data", "conv_other.asnb"] + GM),
    ("m07", "same", "gimax.fa dust_gimax.asnb",
     ["-in", "gimax.fa", "-dbtype", "nucl", "-parse_seqids", "-out", "db", "-title", "t",
      "-max_file_sz", "100", "-mask_data", "dust_gimax.asnb", "-gi_mask", "-gi_mask_name", "gx"]),
    ("m08", "same", DN, ["-in", "nucl.fa", "-dbtype", "nucl", "-parse_seqids", "-out", "sub/db",
                         "-title", "t", "-max_file_sz", "400", "-mask_data", "dust_nucl.asnb",
                         "-gi_mask", "-gi_mask_name", "gm1"]),
    ("m09", "same", "prot.fa seg_prot.asnb",
     P + ["-max_file_sz", "200", "-mask_data", "seg_prot.asnb", "-gi_mask", "-gi_mask_name", "pm1"]),
    ("m10", "same", DN2, NM + ["-blastdb_version", "4", "-mask_data",
                               "dust_nucl.asnb,dust10_nucl.asnb", "-gi_mask", "-gi_mask_name",
                               "m1,m2"]),
    ("c05", "same", "nucl.fa nogi.fa dust_nucl.asnb dust_nogi.asnb",
     ["-in", "nucl.fa nogi.fa", "-dbtype", "nucl", "-parse_seqids", "-out", "db", "-title", "t",
      "-blastdb_version", "4", "-max_file_sz", "400", "-mask_data",
      "dust_nucl.asnb,dust_nogi.asnb", "-gi_mask", "-gi_mask_name", "m1,m2"]),
    ("d01", "same", "dupgi.fa dust_dupgi.asnb",
     ["-in", "dupgi.fa", "-dbtype", "nucl", "-parse_seqids", "-out", "db", "-title", "t",
      "-blastdb_version", "4", "-max_file_sz", "200", "-mask_data", "dust_dupgi.asnb",
      "-gi_mask", "-gi_mask_name", "dg"]),
    ("d02_10", "same", "nucl.fa dust10_nucl.asnb",
     N + ["-max_file_sz", "10", "-mask_data", "dust10_nucl.asnb", "-gi_mask", "-gi_mask_name", "m"]),
    ("d02_20", "same", "nucl.fa dust10_nucl.asnb",
     N + ["-max_file_sz", "20", "-mask_data", "dust10_nucl.asnb", "-gi_mask", "-gi_mask_name", "m"]),
    ("d02_30", "same", "nucl.fa dust10_nucl.asnb",
     N + ["-max_file_sz", "30", "-mask_data", "dust10_nucl.asnb", "-gi_mask", "-gi_mask_name", "m"]),
    ("d03", "same", "nucl.fa conv_other.asnb",
     NM + ["-mask_data", "conv_other.asnb", "-gi_mask", "-gi_mask_name", "abcdefgh"]),
    ("d04", "same", "nucl.fa conv_other.asnb",
     NM + ["-mask_data", "conv_other.asnb", "-gi_mask", "-gi_mask_name", L60]),
    ("g01", "same", "dupgi.fa acc_dup.asnt",
     ["-in", "dupgi.fa", "-dbtype", "nucl", "-parse_seqids", "-out", "db", "-title", "t",
      "-blastdb_version", "4", "-max_file_sz", "200", "-mask_data", "acc_dup.asnt", "-gi_mask",
      "-gi_mask_name", "ad"]),
    ("g02", "same", DN, NM + ["-mask_data", "dust_nucl.asnb", "-gi_mask", "-gi_mask_name",
                              "{C}/abs_gm"]),
]


def setup(case_dir, tokens, inputs):
    """Makes the case directory ready; returns the paths (relative) of the inputs it put there."""
    put = []
    for token in tokens:
        if token.endswith("/"):
            (case_dir / token).mkdir(parents=True, exist_ok=True)
            continue
        if "=>" in token:
            src, _, dst = token.partition("=>")
            text = (inputs / src).read_bytes().replace(b"{C}", str(case_dir).encode())
            (case_dir / dst).write_bytes(text)
            continue
        src, _, dst = token.partition(">")
        dst = dst or src
        src, _, count = src.partition("^")
        if count:
            records = (inputs / src).read_text().split(">")[1:int(count) + 1]
            (case_dir / dst).write_text("".join(">" + r for r in records))
        else:
            shutil.copyfile(inputs / src, case_dir / dst)
        put.append(dst)
    return put


RECURSION = b"BLAST Database error: Illegal configuration: DB alias files are mutually recursive.\n"


def write_source_db(bin_dir, d, env):
    """`in/src4.*`: a version 4 database of nucl.fa with the mask column of dust_nucl.asnb, the
    input of the cases of `-input_type blastdb`. Returns the names of its files (`@SRC4`)."""
    words = N[:-4] + ["-blastdb_version", "4", "-mask_data", "dust_nucl.asnb", "-out", "src4",
                      "-title", "src"]
    result = subprocess.run([str(bin_dir / "makeblastdb")] + words, cwd=d, capture_output=True,
                            env=env)
    if result.returncode != 0 or result.stderr:
        raise SystemExit(f"src4: {result.returncode} {result.stderr!r}")
    return " ".join(sorted(p.name for p in d.glob("src4.*")))


def metadata_sizes(case_dir):
    """The sizes of the metadata files under `case_dir`, taken before a run."""
    return {path: path.stat().st_size for path in case_dir.rglob("*.[np]js")}


def check_bytes_total(case_dir, stamp, before):
    """The "bytes-total" of every metadata file that the run wrote is the sum of the sizes of the
    files it lists, each time it lists one, the metadata file itself at its size before the run
    (`before`: makeblastdb lists the files before it writes the new one, so an old one that the
    run did not remove counts as it was). The Rust test compares LOSAT's distance from that sum
    with 0."""
    for path in sorted(case_dir.rglob("*.[np]js")):
        if path.stat().st_mtime_ns < stamp.stat().st_mtime_ns:
            continue
        meta = json.loads(path.read_text())
        total = 0
        for name in meta["files"]:
            listed = path.parent / name
            if listed == path:
                if path not in before:
                    raise SystemExit(f"{path}: lists itself, which did not exist before the run")
                total += before[path]
            else:
                total += listed.stat().st_size
        if meta["bytes-total"] != total:
            raise SystemExit(f"{path}: bytes-total {meta['bytes-total']}, files {total}")


def write_ncbi_single_alias(bin_dir, inputs, work, env):
    """`in/ncbi_single.nal`: the alias file that NCBI's makeblastdb writes for a single-volume
    database built with `-gi_mask` (DBLIST = the database's absolute path, `{C}/db` here), as an
    earlier run of NCBI's leaves it in a directory."""
    run_dir = work / "ncbi_single"
    run_dir.mkdir()
    for name in ("nucl.fa", "dust_nucl.asnb"):
        shutil.copyfile(inputs / name, run_dir / name)
    words = N + ["-mask_data", "dust_nucl.asnb", "-gi_mask", "-gi_mask_name", "gm1"]
    result = subprocess.run([str(bin_dir / "makeblastdb")] + words, cwd=run_dir,
                            capture_output=True, env=env)
    if result.returncode != 2 or result.stderr != RECURSION:
        raise SystemExit(f"ncbi_single: {result.returncode} {result.stderr!r}")
    alias = (run_dir / "db.nal").read_bytes()
    (inputs / "ncbi_single.nal").write_bytes(alias.replace(str(run_dir).encode(), b"{C}"))


def single_volume_exception(bin_dir, case_dir, argv, rc, stderr, stamp, env):
    """The one place where an expectation leaves NCBI's bytes (PD-LOSAT-NCBI-DEFECTS approved
    exception 23; the module's docstring). Where this run of makeblastdb with `-gi_mask` wrote the
    alias file of a single-volume database (DBLIST = the database's absolute path), that line
    becomes the database's file name without a directory; where the run then ended on "DB alias
    files are mutually recursive" (exit status 2), the line goes from the standard error or the
    `-logfile` file, the exit status is 0, and the metadata file is NCBI's `blastdbcmd -metadata`
    of the database. `stamp` is a file touched before the run. Returns the exit status and the
    standard error."""
    if "-gi_mask" not in argv or "-dbtype" not in argv:
        # Without -dbtype the argument reader stops before anything is written (e10, e11).
        return rc, stderr
    out = argv[argv.index("-out") + 1]
    dbtype = argv[argv.index("-dbtype") + 1]
    alias = case_dir / (out + (".pal" if dbtype == "prot" else ".nal"))
    # An alias file that this run did not write (one that was there before) stays as it is.
    if not alias.exists() or alias.stat().st_mtime_ns <= stamp.stat().st_mtime_ns:
        return rc, stderr
    lines = alias.read_bytes().split(b"\n")
    at = next(i for i, line in enumerate(lines) if line.startswith(b"DBLIST "))
    path = lines[at][len(b"DBLIST "):]
    if not os.path.isabs(path):
        # Several volumes: NCBI's list of their file names, which NCBI's reader opens.
        return rc, stderr
    lines[at] = b"DBLIST " + os.path.basename(path)
    alias.write_bytes(b"\n".join(lines))
    if rc != 2:
        # The run failed before its metadata step (a duplicate GI, a range past a sequence's end):
        # only the alias file differs.
        return rc, stderr
    if "-logfile" in argv:
        log = case_dir / argv[argv.index("-logfile") + 1]
        text = log.read_bytes()
        if RECURSION not in text:
            raise SystemExit(f"{case_dir}: exit status 2 without the alias error in {log}")
        log.write_bytes(text.replace(RECURSION, b"", 1))
    else:
        if RECURSION not in stderr:
            raise SystemExit(f"{case_dir}: exit status 2 without the alias error: {stderr!r}")
        stderr = stderr.replace(RECURSION, b"", 1)
    meta = subprocess.run([str(bin_dir / "blastdbcmd"), "-db", out, "-dbtype", dbtype, "-metadata"],
                          cwd=case_dir, capture_output=True, env=env)
    if meta.returncode != 0 or meta.stderr:
        raise SystemExit(f"{case_dir}: blastdbcmd -metadata: {meta.returncode} {meta.stderr!r}")
    (case_dir / (out + (".pjs" if dbtype == "prot" else ".njs"))).write_bytes(meta.stdout)
    return 0, stderr


def masked_text(text, case_dir):
    """A stream or a text file of the run with the case directory, the time and the seconds
    masked."""
    text = text.replace(str(case_dir), "<C>")
    text = re.sub(r"^(Building a new DB, current time: ).*$", r"\1<T>", text, flags=re.M)
    text = re.sub(r"^(Adding sequences from (?:FASTA|raw db source); added \d+ sequences in )\S+"
                  r"( seconds\.)$",
                  r"\1<S>\2", text, flags=re.M)
    return text


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ncbi-bin", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--out", default=str(Path(__file__).resolve().parent / "gi_mask"))
    args = ap.parse_args()
    bin_dir = Path(args.ncbi_bin).resolve()
    work = Path(args.work).resolve()
    out = Path(args.out)
    if work.exists():
        raise SystemExit(f"{work} exists; give a new work directory")
    inputs = work / "in"
    inputs.mkdir(parents=True)
    env = dict(os.environ, BLAST_USAGE_REPORT="false")
    write_inputs(inputs)
    write_masks(bin_dir, inputs)
    write_ncbi_single_alias(bin_dir, inputs, work, env)
    src4 = write_source_db(bin_dir, inputs, env)
    if out.exists():
        shutil.rmtree(out)
    for sub in ("in", "expected", "db"):
        (out / sub).mkdir(parents=True)
    for f in sorted(inputs.iterdir()):
        shutil.copyfile(f, out / "in" / f.name)

    rows = []
    case_dir = None
    put = []
    stamp = work / "stamp"
    for name, losat, setup_text, words in CASES:
        setup_text = setup_text.replace("@SRC4", src4)
        tokens = setup_text.split()
        if tokens and tokens[0] == "^":
            tokens = tokens[1:]
        else:
            case_dir = work / "cases" / name
            case_dir.mkdir(parents=True)
            put = []
        put += setup(case_dir, tokens, inputs)
        argv = [w.replace("{C}", str(case_dir)) for w in words]
        before = metadata_sizes(case_dir)
        stamp.touch()
        time.sleep(0.05)
        result = subprocess.run([str(bin_dir / "makeblastdb")] + argv, cwd=case_dir,
                                capture_output=True, env=env)
        rc, stderr = single_volume_exception(bin_dir, case_dir, argv, result.returncode,
                                             result.stderr, stamp, env)
        check_bytes_total(case_dir, stamp, before)
        rows.append("\t".join([name, losat, setup_text] + words))
        exp = out / "expected"
        (exp / f"{name}.rc").write_text(f"{rc}\n")
        (exp / f"{name}.stdout").write_text(masked_text(result.stdout.decode("latin-1"), case_dir),
                                           encoding="latin-1")
        (exp / f"{name}.stderr").write_text(masked_text(stderr.decode("latin-1"), case_dir),
                                           encoding="latin-1")
        files = []
        keep_dir = out / "db" / name
        for path in sorted(case_dir.rglob("*")):
            rel = path.relative_to(case_dir).as_posix()
            if path.is_dir() or rel in put:
                continue
            sized = rel.endswith(DATE_SIZED) or rel == "log.txt"
            files.append(f"{rel}\t{'*' if sized else path.stat().st_size}")
            if rel.endswith(NOT_KEPT) or not (name in FULL or rel.endswith(KEPT)
                                              or rel == "log.txt"):
                continue
            target = keep_dir / rel
            target.parent.mkdir(parents=True, exist_ok=True)
            if rel.endswith(TEXT) or rel == "log.txt":
                target.write_bytes(masked_text(path.read_bytes().decode("latin-1"), case_dir)
                                   .encode("latin-1"))
            else:
                shutil.copyfile(path, target)
        (exp / f"{name}.files").write_text("".join(f + "\n" for f in files))
    (out / "cases.tsv").write_text("".join(r + "\n" for r in rows))
    print(f"{len(rows)} cases")


if __name__ == "__main__":
    main()
