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

Output files:
  in/                         the inputs: FASTA files (deterministic, from fixed seeds) and the mask
                              files that dustmasker, segmasker and convert2blastmask 2.17.0 wrote
                              for them (`-parse_seqids`, ASN.1 binary or text); `acc_dup.asnt` is
                              written by hand
  cases.tsv                   one row per case: name, losat, setup, then the words after
                              `makeblastdb`, one per field. `losat`: `same` (LOSAT does as NCBI),
                              `refuse` (LOSAT refuses the input: two algorithms with one name, whose
                              files NCBI writes through two streams; see gi_mask.rs), `parser` (an
                              argument rule: LOSAT's parser refuses the words, exit status 2, where
                              NCBI's prints its usage). `setup`: space separated, `-` for none:
                              `F` copies in/F, `F>G` copies in/F to G, `F^N>G` writes the first N
                              FASTA records of in/F to G, `D/` makes the directory D, `^` (first)
                              keeps the directory of the case before. `{C}` in a word is the
                              absolute path of the case directory.
  expected/<case>.rc          makeblastdb's exit status
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

Rerunning changes the dates in the `.gmi`/`.gni` files and the alias files (the tests mask them).
"""
import argparse
import os
import random
import re
import shutil
import subprocess
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
    run([bin_dir / "segmasker", "-infmt", "fasta", "-parse_seqids", "-in", "prot.fa", "-outfmt",
         "maskinfo_asn1_bin", "-out", "seg_prot.asnb"], d)
    run([bin_dir / "convert2blastmask", "-in", "lower_nucl.fa", "-masking_algorithm", "mymask",
         "-masking_options", "opt:a", "-outfmt", "maskinfo_asn1_bin", "-parse_seqids", "-out",
         "conv_other.asnb"], d)
    (d / "lower_nucl.fa").unlink()


N = ["-in", "nucl.fa", "-dbtype", "nucl", "-parse_seqids", "-out", "db", "-title", "t"]
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


def masked_text(text, case_dir):
    """A stream or a text file of the run with the case directory, the time and the seconds
    masked."""
    text = text.replace(str(case_dir), "<C>")
    text = re.sub(r"^(Building a new DB, current time: ).*$", r"\1<T>", text, flags=re.M)
    text = re.sub(r"^(Adding sequences from FASTA; added \d+ sequences in )\S+( seconds\.)$",
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
    write_inputs(inputs)
    write_masks(bin_dir, inputs)
    if out.exists():
        shutil.rmtree(out)
    for sub in ("in", "expected", "db"):
        (out / sub).mkdir(parents=True)
    for f in sorted(inputs.iterdir()):
        shutil.copyfile(f, out / "in" / f.name)

    env = dict(os.environ, BLAST_USAGE_REPORT="false")
    rows = []
    case_dir = None
    put = []
    for name, losat, setup_text, words in CASES:
        tokens = setup_text.split()
        if tokens and tokens[0] == "^":
            tokens = tokens[1:]
        else:
            case_dir = work / "cases" / name
            case_dir.mkdir(parents=True)
            put = []
        put += setup(case_dir, tokens, inputs)
        argv = [w.replace("{C}", str(case_dir)) for w in words]
        result = subprocess.run([str(bin_dir / "makeblastdb")] + argv, cwd=case_dir,
                                capture_output=True, env=env)
        rows.append("\t".join([name, losat, setup_text] + words))
        exp = out / "expected"
        (exp / f"{name}.rc").write_text(f"{result.returncode}\n")
        (exp / f"{name}.stdout").write_text(masked_text(result.stdout.decode("latin-1"), case_dir),
                                           encoding="latin-1")
        (exp / f"{name}.stderr").write_text(masked_text(result.stderr.decode("latin-1"), case_dir),
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
