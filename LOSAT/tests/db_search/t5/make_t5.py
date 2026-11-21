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
record with the suffix goes to the version 4 database `psuffix` of its own (written before LOSAT's
searches wrote such a heading; rows `db.t5.psuffix.*` search it since line RG item T5c).

With `--suffix-only` (`python3 tests/db_search/t5/make_t5.py NCBI_BIN --suffix-only`, line RG item
T5c) it writes only the protein databases of titles ending in `]` whose end NCBI's
`x_AdjustProteinTitleSuffix` rewrites from the Bioseq's first source (SUFFIX, SUFFIX_STORED): `psfx`
(version 5, `-parse_seqids`, `-taxid_map`: a titled line and an untitled one, whose outfmt 0 heading
shows the generated title, and `lcl|Subject_9nn` records, whose tabular ids are its first word) and
the version 4 databases `psfs_<taxid>` without `-parse_seqids` (`-taxid`, `psfs_none` without one)
with the alias `psfs` over them (the tabular ids are the first word of the generated title).

With `--word-only` (`python3 tests/db_search/t5/make_t5.py NCBI_BIN --word-only`, fix 1 of line RG
item T5c) it writes only the version 4 protein databases `pword_<taxid>` without `-parse_seqids`
(WORD_STORED; `pword_none` without `-taxid`) and the alias `pword` over them: titles whose generated
title has no first word, so that NCBI keeps the sequence's own `gnl|BL_ORD_ID|n` in the tabular ids
(`BL_ORD_ID:n` for `sacc` and `saccver`): a title of one space, a prefix alone (`TPA:`), and
`MAG: (plastid) [ ]`, which `x_AdjustProteinTitleSuffix` cuts to `MAG:` where the Bioseq has a
source (taxid 77777, which no taxonomy database here names); `  MAG:  ` and `chi, partial [ ]` give a
word.

With `--taxid-only` (`python3 tests/db_search/t5/make_t5.py NCBI_BIN --taxid-only`, fix 2 of line
RG item T5c) it writes only the version 4 databases `ptx_<name>` and `ntx_<name>` of TAXID_ENCODINGS
without `-parse_seqids` (`-taxid 9606`; two records, the letters of `pdup.faa` / `ndup.fna`): the
`taxid [2]` of the first record (OID 0) holds another encoding than makeblastdb's `02 02 25 86`
(the record and the header offsets of the index rewritten): the application INTEGER tag `0x42`
and nine content bytes of sign extension, which NCBI's `ReadInt4` reads as 9606, and 2^31 in five
bytes, which it fails on (`eOverflow`, `TTaxId` is an `int`) wherever it reads the header (probe
`remaining-gaps-20261112/T/taxid_t5c.md` of line RG).

With `--sfu-only` (`python3 tests/db_search/t5/make_t5.py NCBI_BIN --sfu-only`, merge 4 of line RG
item T5c with item T7) it writes only the version 4 protein databases of titles ending in ` [ ]` with
the taxid 100009, whose name (`Bad \x81 name`) NCBI's `x_AdjustProteinTitleSuffix` writes at their end
and `GenerateDefline` then cannot decode (the lookup fails: `Unknown`): `psfu_100009` without
`-parse_seqids` (SFU_STORED; the alias `psfu`) and `psfux` with `-parse_seqids` and `-taxid_map`
(SFUX: an `lcl|Subject_901` record, a titled line followed by an untitled one, whose outfmt 0 heading
NCBI skips, and a control).

With `--k5-only` (`python3 tests/db_search/t5/make_t5.py NCBI_BIN --k5-only`, merge 5 of line RG
item T5c with item K4h) it writes only the taxonomy database `tests/seqdb/taxdb/k5/` of
`tests/seqdb/taxdb/k5.tsv` (the rows of `taxdb.tsv`, 99990 still the last entry, and K5_TAXA: names
with two-byte UTF-8, with control bytes and with a byte that `GenerateDefline` cannot decode in all four
fields), the version 4 protein databases of titles ending in `]` with those taxids: `pk5_<taxid>`
without `-parse_seqids` (K5_STORED; the alias `pk5`), `pk5x` (taxid 5003 alone) and `pk5s` with
`-parse_seqids` and `-taxid_map` (K5S), and the queries `k5_qa` and `k5_qc` (`.faa`, and `.fna` of
the same letters back-translated): q0 without hits, q1 with a title of two-byte UTF-8 (`qa`, NCBI's
outfmt 5 stream leaves its `Iteration` out) or a control byte (`qc`) hitting every record, q2 hitting
only the records of group A (K5_A), so that `BATCH_SIZE=1` with `PRE_FETCH_SEQS_LIMIT` of that count
turns q1's prefetch off and q2's on. makeblastdb stores a title byte outside 0x20..0x7E as `#`, so
these bytes reach a database title only through the name at its rewritten end.

With `--header-only` (`python3 tests/db_search/t5/make_t5.py NCBI_BIN --header-only`, line RG item
T5d) it writes the version 4 databases `phd_<name>` and `nhd_<name>` of HEADER_EDITS (`-taxid
9606`; the three records of `hd.faa` / `hd.fna`, the nucleotides the codons of the proteins, so
that every program searches them with `hd_q.faa` / `hd_q.fna`, whose query n hits record n): the
header record of OID 1 (OID 0 for the `_z` databases, which the queries `hd_q1.*` of record 1
alone never report) edited where NCBI's read of it fails (a taxid tag, data where the taxid's end
is expected, the taxid's tag primitive, 2^31 as the taxid, the title's tag, the record cut before
its last four bytes; with `-parse_seqids` and `lcl|hd<n>` ids, the tag of the first string of
the Seq-ids; a member tag in the long form of no member) or where NCBI reads it (the title or
that string a UTF8String, with NCBI's warning; the title's member tag in the long form), and the
record cut after its first byte (`eof`, which the searches end on with exit status 0). NCBI's
reports of them: `remaining-gaps-20261112/T/SURVEY_T5d.md` and `T/fix1_t5d.md` of line RG. The
names after `--header-only` build only those edits (the other databases keep their dates).

With `--fetch-only` (`python3 tests/db_search/t5/make_t5.py NCBI_BIN --fetch-only`, fix 2 of line RG item
T5d) it writes the databases of NCBI's other reads of a UTF8String header (`fetch_main`): the fetches of the
traceback and of the composition-based redo under a negative seqid list, and blastdbcmd's version check.

With `--dup-only` (fix 3 of line RG item T5d) it writes the databases of NCBI's scope's read of another
sequence that holds a reported sequence's Seq-align id (`dup_main`); with `--dup5-only` (fix 5) those of the
order of such a load to a lookup that throws on its id and to another such load (`dup5_main`). With
`--sw-only` (merge 6 of T5d) it writes the databases of the reads of the composition-based redo's fetches
before its Smith-Waterman redo fails a batch (`sw_main`).

Python 3 standard library only; makeblastdb is NCBI's (the comparison oracle builds the inputs).
"""
import os
import struct
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
    # Line RG item T7: the name NCBI cannot decode before and after a named sequence (rows db.t7.*).
    "nt7_bf": ("n", [100009, 100001]),
    "nt7_fb": ("n", [100001, 100009]),
    "pt7_bf": ("p", [100009, 100001]),
    "pt7_fb": ("p", [100001, 100009]),
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


# The records of `psfx` (line RG item T5c): [(Seq-id, title, taxid or None), ...]. The probe of
# T5c (`remaining-gaps-20261112/T/t5c/make_probe.py`) with three organelle names more.
SUFFIX = [
    [("gb|AC400001.1|", "alpha, partial [Bacillus sp.]", None), ("gb|AD400001.1|", "", 100002)],
    [("gb|AC400002.1|", "alpha, partial [Bacillus sp.]", 100002), ("gb|AD400002.1|", "", None)],
    [("gb|AC400003.1|", "beta [bacillus SP.]", 100002), ("gb|AD400003.1|", "", None)],
    [("gb|AC400004.1|", "gamma (mitochondrion) [Bacillus sp.]", 100002), ("gb|AD400004.1|", "", 100002)],
    [("gb|AC400005.1|", "delta (foo) [Bacillus sp.]", 100002), ("gb|AD400005.1|", "", None)],
    [("gb|AC400006.1|", "epsilon [Homo sapiens]", 100002), ("gb|AD400006.1|", "", None)],
    [("gb|AC400007.1|", "zeta protein, partial [ ]", 100002), ("gb|AD400007.1|", "", None)],
    [("gb|AC400008.1|", "eta protein, partial []", 100002), ("gb|AD400008.1|", "", None)],
    [("gb|AC400009.1|", "theta. [Bacillus sp.]", 100002), ("gb|AD400009.1|", "", None)],
    [("gb|AC400010.1|", "TPA: iota, partial [Bacillus sp.].", 100002), ("gb|AD400010.1|", "", None)],
    [("gb|AC400011.1|", "kappa, partial (chloroplast) [Bacillus sp.]", 100002), ("gb|AD400011.1|", "", None)],
    [("gb|AC400012.1|", "lambda (chloroplast), partial [Bacillus sp.]", 100002), ("gb|AD400012.1|", "", None)],
    [("gb|AC400013.1|", "mu(plasmid) [Bacillus sp.]", 100002), ("gb|AD400013.1|", "", None)],
    [("gb|AC400014.1|", "nu (x)(virus) [Bacillus sp.]", 100002), ("gb|AD400014.1|", "", None)],
    [("gb|AC400015.1|", "x, [Bacillus sp.]", 100002), ("gb|AD400015.1|", "", None)],
    [("gb|AC400016.1|", "omicron, partial  [Bacillus sp.]", 100002), ("gb|AD400016.1|", "", None)],
    [("gb|AC400017.1|", "pi, partial[Bacillus sp.]", 100002), ("gb|AD400017.1|", "", None)],
    [("gb|AC400018.1|", "rho, partial, partial [Bacillus sp.]", 100002), ("gb|AD400018.1|", "", None)],
    [("gb|AC400019.1|", "sigma &amp; tau, partial [Bacillus sp.]", 100002), ("gb|AD400019.1|", "", None)],
    [("gb|AC400020.1|", "upsilon, partial [Bacillus sp.] ~;", 100002), ("gb|AD400020.1|", "", None)],
    [("gb|AC400021.1|", "phi, partial [Bacillus sp.]", 77777), ("gb|AD400021.1|", "", None)],
    [("gb|AC400022.1|", "chi protein, partial [ ]", 77777), ("gb|AD400022.1|", "", None)],
    [("gb|AC400023.1|", "psi protein, partial []", 77777), ("gb|AD400023.1|", "", None)],
    [("gb|AC400024.1|", "omega, partial [uncultured bacterium]", None), ("gb|AD400024.1|", "", 100001)],
    [("gb|AC400025.1|", "aa, partial [a &amp; b &lt;x&gt; &#65;]", 100004), ("gb|AD400025.1|", "", None)],
    [("gb|AC400026.1|", "bb, partial [  spaced   NAME  ]", 100003), ("gb|AD400026.1|", "", None)],
    [("gb|AC400027.1|", "cc, partial [TPA: tpa organism]", 100005), ("gb|AD400027.1|", "", None)],
    [("gb|AC400028.1|", "ee, partial [only common]", 100010), ("gb|AD400028.1|", "", None)],
    [("gb|AC400029.1|", "ff, partial [Bacillus sp.]", 100002), ("gb|AD400029.1|", "", 100001),
     ("gb|AE400029.1|", "other title", 100001)],
    [("gb|AC400030.1|", ",, [Bacillus sp.]", 100002), ("gb|AD400030.1|", "", None)],
    [("gb|AC400031.1|", "gg, partial [Bacillus sp.]", None), ("gb|AD400031.1|", "", None)],
    [("gb|AC400032.1|", "hh (mitochondrion), partial [Bacillus sp.]", 100002), ("gb|AD400032.1|", "", None)],
    [("gb|AC400033.1|", "ii protein (plastid) [ ]", 100002), ("gb|AD400033.1|", "", None)],
    [("gb|AC400035.1|", "MAG: kk, partial [Bacillus sp.]", 100002), ("gb|AD400035.1|", "", None)],
    [("gb|AC400036.1|", "ll (chromosome) [Bacillus sp.]", 100002), ("gb|AD400036.1|", "", None)],
    [("gb|AC400037.1|", "mm (protoplast) [Bacillus sp.]", 100002), ("gb|AD400037.1|", "", None)],
    [("gb|AC400038.1|", "nn (extrachromosomal) [Bacillus sp.]", 100002), ("gb|AD400038.1|", "", None)],
    [("lcl|Subject_901", "alpha, partial [Bacillus sp.]", 100002)],
    [("lcl|Subject_902", "x, [Bacillus sp.]", 100002)],
    [("lcl|Subject_903", "a, [ ]", 77777)],
    [("lcl|Subject_904", "theta. [Bacillus sp.]", 100002)],
    [("lcl|Subject_905", "beta, partial [Bacillus sp.]", None), ("gb|AF400005.1|", "", 100002)],
    [("lcl|Subject_906", "kappa, partial (chloroplast) [Bacillus sp.]", 100002)],
]
# The titles of `psfs_<taxid>` (version 4, without -parse_seqids; None: no -taxid).
SUFFIX_STORED = {
    100002: ["alpha, partial [Bacillus sp.]", "x, [Bacillus sp.]", "theta. [Bacillus sp.]",
             "zeta protein, partial [ ]", "kappa, partial (chloroplast) [Bacillus sp.]",
             "TPA: iota, partial [Bacillus sp.].", "beta [bacillus SP.]", "omicron, partial  [Bacillus sp.]",
             "epsilon, [Homo sapiens]", "sigma &amp; tau, partial [Bacillus sp.]", "ab, [Bacillus sp.]",
             "pi, partial[Bacillus sp.]", "nu (x)(virus) [Bacillus sp.]", "ll (chromosome) [Bacillus sp.]",
             "mm, partial (protoplast) [Bacillus sp.]"],
    77777: ["a, [ ]", "phi, partial [Bacillus sp.]", "psi protein, partial []", "chi, partial [ ]"],
    None: ["alpha, partial [Bacillus sp.]", "chi, partial [ ]"],
}

# The titles of `pword_<taxid>` (version 4, without -parse_seqids; None: no -taxid).
WORD_STORED = {
    77777: ["MAG: (plastid) [ ]", " ", "TPA:", "chi, partial [ ]", "  MAG:  "],
    None: ["MAG: (plastid) [ ]", " ", "chi, partial [ ]"],
}

# The titles of `psfu_<taxid>` (merge 4 of line RG item T5c: the rewritten end takes a name that
# NCBI cannot decode; `kappa, partial [Homo sapiens]` keeps its end).
SFU_STORED = {
    100009: ["chi, partial [ ]", "hypothetical protein [ ]", "kappa, partial [Homo sapiens]"],
}
# The records of `psfux` (version 4, -parse_seqids, -taxid_map), as SUFFIX.
SFUX = [
    [("lcl|Subject_901", "chi, partial [ ]", 100009)],
    [("gb|AC400031.1|", "omega, partial [ ]", None), ("gb|AD400031.1|", "", 100009)],
    [("lcl|Subject_904", "delta, partial [Homo sapiens]", 100001)],
]

# The inner bytes of the taxid `[2]` of the first record of `ptx_<name>` and `ntx_<name>` (fix 2 of
# line RG item T5c); makeblastdb writes `02 02 25 86` (INTEGER 9606).
TAXID_ENCODINGS = {
    # The application INTEGER tag, which NCBI's ExpectIntegerTag takes as the INTEGER tag.
    "app42": "42022586",
    # Nine content bytes, the first seven sign extension, which ReadStdSigned skips.
    "sx9": "0209" + "00" * 7 + "2586",
    # 2^31 in five bytes: outside an int, NCBI's "overflow error".
    "big5": "02050080000000",
}
TAXID_TITLES = ["patched record, its taxid in another encoding", "plain record as makeblastdb writes it"]

# The names that `--k5-only` adds to `taxdb.tsv`'s rows for `tests/seqdb/taxdb/k5` (merge 5 of line RG
# item T5c with item K4h): (taxid, scientific, common, blast, kingdom), below 99990 so that it stays
# the last entry (`offset_last_past_end`'s offset error).
K5_TAXA = [
    (5001, b"Caf\xc3\xa9 sp.", b"caf\xc3\xa9 common", b"b\xc3\xa9blast", b"K\xc3\xa9ngdom"),
    (5002, b"Ctl\x07 name", b"c\x1bcommon", b"b\x7fblast", b"K\x01ngdom"),
    (5003, b"Bad \x81 name", b"x\x81common", b"ybl\x81st", b"Bacteria"),
]
# The records of `pk5_<taxid>`: (title, group); group A holds K5_A, group B K5_B.
K5_STORED = {
    99990: [("chi, partial [ ]", "A"), ("alpha, partial [Homo sapiens]", "A"), ("hypothetical protein [ ]", "A"),
            ("bchi, partial [ ]", "B")],
    5001: [("eta, partial [ ]", "A"), ("MAG: (plastid) [ ]", "A"), ("beta2, partial [ ]", "B")],
    5002: [("theta, partial [ ]", "A"), ("iota, partial [Homo sapiens]", "A"), ("btheta, partial [ ]", "B")],
    77777: [("lambda, partial [ ]", "A"), ("blambda [ ]", "B")],
    9606: [("mu, partial [ ]", "A"), ("nu, partial [Homo sapiens]", "A"), ("bmu, partial [ ]", "B")],
}
K5X_STORED = [("kappa, partial [ ]", "A"), ("zeta, partial [Homo sapiens]", "A"), ("bkappa, partial [ ]", "B")]
# The records of `pk5s`: ([(Seq-id, title, taxid or None)], group).
K5S = [
    ([("lcl|Subject_901", "chi, partial [ ]", 99990)], "A"),
    ([("lcl|Subject_902", "eta, partial [ ]", 5001)], "A"),
    ([("lcl|Subject_903", "theta, partial [ ]", 5002)], "A"),
    ([("lcl|Subject_904", "lambda, partial [ ]", 77777)], "A"),
    ([("lcl|Subject_905", "alpha, partial [Homo sapiens]", 9606)], "A"),
    ([("gb|AC500001.1|", "omega, partial [ ]", None), ("gb|AD500001.1|", "", 99990)], "A"),
    ([("gb|AC500002.1|", "sigma, partial [ ]", None), ("gb|AD500002.1|", "", 5001)], "A"),
    ([("gb|AC500003.1|", "tau, partial [ ]", None), ("gb|AD500003.1|", "", 5002)], "A"),
    ([("gb|AC500004.1|", "upsilon, partial [ ]", None), ("gb|AD500004.1|", "", 77777)], "A"),
    ([("lcl|Subject_911", "bchi, partial [ ]", 99990)], "B"),
    ([("lcl|Subject_912", "btheta, partial [ ]", 5002)], "B"),
    ([("gb|AC500011.1|", "bomega, partial [ ]", None), ("gb|AD500011.1|", "", 5001)], "B"),
]
K5_A = ("MWMEHSVLKANFAVLRKVMDDPKIRDGPMAVMKGWENYCWPANMSWKVLMKYEFSSDIHPCNLWICWIPKEMATPGWWNTDTLHMSEMGLWMHK"
        "YPKCMCVMVRPTDQSLRHWKIMHDANFEKLKIFHYM")
K5_Q2 = ("SVLKANFAVLRKVMDDPKIRDGPMAVMKGWENYCWPANMSWKVLMKYEFSSDIHPCNLWICWIPKEMATPGWWNTDTLHMSEMGLWMHKYPKCM"
         "CVMVRPTDQSLRHWKIMHDANFEKLKIFHYMVGHGRTASEG")
K5_B = ("LPVASKCGEPTKQYEKAIRLGQGDFSFFAAIIGGMNIYIGIQMAPRGFLDNMANDMPMTNGTTGCLAPQAYRPQASCGIEKSPVPVLSEPMCR"
        "DINVPFNLYDMNMGDFMTGCDYQCKPL")
K5_NOHIT = "WWWWWWWWWWHHHHHHHHHHWWWWWWWWWW"
K5_QUERIES = {"qa": b"q1 a\xc3\xa9b z", "qc": b"q1 c\x07d z"}
K5_CODONS = dict(A="GCT", R="CGT", N="AAT", D="GAT", C="TGT", Q="CAA", E="GAA", G="GGT", H="CAT", I="ATT",
                 L="CTG", K="AAA", M="ATG", F="TTT", P="CCG", S="TCT", T="ACT", W="TGG", Y="TAT", V="GTT")


def tlv_end(data, at):
    """The end of the BER value at `at` (indefinite constructed values walked)."""
    at += 1
    first = data[at]
    at += 1
    if first == 0x80:
        while data[at:at + 2] != b"\x00\x00":
            at = tlv_end(data, at)
        return at + 2
    if first < 0x80:
        return at + first
    count = first - 0x80
    return at + count + int.from_bytes(data[at:at + count], "big")


def with_taxid(record, inner):
    """`record` (a `Blast-def-line-set` of one line, as makeblastdb writes it) with the inner bytes of
    its indefinite `taxid [2]` replaced by `inner`."""
    assert record[:4] == b"\x30\x80\x30\x80", record[:4].hex()
    at = 4
    while record[at] != 0xA2:
        at = tlv_end(record, at)
    end = tlv_end(record, at)
    assert record[at:at + 2] == b"\xa2\x80" and record[end - 2:end] == b"\x00\x00", record[at:end].hex()
    return record[:at + 2] + bytes.fromhex(inner) + record[end - 2:]


def patch_first_taxid(base, ext, inner):
    """Rewrites the header file of the version 4 database `base` with the first record's taxid
    `inner` (`with_taxid`) and the header offsets of its index file."""
    index = bytearray(open(f"{base}.{ext}in", "rb").read())
    version, = struct.unpack(">I", index[0:4])
    assert version == 4, version
    at = 8
    for _ in range(2):  # title, date
        at += 4 + struct.unpack(">I", index[at:at + 4])[0]
    count, = struct.unpack(">I", index[at:at + 4])
    at += 4 + 8 + 4  # count, total length, longest sequence
    offsets = list(struct.unpack(f">{count + 1}I", index[at:at + 4 * (count + 1)]))
    headers = open(f"{base}.{ext}hr", "rb").read()
    records = [headers[offsets[k]:offsets[k + 1]] for k in range(count)]
    records[0] = with_taxid(records[0], inner)
    offsets = [0]
    for record in records:
        offsets.append(offsets[-1] + len(record))
    index[at:at + 4 * (count + 1)] = struct.pack(f">{count + 1}I", *offsets)
    open(f"{base}.{ext}hr", "wb").write(b"".join(records))
    open(f"{base}.{ext}in", "wb").write(bytes(index))


def taxid_main(ncbi_bin):
    here = os.path.join("tests", "db_search", "t5")
    work = os.path.join(here, "work")
    os.makedirs(work, exist_ok=True)
    env = {"PATH": "/usr/bin:/bin", "HOME": os.path.abspath(work), "BLAST_USAGE_REPORT": "false"}
    for kind, dbtype, ext, source in (("p", "prot", "p", "pdup.faa"), ("n", "nucl", "n", "ndup.fna")):
        seq = letters(os.path.join("tests", "db_search", "ids", source))
        fasta = os.path.join(work, f"{kind}tx.fa")
        with open(fasta, "w") as f:
            for title in TAXID_TITLES:
                f.write(">" + title + "\n" + seq + "\n")
        for name, inner in TAXID_ENCODINGS.items():
            out = os.path.join(here, f"{kind}tx_{name}")
            subprocess.run([os.path.join(ncbi_bin, "makeblastdb"), "-in", fasta, "-dbtype", dbtype,
                            "-blastdb_version", "4", "-taxid", "9606", "-title", f"{kind}tx_{name}",
                            "-out", out], check=True, stdout=subprocess.DEVNULL, env=env)
            # makeblastdb's version 4 JSON metadata file is not read by a search.
            if os.path.exists(f"{out}.{ext}js"):
                os.remove(f"{out}.{ext}js")
            patch_first_taxid(out, ext, inner)
    for name in os.listdir(work):
        os.remove(os.path.join(work, name))
    os.rmdir(work)


# The edits of `phd_<name>` / `nhd_<name>` (line RG item T5d): each takes the header record (a
# `Blast-def-line-set` of one line as makeblastdb writes it) and gives the edited record; `parsed`
# databases are made with `-parse_seqids`.
def with_first_tag(record, inner_of, old, new):
    """`record` with the tag byte `old` that first follows the member tag `inner_of` set to `new`."""
    at = record.index(inner_of)
    at = record.index(bytes([old]), at + len(inner_of))
    return record[:at] + bytes([new]) + record[at + 1:]


def with_long_tag(record, member, number):
    """`record` with the first member tag `member` in the long form: 0xBF, then `number` (one byte)."""
    at = record.index(member)
    return record[:at] + b"\xbf" + bytes([number]) + record[at + 1:]


def with_primitive_taxid(record):
    """`record` with the taxid `[2]` as the primitive tag 0x82 before the INTEGER."""
    at = 4
    while record[at] != 0xA2:
        at = tlv_end(record, at)
    end = tlv_end(record, at)
    return record[:at] + b"\x82" + record[at + 2:end - 2] + record[end:]


HEADER_EDITS = {
    # NCBI: `unexpected tag: VisibleString (26), should be: Integer (2)`, `error while reading taxid`.
    "vstr": (False, lambda r: with_taxid(r, "1a024142")),
    # NCBI: `EndOfContentsByte expected` at the taxid's end (one exception).
    "extra": (False, lambda r: with_taxid(r, "02022586020105")),
    # NCBI: `unexpected tag: contextspecific/2 (130), should be contextspecific/constructed/0 (160)`.
    "prim": (False, with_primitive_taxid),
    # NCBI: `overflow error` (eOverflow).
    "big5": (False, lambda r: with_taxid(r, "02050080000000")),
    # NCBI: `unexpected tag: OctetString (4), should be: VisibleString (26)`, `error while reading title`.
    "title": (False, lambda r: with_first_tag(r, b"\xa0\x80", 0x1A, 0x04)),
    # NCBI: the end of the record in the taxid's member frame (`HandleEOF`, eEOF).
    "cut": (False, lambda r: r[:-4]),
    # NCBI: the string of `local str` (`( at [].[].local.str)`), `error while reading seqid`.
    "seqid": (True, lambda r: with_first_tag(r, b"\xa1\x80", 0x1A, 0x04)),
    # Records that NCBI reads (RG-T5d fix 1): the title a UTF8String, which NCBI reads as the
    # VisibleString with `Warning: (803.10) [<prog>] CObjectIStreamAsnBinary: UTF8String data for
    # VisibleString member Blast-def-line-set.E.title, ...` once per process (none in the maskers);
    "utf8": (False, lambda r: with_first_tag(r, b"\xa0\x80", 0x1A, 0x0C)),
    # the string of `local str` a UTF8String (`Blast-def-line-set.E.seqid.E.local.str`);
    "utf8id": (True, lambda r: with_first_tag(r, b"\xa1\x80", 0x1A, 0x0C)),
    # the title's member tag in the long form (`bf 00`, NCBI's `PeekLongTag`), with no message.
    "long": (False, lambda r: with_long_tag(r, b"\xa0\x80", 0x00)),
    # NCBI: a member tag in the long form whose number is no member's (`bf 20`: `unexpected
    # member: [32]` at the tag's first byte).
    "longx": (False, lambda r: with_long_tag(r, b"\xa0\x80", 0x20)),
    # NCBI: the record cut after its first byte, before the set's elements: the `CEofException`
    # that `HandleEOF` rethrows, which the search programs' `CATCH_ALL` does not post (exit 0).
    "eof": (False, lambda r: r[:1]),
}
# The edits of OID 0 too (`<kind>hd_<name>_z`).
HEADER_EDITS_Z = ("title", "utf8")
HEADER_SEEDS = (11, 22, 33)


def header_letters():
    """The three proteins of `hd.faa` (120 letters each, a linear congruential stream of each seed)
    and their nucleotides, one codon per letter."""
    amino = "ACDEFGHIKLMNPQRSTVWY"
    codon = dict(zip(amino, ["GCT", "TGT", "GAT", "GAA", "TTT", "GGT", "CAT", "ATT", "AAA", "CTG",
                             "ATG", "AAT", "CCG", "CAG", "CGT", "AGC", "ACC", "GTG", "TGG", "TAT"]))
    proteins = []
    for seed in HEADER_SEEDS:
        state, letters_ = seed, []
        for _ in range(120):
            state = (state * 1103515245 + 12345) % 2 ** 31
            letters_.append(amino[(state >> 16) % 20])
        proteins.append("".join(letters_))
    return proteins, ["".join(codon[c] for c in p) for p in proteins]


def patch_record(base, ext, oid, edit):
    """Rewrites the header file of the version 4 database `base` with the record of `oid` edited
    by `edit`, and the header offsets of its index file."""
    index = bytearray(open(f"{base}.{ext}in", "rb").read())
    version, = struct.unpack(">I", index[0:4])
    assert version == 4, version
    at = 8
    for _ in range(2):  # title, date
        at += 4 + struct.unpack(">I", index[at:at + 4])[0]
    count, = struct.unpack(">I", index[at:at + 4])
    at += 4 + 8 + 4  # count, total length, longest sequence
    offsets = list(struct.unpack(f">{count + 1}I", index[at:at + 4 * (count + 1)]))
    headers = open(f"{base}.{ext}hr", "rb").read()
    records = [headers[offsets[k]:offsets[k + 1]] for k in range(count)]
    records[oid] = edit(records[oid])
    offsets = [0]
    for record in records:
        offsets.append(offsets[-1] + len(record))
    index[at:at + 4 * (count + 1)] = struct.pack(f">{count + 1}I", *offsets)
    open(f"{base}.{ext}hr", "wb").write(b"".join(records))
    open(f"{base}.{ext}in", "wb").write(bytes(index))


def header_main(ncbi_bin, only=()):
    here = os.path.join("tests", "db_search", "t5")
    work = os.path.join(here, "work")
    os.makedirs(work, exist_ok=True)
    proteins, nucleotides = header_letters()
    for kind, seqs in (("p", proteins), ("n", nucleotides)):
        ext = "faa" if kind == "p" else "fna"
        with open(os.path.join(here, f"hd.{ext}"), "w") as f:
            for n, seq in enumerate(seqs):
                f.write(f">hd{n} record {n} of the T5d databases\n{seq}\n")
        with open(os.path.join(here, f"hd_q.{ext}"), "w") as f:
            for n, seq in enumerate(seqs):
                f.write(f">q{n}\n{seq}\n")
        with open(os.path.join(here, f"hd_q1.{ext}"), "w") as f:
            f.write(f">q1\n{seqs[1]}\n")
        with open(os.path.join(work, f"hdp.{ext}"), "w") as f:
            for n, seq in enumerate(seqs):
                f.write(f">lcl|hd{n} record {n} of the T5d databases\n{seq}\n")
    for kind, dbtype, ext in (("p", "prot", "p"), ("n", "nucl", "n")):
        fext = "faa" if kind == "p" else "fna"
        for name, (parsed, edit) in HEADER_EDITS.items():
            if only and name not in only:
                continue
            for suffix, oid in (("", 1), ("_z", 0)):
                if suffix and name not in HEADER_EDITS_Z:
                    continue
                fasta = os.path.join(work, f"hdp.{fext}") if parsed else os.path.join(here, f"hd.{fext}")
                out = os.path.join(here, f"{kind}hd_{name}{suffix}")
                makeblastdb(ncbi_bin, fasta, dbtype, 9606, out, parsed)
                patch_record(out, ext, oid, edit)
    for name in os.listdir(work):
        os.remove(os.path.join(work, name))
    os.rmdir(work)


def patch_record_any(base, ext, oid, edit):
    """`patch_record` of a database of either version (version 5's index holds the LMDB file's name
    after the title)."""
    index = bytearray(open(f"{base}.{ext}in", "rb").read())
    version, = struct.unpack(">I", index[0:4])
    at = 8 + (4 if version == 5 else 0)
    for _ in range(3 if version == 5 else 2):  # title, (LMDB name,) date
        at += 4 + struct.unpack(">I", index[at:at + 4])[0]
    count, = struct.unpack(">I", index[at:at + 4])
    at += 4 + 8 + 4
    offsets = list(struct.unpack(f">{count + 1}I", index[at:at + 4 * (count + 1)]))
    headers = open(f"{base}.{ext}hr", "rb").read()
    records = [headers[offsets[k]:offsets[k + 1]] for k in range(count)]
    records[oid] = edit(records[oid])
    offsets = [0]
    for record in records:
        offsets.append(offsets[-1] + len(record))
    index[at:at + 4 * (count + 1)] = struct.pack(f">{count + 1}I", *offsets)
    open(f"{base}.{ext}hr", "wb").write(b"".join(records))
    open(f"{base}.{ext}in", "wb").write(bytes(index))


def fetch_main(ncbi_bin):
    """`--fetch-only` (line RG item T5d fix 2): the databases of NCBI's other reads of a header with a
    UTF8String title. The version 5 databases `pfx_utf8`, `pfx_utf8r` and `nfx_utf8` with `-parse_seqids`
    (`NSW0.1` the protein hd0 of `hd.faa`, which the query `hd_q0.faa` matches, `NSW1.1` hd2, whose chance
    alignment with it has an E-value of about 6, with the UTF8String title, and `NSW2.1` hd1; `pfx_utf8r` the
    first two in the other order; `nfx_utf8` their codons as `NSWN<n>.1`; `pfx_copy` and `nfx_copy` with
    `NSC1.1` / `NSCN1.1` a weaker copy of hd0 instead, every sixth letter changed, which the redo fetches),
    and the negative seqid lists `fx_neg1.txt` (the weak sequence) and `fx_neg2.txt` (another; `nfx_*`,
    `cfx_*` and `ncfx_*` of the other databases): under them NCBI's
    traceback and composition-based redo read the header of each subject they fetch, the redo only those its
    early termination test leaves (`BlastCompo_EarlyTermination`: with `-max_target_seqs 1` the strong match
    fills the heap, and the weak one is not fetched where the redo takes it after). The version 4 database
    `pva_utf8` with `-parse_seqids` (`gb|HDA<n>.1|`, the UTF8String title at OID 1) and `pva_batch.txt`:
    blastdbcmd's version check of `HDA1.2` reads HDA1.1's header."""
    here = os.path.join("tests", "db_search", "t5")
    work = os.path.join(here, "work")
    os.makedirs(work, exist_ok=True)
    proteins, nucleotides = header_letters()
    for ext, seqs in (("faa", proteins), ("fna", nucleotides)):
        with open(os.path.join(here, f"hd_q0.{ext}"), "w") as f:
            f.write(f">q0\n{seqs[0]}\n")
    utf8 = lambda r: with_first_tag(r, b"\xa0\x80", 0x1A, 0x0C)
    amino = "ACDEFGHIKLMNPQRSTVWY"
    # A weaker copy of hd0: every sixth letter the next of the alphabet (83 % identity), and its codons.
    copy = "".join(amino[(amino.index(c) + 1) % 20] if i % 6 == 3 else c for i, c in enumerate(proteins[0]))
    codons = {c: nucleotides[0][3 * i:3 * i + 3] for i, c in enumerate(proteins[0])}
    for k, c in enumerate(proteins[1] + proteins[2]):
        codons.setdefault(c, (nucleotides[1] + nucleotides[2])[3 * k:3 * k + 3])
    copied = "".join(codons[c] for c in copy)
    for name, kind, order, prefix, weak in (("pfx_utf8", "p", (0, 2, 1), "NSW", None),
                                             ("pfx_utf8r", "p", (2, 0, 1), "NSW", None),
                                             ("nfx_utf8", "n", (0, 2, 1), "NSWN", None),
                                             ("pfx_copy", "p", (0, 2, 1), "NSC", copy),
                                             ("nfx_copy", "n", (0, 2, 1), "NSCN", copied)):
        seqs = list(proteins if kind == "p" else nucleotides)
        if weak is not None:
            seqs[2] = weak
        fasta = os.path.join(work, f"{name}.fa")
        names = {0: "0.1 strong", 2: "1.1 weak", 1: "2.1 other"}
        with open(fasta, "w") as f:
            for k in order:
                f.write(f">{prefix}{names[k]}\n{seqs[k]}\n")
        out = os.path.join(here, name)
        cmd = [os.path.join(ncbi_bin, "makeblastdb"), "-in", fasta, "-dbtype", "prot" if kind == "p" else "nucl",
               "-blastdb_version", "5", "-taxid", "9606", "-title", name, "-out", out, "-parse_seqids"]
        env = {"PATH": "/usr/bin:/bin", "HOME": os.path.abspath(here), "BLAST_USAGE_REPORT": "false"}
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, env=env)
        patch_record_any(out, kind, order.index(2), utf8)
    for kind, prefix in (("", "NSW"), ("n", "NSWN"), ("c", "NSC"), ("nc", "NSCN")):
        for k in (1, 2):
            with open(os.path.join(here, f"{kind}fx_neg{k}.txt"), "w") as f:
                f.write(f"{prefix}{k}.1\n")
    fasta = os.path.join(work, "pva.fa")
    with open(fasta, "w") as f:
        for n, seq in enumerate(proteins):
            f.write(f">gb|HDA{n}.1| record {n} of the T5d databases\n{seq}\n")
    out = os.path.join(here, "pva_utf8")
    makeblastdb(ncbi_bin, fasta, "prot", 9606, out, True)
    patch_record(out, "p", 1, utf8)
    with open(os.path.join(here, "pva_batch.txt"), "w") as f:
        f.write("HDA1.2\nHDA2.1\n")
    for name in os.listdir(work):
        os.remove(os.path.join(work, name))
    os.rmdir(work)


def sw_main(ncbi_bin):
    """`--sw-only` (merge 6 of line RG item T5d, with the rows of line RG item X4-excl): the version 5
    databases `psw_utf8` and `psw_utf8x` with `-parse_seqids` of the records of `../x4excl/rev`
    (sx_strong, `gnl|ti|7100`, which `../x4excl/neg.txt` keeps without a definition line, so the fetch of
    the composition-based redo excludes it, then sfail, whose Smith-Waterman redo against tks_W of
    `tests/fasta/rps/tks_w.faa` fails at the batch's last match), the title of sfail's record (`psw_utf8`,
    OID 1) or of sx_strong's (`psw_utf8x`, OID 0) a UTF8String: under the list NCBI's redo reads the header
    of each subject it fetches, both before its failure (`BLAST engine error: Traceback failed`), so its
    warning comes before its error."""
    here = os.path.join("tests", "db_search", "t5")
    work = os.path.join(here, "work")
    os.makedirs(work, exist_ok=True)

    def records_of(path):
        seqs, name = {}, None
        for line in open(path):
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = ""
            elif name:
                seqs[name] += line
        return seqs

    tks_w = records_of(os.path.join("tests", "fasta", "rps", "tks_w.faa"))["tks_W"]
    sfail = records_of(os.path.join("tests", "fasta", "x4", "sw_subjects.faa"))["sfail"]
    utf8 = lambda r: with_first_tag(r, b"\xa0\x80", 0x1A, 0x0C)
    for name, oid in (("psw_utf8", 1), ("psw_utf8x", 0)):
        fasta = os.path.join(work, f"{name}.fa")
        with open(fasta, "w") as f:
            f.write(f">gnl|ti|7100 sx_strong tks_W itself\n{tks_w}\n")
            f.write(f">dbj|BAA200102.1| sfail tks_M with D first, twice\n{sfail}\n")
        out = os.path.join(here, name)
        cmd = [os.path.join(ncbi_bin, "makeblastdb"), "-in", fasta, "-dbtype", "prot", "-parse_seqids",
               "-blastdb_version", "5", "-title", name, "-out", out]
        env = {"PATH": "/usr/bin:/bin", "HOME": os.path.abspath(here), "BLAST_USAGE_REPORT": "false"}
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, env=env)
        patch_record_any(out, "p", oid, utf8)
    for name in os.listdir(work):
        os.remove(os.path.join(work, name))
    os.rmdir(work)


def dup_main(ncbi_bin):
    """`--dup-only` (line RG item T5d fix 3): databases whose lists hold `gb|DUP00001.1|` twice. `pds` and
    `pda` (version 5, protein, `-parse_seqids`): the subject `gb|DUP00001.1| subject record`, hd0 of `hd.faa`
    (which `hd_q0.faa` matches); `pdo_same`: that record with its title a UTF8String; `pdo_none`:
    `gb|DUP00001.1| other record`, hd1 (no match), with its title a UTF8String. NCBI orders a list's volumes by
    path, so `pdo_* pds` puts the other sequence at OID 0, whose header NCBI's scope reads where it resolves the
    subject's Seq-align id (the data loader's `SeqidToOid` gives it first), and `pda pdo_*` the subject first.
    `nds`, `ndo_same` and `ndo_none`: the same with their codons (nucleotide); `pds4` and `pdo4_same`: the
    protein ones of version 4. `dup_batch.txt`: `DUP00001.1` (blastdbcmd `-entry_batch`)."""
    here = os.path.join("tests", "db_search", "t5")
    work = os.path.join(here, "work")
    os.makedirs(work, exist_ok=True)
    proteins, nucleotides = header_letters()
    utf8 = lambda r: with_first_tag(r, b"\xa0\x80", 0x1A, 0x0C)
    subject = ("gb|DUP00001.1| subject record", 0)
    for name, (title, letters), version, edit in (
            ("pds", subject, 5, None), ("pda", subject, 5, None), ("pdo_same", subject, 5, utf8),
            ("pdo_none", ("gb|DUP00001.1| other record", 1), 5, utf8),
            ("nds", subject, 5, None), ("ndo_same", subject, 5, utf8),
            ("ndo_none", ("gb|DUP00001.1| other record", 1), 5, utf8),
            ("pds4", subject, 4, None), ("pdo4_same", subject, 4, utf8)):
        kind = name[0]
        fasta = os.path.join(work, f"{name}.fa")
        with open(fasta, "w") as f:
            f.write(f">{title}\n{(proteins if kind == 'p' else nucleotides)[letters]}\n")
        out = os.path.join(here, name)
        if version == 4:
            makeblastdb(ncbi_bin, fasta, "prot", 9606, out, True)
        else:
            cmd = [os.path.join(ncbi_bin, "makeblastdb"), "-in", fasta, "-dbtype", "prot" if kind == "p" else "nucl",
                   "-blastdb_version", "5", "-taxid", "9606", "-title", name, "-out", out, "-parse_seqids"]
            env = {"PATH": "/usr/bin:/bin", "HOME": os.path.abspath(here), "BLAST_USAGE_REPORT": "false"}
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, env=env)
        if edit is not None:
            patch_record_any(out, kind, 0, edit)
    with open(os.path.join(here, "dup_batch.txt"), "w") as f:
        f.write("DUP00001.1\n")
    for name in os.listdir(work):
        os.remove(os.path.join(work, name))
    os.rmdir(work)


def dup5_main(ncbi_bin):
    """`--dup5-only` (line RG item T5d fix 5): version 5 databases (`-parse_seqids`, taxid 9606) whose
    lists hold a reported sequence's Seq-align id twice, the copy first (NCBI orders a list's volumes by
    path), beside a sequence whose id NCBI's lookup throws on (`SeqDB_SimplifyAccession` of a `PIG` string
    tag; NCBI's prefetch catches it with `Error pre-fetching sequence data`). `pdo5_same` and `pds5`:
    `gb|DUP00005.1| subject record`, hd0 (which `q0` matches), the copy with its title a UTF8String;
    `ppig1`: `gnl|PIG|x001 throwing record`, hd1 (`q1`), whose Seq-id handle (packed number 1) NCBI's
    prefetch loads before DUP00005's (5); `ppig5`: `gnl|PIG|x005 throwing tie`, hd1, one packed number
    with DUP00005's (NCBI orders the two handles by the addresses of their infos); `pdo5b_same`,
    `pdo5b_id` and `pds5b`: `gb|DUQ00005.1| subject record`, hd2 (`q2`), packed 5 too, the copies with
    their title or their accession a UTF8String; `pgo5` and `pgs5`: `gi|5 gi subject record`, hd0, the
    copy's title a UTF8String (a gi's lookup reads no accession table; packed 5 too), and `pplain` (hd2,
    without `-parse_seqids`), whose list's lookup of a text id stops NCBI's prefetch (`DB contains no
    accession info`). `nds5`, `ndo5_same` and `npig1`: the same with their codons (nucleotide). Queries `hd_q10.faa` (`q1`, `q0`, `q0b`: hd1, hd0, hd0) and `hd_q02.faa` (`q0`,
    `q2`), and their codons (`.fna`)."""
    here = os.path.join("tests", "db_search", "t5")
    work = os.path.join(here, "work")
    os.makedirs(work, exist_ok=True)
    proteins, nucleotides = header_letters()
    title = lambda r: with_first_tag(r, b"\xa0\x80", 0x1A, 0x0C)
    accession = lambda r: with_first_tag(r, b"\xa1\x80", 0x1A, 0x0C)
    dup = "gb|DUP00005.1| subject record"
    duq = "gb|DUQ00005.1| subject record"
    for name, header, letters, edit in (
            ("pdo5_same", dup, 0, title), ("pds5", dup, 0, None),
            ("ppig1", "gnl|PIG|x001 throwing record", 1, None), ("ppig5", "gnl|PIG|x005 throwing tie", 1, None),
            ("pdo5b_same", duq, 2, title), ("pdo5b_id", duq, 2, accession), ("pds5b", duq, 2, None),
            ("pgo5", "gi|5 gi subject record", 0, title), ("pgs5", "gi|5 gi subject record", 0, None),
            ("pplain", "plain record", 2, None),
            ("ndo5_same", dup, 0, title), ("nds5", dup, 0, None),
            ("npig1", "gnl|PIG|x001 throwing record", 1, None)):
        kind = name[0]
        fasta = os.path.join(work, f"{name}.fa")
        with open(fasta, "w") as f:
            f.write(f">{header}\n{(proteins if kind == 'p' else nucleotides)[letters]}\n")
        out = os.path.join(here, name)
        cmd = [os.path.join(ncbi_bin, "makeblastdb"), "-in", fasta, "-dbtype", "prot" if kind == "p" else "nucl",
               "-blastdb_version", "5", "-taxid", "9606", "-title", name, "-out", out]
        cmd += [] if name == "pplain" else ["-parse_seqids"]
        env = {"PATH": "/usr/bin:/bin", "HOME": os.path.abspath(here), "BLAST_USAGE_REPORT": "false"}
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, env=env)
        if edit is not None:
            patch_record_any(out, kind, 0, edit)
    for ext, seqs in (("faa", proteins), ("fna", nucleotides)):
        with open(os.path.join(here, f"hd_q10.{ext}"), "w") as f:
            f.write(f">q1\n{seqs[1]}\n>q0\n{seqs[0]}\n>q0b\n{seqs[0]}\n")
        with open(os.path.join(here, f"hd_q02.{ext}"), "w") as f:
            f.write(f">q0\n{seqs[0]}\n>q2\n{seqs[2]}\n")
    for name in os.listdir(work):
        os.remove(os.path.join(work, name))
    os.rmdir(work)


def stored_main(ncbi_bin, prefix, stored):
    """The version 4 protein databases `<prefix>_<taxid>` of the titles `stored` (the letters of
    `pdup.faa`, `-taxid <taxid>`, none for `None`) and the alias `<prefix>` over them."""
    here = os.path.join("tests", "db_search", "t5")
    seq = letters(os.path.join("tests", "db_search", "ids", "pdup.faa"))
    work = os.path.join(here, "work")
    os.makedirs(work, exist_ok=True)
    env = {"PATH": "/usr/bin:/bin", "HOME": os.path.abspath(work), "BLAST_USAGE_REPORT": "false"}
    names = []
    for taxid, titles in stored.items():
        name = f"{prefix}_{taxid if taxid is not None else 'none'}"
        names.append(name)
        fasta = os.path.join(work, name + ".fa")
        with open(fasta, "w") as f:
            for title in titles:
                f.write(">" + title + "\n" + seq + "\n")
        cmd = [os.path.join(ncbi_bin, "makeblastdb"), "-in", fasta, "-dbtype", "prot",
               "-blastdb_version", "4", "-title", name, "-out", os.path.join(here, name)]
        if taxid is not None:
            cmd += ["-taxid", str(taxid)]
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, env=env)
        # makeblastdb's version 4 JSON metadata file is not read by a search.
        if os.path.exists(os.path.join(here, name + ".pjs")):
            os.remove(os.path.join(here, name + ".pjs"))
    with open(os.path.join(here, f"{prefix}.pal"), "w") as f:
        f.write(f"TITLE {prefix}\nDBLIST {' '.join(names)} \n")
    for name in os.listdir(work):
        os.remove(os.path.join(work, name))
    os.rmdir(work)


def suffix_main(ncbi_bin):
    here = os.path.join("tests", "db_search", "t5")
    seq = letters(os.path.join("tests", "db_search", "ids", "pdup.faa"))
    work = os.path.join(here, "work")
    os.makedirs(work, exist_ok=True)
    env = {"PATH": "/usr/bin:/bin", "HOME": os.path.abspath(work), "BLAST_USAGE_REPORT": "false"}
    mk = os.path.join(ncbi_bin, "makeblastdb")
    fasta, taxmap = os.path.join(work, "psfx.fa"), os.path.join(work, "psfx.map")
    with open(fasta, "w") as f, open(taxmap, "w") as m:
        for lines in SUFFIX:
            parts = []
            for seqid, title, taxid in lines:
                parts.append(seqid + (" " + title if title else ""))
                if taxid is not None:
                    m.write(f"{map_key(seqid)}\t{taxid}\n")
            f.write(">" + CTRL_A.join(parts) + "\n" + seq + "\n")
    subprocess.run([mk, "-in", fasta, "-dbtype", "prot", "-parse_seqids", "-blastdb_version", "5",
                    "-taxid_map", taxmap, "-title", "psfx", "-out", os.path.join(here, "psfx")],
                   check=True, stdout=subprocess.DEVNULL, env=env)
    names = []
    for taxid, titles in SUFFIX_STORED.items():
        name = f"psfs_{taxid if taxid is not None else 'none'}"
        names.append(name)
        fasta = os.path.join(work, name + ".fa")
        with open(fasta, "w") as f:
            for title in titles:
                f.write(">" + title + "\n" + seq + "\n")
        cmd = [mk, "-in", fasta, "-dbtype", "prot", "-blastdb_version", "4", "-title", name,
               "-out", os.path.join(here, name)]
        if taxid is not None:
            cmd += ["-taxid", str(taxid)]
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, env=env)
        # makeblastdb's version 4 JSON metadata file is not read by a search.
        if os.path.exists(os.path.join(here, name + ".pjs")):
            os.remove(os.path.join(here, name + ".pjs"))
    with open(os.path.join(here, "psfs.pal"), "w") as f:
        f.write(f"TITLE psfs\nDBLIST {' '.join(names)} \n")
    for name in os.listdir(work):
        os.remove(os.path.join(work, name))
    os.rmdir(work)


def sfu_main(ncbi_bin):
    """`psfu_100009` with its alias `psfu` (`stored_main`) and `psfux` (SFUX)."""
    stored_main(ncbi_bin, "psfu", SFU_STORED)
    here = os.path.join("tests", "db_search", "t5")
    seq = letters(os.path.join("tests", "db_search", "ids", "pdup.faa"))
    work = os.path.join(here, "work")
    os.makedirs(work, exist_ok=True)
    env = {"PATH": "/usr/bin:/bin", "HOME": os.path.abspath(work), "BLAST_USAGE_REPORT": "false"}
    fasta, taxmap = os.path.join(work, "psfux.fa"), os.path.join(work, "psfux.map")
    with open(fasta, "w") as f, open(taxmap, "w") as m:
        for lines in SFUX:
            parts = []
            for seqid, title, taxid in lines:
                parts.append(seqid + (" " + title if title else ""))
                if taxid is not None:
                    m.write(f"{map_key(seqid)}\t{taxid}\n")
            f.write(">" + CTRL_A.join(parts) + "\n" + seq + "\n")
    out = os.path.join(here, "psfux")
    subprocess.run([os.path.join(ncbi_bin, "makeblastdb"), "-in", fasta, "-dbtype", "prot", "-parse_seqids",
                    "-blastdb_version", "4", "-taxid_map", taxmap, "-title", "psfux", "-out", out],
                   check=True, stdout=subprocess.DEVNULL, env=env)
    # makeblastdb's version 4 JSON metadata file is not read by a search.
    if os.path.exists(out + ".pjs"):
        os.remove(out + ".pjs")
    for name in os.listdir(work):
        os.remove(os.path.join(work, name))
    os.rmdir(work)


def k5_main(ncbi_bin):
    """`tests/seqdb/taxdb/k5`, `pk5_<taxid>` with the alias `pk5`, `pk5x`, `pk5s` and the queries
    `k5_q*` (K5_TAXA, K5_STORED, K5X_STORED, K5S, K5_QUERIES)."""
    here = os.path.join("tests", "db_search", "t5")
    taxdir = os.path.join("tests", "seqdb", "taxdb")
    rows = [line for line in open(os.path.join(taxdir, "taxdb.tsv"), "rb").read().split(b"\n") if line]
    rows += [b"\t".join([str(row[0]).encode()] + list(row[1:])) for row in K5_TAXA]
    rows.sort(key=lambda row: int(row.split(b"\t")[0]))
    tsv = os.path.join(taxdir, "k5.tsv")
    with open(tsv, "wb") as f:
        f.write(b"\n".join(rows) + b"\n")
    subprocess.run([sys.executable, os.path.join(taxdir, "make_taxdb.py"), "write", tsv,
                    os.path.join(taxdir, "k5")], check=True)
    work = os.path.join(here, "work")
    os.makedirs(work, exist_ok=True)
    env = {"PATH": "/usr/bin:/bin", "HOME": os.path.abspath(work), "BLAST_USAGE_REPORT": "false"}
    seqs = {"A": K5_A, "B": K5_B}

    def make(name, records, taxid, parse, taxmap=None):
        fasta = os.path.join(work, name + ".fa")
        with open(fasta, "w") as f:
            for defline, group in records:
                f.write(">" + defline + "\n" + seqs[group] + "\n")
        cmd = [os.path.join(ncbi_bin, "makeblastdb"), "-in", fasta, "-dbtype", "prot",
               "-blastdb_version", "4", "-title", name, "-out", os.path.join(here, name)]
        if taxid is not None:
            cmd += ["-taxid", str(taxid)]
        if parse:
            cmd += ["-parse_seqids", "-taxid_map", taxmap]
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, env=env)
        # makeblastdb's version 4 JSON metadata file is not read by a search.
        if os.path.exists(os.path.join(here, name + ".pjs")):
            os.remove(os.path.join(here, name + ".pjs"))

    for taxid, records in K5_STORED.items():
        make(f"pk5_{taxid}", records, taxid, False)
    with open(os.path.join(here, "pk5.pal"), "w") as f:
        f.write(f"TITLE pk5\nDBLIST {' '.join(f'pk5_{taxid}' for taxid in K5_STORED)} \n")
    make("pk5x", K5X_STORED, 5003, False)
    taxmap = os.path.join(work, "pk5s.map")
    records = []
    with open(taxmap, "w") as m:
        for lines, group in K5S:
            records.append((CTRL_A.join(seqid + (" " + title if title else "") for seqid, title, _ in lines),
                            group))
            for seqid, _, taxid in lines:
                if taxid is not None:
                    m.write(f"{map_key(seqid)}\t{taxid}\n")
    make("pk5s", records, None, True, taxmap)
    for name in os.listdir(work):
        os.remove(os.path.join(work, name))
    os.rmdir(work)
    for name, title in K5_QUERIES.items():
        queries = [(b"q0 none", K5_NOHIT), (title, K5_A + K5_B), (b"q2 two", K5_Q2)]
        with open(os.path.join(here, f"k5_{name}.faa"), "wb") as f:
            for t, letters_ in queries:
                f.write(b">" + t + b"\n" + letters_.encode() + b"\n")
        with open(os.path.join(here, f"k5_{name}.fna"), "wb") as f:
            for t, letters_ in queries:
                f.write(b">" + t + b"\n" + "".join(K5_CODONS[c] for c in letters_).encode() + b"\n")


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
    elif "--suffix-only" in sys.argv[2:]:
        suffix_main(sys.argv[1])
    elif "--word-only" in sys.argv[2:]:
        stored_main(sys.argv[1], "pword", WORD_STORED)
    elif "--taxid-only" in sys.argv[2:]:
        taxid_main(sys.argv[1])
    elif "--fetch-only" in sys.argv[2:]:
        fetch_main(sys.argv[1])
    elif "--sw-only" in sys.argv[2:]:
        sw_main(sys.argv[1])
    elif "--dup-only" in sys.argv[2:]:
        dup_main(sys.argv[1])
    elif "--dup5-only" in sys.argv[2:]:
        dup5_main(sys.argv[1])
    elif "--header-only" in sys.argv[2:]:
        header_main(sys.argv[1], [name for name in sys.argv[2:] if not name.startswith("--")])
    elif "--sfu-only" in sys.argv[2:]:
        sfu_main(sys.argv[1])
    elif "--k5-only" in sys.argv[2:]:
        k5_main(sys.argv[1])
    else:
        main()
