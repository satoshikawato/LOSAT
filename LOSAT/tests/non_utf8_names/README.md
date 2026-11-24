# Fixtures of names and words that are not UTF-8 (`LOSAT/tests/non_utf8_names.rs`)

What NCBI BLAST+ 2.17.0 (the release build on the developer machine, `$NCBI_BIN`) wrote for the
cases of `cases.tsv`: a file name of BLASTX's and RPSTBLASTN's `-query`, `-subject` and `-out`
that does not open (`R1*`), BLASTX's `-subject` name in the `Database:` lines of outfmt 7 and 0
(`R2*`), the words of a text seqid list that `CSeq_id` does not parse (`R3*`, `-seqidlist` and
`-negative_seqidlist` on `db_search/q1_sq1`), a text `-mask_data` file of makeblastdb whose
member name is the byte 0xE9 after the 8,191 bytes that the format guess reads (`R4`,
`mask_e9.asnt`; NCBI's `unexpected member` message keeps the byte), `blastdbcmd -list` over a database and a directory
with such a name (`R6*`), a directory of `NCBI_CONFIG_PATH` with such a name whose `.ncbirc`
names `BLASTDB` (`R7`) and a version 5 volume that names a missing LMDB file of such a name
(`R8`, its `.pin` patched) and `blastdbcmd -metadata` of a volume whose LMDB file and its
`.pot`, `.ptf` and `.pto` have such a name (`R8m`; NCBI counts and lists them), `blastdbcmd
-taxids` over such files of `db_search/ids/prot_v5` (`R8t` a taxid of the database, `R8u` one not in
it, whose message names the LMDB file), each name with the
byte 0xE9, and a `-matrix` name of 1,000 `a` and 30
`é` that NCBI's 1,024-byte message buffer cuts inside a character, after its byte 0xC3 (`R9`). Frozen on 2026-11-20 by line LIB-B (`R4` on 2026-11-20 22:56 JST, added 2026-11-21)
(`make_fixtures.py --ncbi-bin "$NCBI_BIN"`, which writes `expected/<case>.stdout`, `.stderr` and
`.rc`, the case's directory as `{D}`); the same cases on the workstation's 2.17.0 gave the same
bytes (line LIB, survey 2 of the API messages).

The `K10_*` cases (line RG item K10a, frozen on 2026-11-24 with the same command, `R6a` again
with them) are the values of dustmasker (`dm`), segmasker (`sm`), windowmasker (`wm`),
makembindex (`mi`), convert2blastmask (`cb`), blastdbcmd (`bd`) and makeprofiledb (`mp`) with the
byte 0xE9: names of input files, databases (`-infmt blastdb`, `-db`, a `-fa_list` line, `-list`),
output files and prefixes, missing ones, `BLASTDB` of the environment and of `.ncbirc`, and words
that name a format, an algorithm, its options, an entry, a range or a mask (a `-taxids` word that is
not a number stays LOSAT's refusal, as `-taxids abc`: NCBI's exception names its own source files). NCBI reads
each as the bytes of its `std::string`; `R6a` (`-list_outfmt '%f %t'`), once LOSAT's refusal, is
now NCBI's bytes too. On the same day every case got `expected/<case>.files`, the names of the
files and directories its run made (one per line in byte order, a byte outside printable ASCII or
the backslash as `\xHH`), which the test compares too; the earlier cases' standard output, error and
exit status came out the same again, but for the time of `R4`'s makeblastdb line and the build path
in its exception message (kept as frozen on 2026-11-20). The `K10_*_gimask_*` cases (frozen the same day) read
GI-indexed masks whose alias file's `MASKLIST` names them with the byte 0xE9 (the database and the
mask files of `LOSAT/tests/seqdb/gi_mask_read/db/sv4`, the masks copied under that name):
blastdbcmd `-info`, `%m0`, `-mask_sequence_with` by that name and by an unknown one (NCBI's table of
names keeps the byte), a `MASKLIST` name without files, makembindex `-show_filters` and `-db_mask`,
and blastn `-db_soft_mask`.
