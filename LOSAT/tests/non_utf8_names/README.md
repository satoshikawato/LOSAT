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
and blastn `-db_soft_mask`. Then (round 2 of the item's review, the same day) each line of `.files` also holds the file's
size and the FNV-1a 64 of its bytes (the case's directory in them as `{D}`), but for the files that
hold their time of writing, whose names end with `.nin`, `.pin`, `.nal`, `.pal`, `.njs` or `.pjs`
(`*`; two freezes a minute apart gave the same `.files` for every case), and the test masks the
seconds of makeblastdb's `added <n> sequences in <seconds> seconds.`; the setup steps gained
`symlink` and `cwd` (the directory where the program runs), and eleven cases: windowmasker
`-exclude_ids` of a file and masking a database (`-infmt blastdb`), makeprofiledb `-taxid_map` of a
file, an `-in` list entry and a run in a directory, makeblastdb in a directory (`K10_mb_cwd`: the
database name is the directory's bytes and `-out`), blastdbcmd `-taxidlist` and `-ipg_batch` of a
file, an alias file's `DBLIST` entry, a database through a symbolic link to a directory, and a
volume without its `.nsq` (NCBI's `File (...) not found` keeps the byte), each name with 0xE9.
Round 3 of the review (the same day, all cases frozen again) replaced that rule: `.files` now lists
every path the run made, changed or removed (a removed one with `gone`), and a file's size and
digest are those of its bytes with the build's dates as `<DATE>` (`undated`): the date field of an
index file's header (`.nin`, `.pin`: its length and its bytes with their padding), every
`MM/DD/YYYY hh:mm:ss` (an alias file's comment, the index file of a column or of a GI mask) and
`YYYY-MM-DDThh:mm:ss` (a metadata file's `last-updated`), and a metadata file's `bytes-total`
counted without the date fields of the index files it lists (their length changes with the date).
So the titles and names inside index, alias and metadata files are compared too (`K10_mp_title`,
`K10_mp_in_file`). Eight cases came with it: `makeblastdb -input_type blastdb` copying the
`MASKLIST` database (its mask column), `-db_hard_mask` by that mask in blastn, tblastn and tblastx
(`sv4`) and blastp (`svp`, the mask `pm1` under the name `p\xe91`), a second makeblastdb in
`d\xe9` over an existing database (`Deleted existing ... named` with the byte, the old volume's
extra files gone), and rpsblast and rpstblastn on a profile database found through `BLASTDB` of
`d\xe9` without its `.aux` (NCBI's engine error names the path with the byte).

The `K10d_*` cases (line RG item K10d, frozen on 2026-11-25 01:56 JST with the same command, and
again at 02:06 with the same bytes) are the values of blastdb_aliastool with the byte 0xE9 (and 0xFF
in a title, an output name and a `-vdblist_file` line), in each of its modes, on the databases and
lists of `LOSAT/tests/fixtures/blastdb_aliastool/in`: `-dblist` (one name, two, a missing one),
`-dblist_file` (its name, its lines, a missing file), `-db` with `-gilist` (a text list converted to
`<out>.p.gil`, a binary one, a missing one, an output directory that does not exist), `-seqidlist`,
`-taxidlist` and `-oid_masks` (existing and missing), the default title `<db> limited by <list>`,
`-title`, `-out`, `-num_volumes`, `-vdblist` and `-vdblist_file` (a new `.nvl`, and the line added to
a `-dblist` alias file), `-gi_file_in` and `-gi_file_out` (the default name `<in>.bgl`, a directory
that does not exist, and a name of 256 bytes that cannot be written: NCBI's `... not written`),
`-seqid_file_in`, `-seqid_file_out`, `-seqid_title` (the default title is the input's name without its
directory), `-seqid_db` (in a directory, and a version 5 database renamed `p\xe9` whose `.pin` and
LMDB file were patched to the new name: `patch`), `-seqid_file_info` (of `ids_db_e9.bsl`, which NCBI
2.17.0 wrote on 2026-11-25 with `-seqid_db` of that renamed database, so its `DB Vols:` line holds the
byte), `-logfile`, and constrained values (`-dbtype`, `-seqid_dbtype`) and numbers (`-num_volumes`,
`-oid_masks`, added and frozen at 04:0x after review F3) that NCBI's reader refuses with the bytes.
With them the date of a binary seqid list (`.bsl`: its date and the date's length byte) is masked
too (`undated`), and its 8-byte size, which counts the date, is zeroed when it holds the length of
the file, so a wrong size still differs (review F1; `undated_seqid_list` of `tests/temp_support`).
