# Fixtures of the database reader (`LOSAT/src/seqdb`)

Small NCBI BLAST databases for the unit tests in `LOSAT/src/seqdb/tests.rs`.

- `prot.faa`, `nucl.fna`: the FASTA files the databases were built from.
- `db/`: databases written by makeblastdb 2.17.0 (versions 4 and 5, single and in several
  volumes, with and without `-parse_seqids`), and alias files written by hand or by
  blastdb_aliastool. Kept: the files LOSAT reads (`.pin .psq .phr`, `.nin .nsq .nhr`, `.pal`,
  `.nal`) and the LMDB file of each version 5 database (`.pdb`, `.ndb`), which must exist.
- `expected/NAME.info.tsv`, `expected/NAME.entries.tsv`: what NCBI's blastdbcmd reads back
  from each database (`-info`; `-entry all -outfmt "%o %l %i %s %t"`).
- `expected/date_pairs.tsv`: the date blastdbcmd reports for a database of two volumes with
  the given date strings (first, second, reported or `ERROR`).
- `expected/search_path/`: what NCBI's blastp and blastdbcmd 2.17.0 do with the name given to
  `-db` (tests `a_name_is_found_where_ncbi_finds_it`,
  `relative_directories_are_used_as_they_are_written`,
  `files_that_cannot_be_opened_give_ncbis_errors` in `LOSAT/src/seqdb/tests.rs` and
  `the_registry_value_of_blastdb_is_ncbis` in `LOSAT/src/blastinput/ncbi_environment.rs`).
  - `layout.tsv`, `open_layout.tsv`: a directory tree made of the databases of `db/` (lines
    `dir`, `copy`, `text`, `symlink`, `patch`, `trunc`, `mode`); the tests build it in a
    temporary directory. `{R}` stands for the root of the tree.
  - `cases.tsv`: one line per run of `blastp -db NAME -query Q -num_threads 1 -outfmt 6`: the
    current directory, `BLASTDB` (`<unset>`: not set), the registry value `[BLAST] BLASTDB`
    (`<none>`: no entry), the name, and NCBI's result. `ok`: title, number of sequences and the
    "Volumes:" list of `blastdbcmd -info` (names separated by ` | `). `error`: the text after
    `BLAST Database error: ` and the exit status. `other`: NCBI failed in another way (its
    message). Column `losat`: what LOSAT's reader does where it refuses the database on purpose
    (`unsupported:TEXT`, `format:TEXT`), `parsed` where it opens a database with parsed Seq-ids
    (which the search programs read), empty where it must do as NCBI.
  - `relative.tsv`: runs from the crate root with relative directories in `BLASTDB` and in the
    registry value (`{P}`: the crate root).
  - `registry.tsv`: the value of `[BLAST] BLASTDB` that NCBI's registry gives for an environment
    variable, a `blastp.ini` and a `.ncbirc` (read off the search path that the "No alias or index
    file found" message prints).
  - `open_cases.tsv`: files without read permission, a directory that cannot be searched, alias
    files that name each other through the search path, and a database of six volumes under
    `ulimit -n` (blastp and blastdbcmd).
  Fields use backslash escapes (`\\`, `\t`, `\n`), `\s` for a space at either end and `\e` for an
  empty field. The tables
  were written by `db3_path_cases.py` and `db3_path_open_cases.py` of the database search
  campaign's task folder (outside the repository), which hold the list of cases.

Mask data (`makeblastdb -mask_data`), read by `LOSAT/src/seqdb/column.rs` and `mask.rs`; the
tests are in `LOSAT/src/seqdb/mask_tests.rs`:

- `mask_nucl.fna`, `mask_prot.faa`: the FASTA files of the masked databases (no
  `-parse_seqids`: the masks name the sequences by their number, `local id N`).
- `mask_wm.asnt`, `mask_wm_opts.asnt`, `mask_wm_colon.asnt`: hand-written windowmasker mask files
  (`Blast-db-mask-info` in ASN.1 text; the second states the option `-t_extend 3`, the third
  `-x a:b`).
- `mask_dust.asnb`, `mask_seg.asnb`: the dust and seg masks of the same databases, written by
  dustmasker and segmasker (`-outfmt maskinfo_asn1_bin`, ASN.1 binary).
- `db/mask_*`: databases written by makeblastdb 2.17.0 with mask data: `mask_nucl_v4`,
  `mask_nucl_v5` (needs its `.ndb`), `mask_nucl_mv` (an alias of two volumes), `mask_nucl_ids`
  (ids 100 and 101 from `-mask_id`, a description with `:`), `mask_nucl_z`, `mask_nucl_y`,
  `mask_nucl_opts` (ids 30 and 31), `mask_prot_v4` (seg, id 21), `mask_remap` (an alias of five
  volumes that give the same algorithms different ids; `mask_nucl_plain` has no column),
  `mask_nucl_colon` (id 31, options `-x a:b`) and `mask_nucl_colon_id` (id 100, name `c:id` from
  `-mask_id c:id`): a `:` in a printed field is stored as 0x01 and restored by NCBI.
  A column is the files `.?aa` (index), `.?ab` (big-endian data) and `.?ac` (little-endian data);
  NCBI reads the `.?ac` when it exists, else the `.?ab` (the ranges of masks then come out with the
  bytes of each number reversed).
- `db/user-column-db.{pin,psq,phr,paa,pab}`: a copy of NCBI's
  `c++/src/algo/blast/unit_tests/seqdb_reader/data/user-column-db` (598d8ae6; written by NCBI's
  `CWriteDB::CreateUserColumn` in 2007): six proteins and the user column `comedy`, whose data
  file is the big-endian `.pab` alone; no masks. Its `.pnd .pni .psd .psi` and `.pog` files are
  left out.
- `db/seqn.{nin,nsq,nhr,nni,nnd}` and `db/seqp.{pin,psq,phr,pni,pnd}`: copies of NCBI's
  `c++/src/algo/blast/unit_tests/seqdb_reader/data/seqn` and `seqp` (598d8ae6; version 4, GIs with
  DDBJ, GenBank and other accessions), whose `GetBioseq` and `GetHdr` texts NCBI's cases GetBioseqN,
  GetBioseqP, GetHdrN and GetHdrP expect (`LOSAT/tests/blastdbcmd_retrieve.rs`, `-outfmt %b` and
  `%d`). Their string-id ISAM files (`.nsd .nsi .psd .psi`), the `.nhd .nhi .phd .phi` and `.nog .pog`
  files are left out: the cases look the sequences up by GI.
- `expected/mask_*.algorithms.txt`: the table that NCBI prints for the algorithms of the
  database (`GetAvailableMaskAlgorithmDescriptions`, with its leading empty line).
- `expected/mask_*.masks.tsv`: rows `algorithm<TAB>oid<TAB>ranges` from
  `blastdbcmd -outfmt "%o\t%m<algorithm>"`; ranges are `start-end;` repeated (end exclusive) or
  `none`. In `mask_remap`, NCBI's rows for an algorithm in a volume that has a mask column but
  not that algorithm come from a stale cache; LOSAT answers `AlgorithmNotInVolume` there.
- `make_mask_fixtures.py` rebuilds these (`--ncbi-bin DIR --work DIR`); a rebuild changes the
  dates stored in the index and column files. The two `mask_nucl_colon*` databases were added
  later by running it into a scratch directory and copying only their files, so that the dates
  of the other fixtures stay.

Mask files (the input of `makeblastdb -mask_data`, read by `LOSAT/src/writedb/mask_info/`; the tests
are in `mask_info/tests.rs`). `mask_parse/` has files that NCBI's tools write and files that
exercise NCBI's reader, with `mask_parse/expected.tsv` recording what `makeblastdb` 2.17.0 does with
each of them:

- `mask_parse/{dust,seg,wm,c2m}.{asnt,asnb}`: dustmasker, segmasker, windowmasker and
  convert2blastmask output (`-outfmt maskinfo_asn1_text` / `maskinfo_asn1_bin`) for the same input;
  `*_ps.*` the same with `-parse_seqids` (Seq-ids of every kind: `local str`, `general`, `other`,
  `genbank`, `pdb`, `patent`, `gi`, ...). Inputs: `mask_nucl.fna`, `mask_prot.faa`,
  `mask_parse/{ids_nucl.fna,ids_prot.faa,lower_nucl.fna,wm_nucl.fna}`.
- `mask_parse/h_*`: hand-written text files (comments, strings, two objects, `more TRUE` with a
  continuation, empty `masks`, different Seq-ids in one Seq-loc, enumerations as numbers,
  every Seq-loc alternative and the combination of several Seq-locs of one sequence with empty and
  whole ranges (`h_seqloc_*`, `h_unsupported_*`), files cut inside an object, garbage after the
  last object, files that are neither text nor binary ASN.1, NCBI's pointer syntax `NULL`, `@n`
  and `:Type` (`h_pointer_*`), ...); `mask_parse/b_*`: damaged and hand-made binary files
  (`b_seqloc_*`: the Seq-loc alternatives in BER; `b_pointer_*`: the pointer syntax in BER;
  `b_length_long_wrap_*`: 8-byte lengths whose end wraps past 2^64).
- `mask_parse/expected.tsv`: per file the exit status and standard error of
  `makeblastdb -in <fasta> -dbtype <type> -mask_data <file>` (no `-parse_seqids`), and the masks that
  `blastdbcmd -outfmt '%o\t%m<algorithm>'` reads back from the database it wrote (`RUN`, `NODB`
  and `MASK` rows; the last field of `RUN` is the algorithm's options as `blastdbcmd -info` prints them;
  the header of the file explains them).
- `make_mask_parse_fixtures.py` rebuilds these (`--ncbi-bin DIR --work DIR`); the tools and the
  inputs are deterministic, so a rebuild changes nothing.
- `stl_sort/cases.txt`: orders in which C++ `std::sort` leaves ranges under NCBI's
  `CRange_Less` / `CRange_ReverseLess` (`stl_sort/README.md`), for the test of the sort in the port
  of `CSeq_loc::Add`.

`make_fixtures.py` rebuilds everything but `expected/search_path/`. NCBI BLAST+ is the comparison oracle only; the tests
read the committed files and run no NCBI program. A rebuild changes the dates stored in the
index files and in `expected/`.

`lmdb/` and `taxsqlite/` hold the fixtures of the read-only LMDB, taxid-lookup and taxonomy readers
(`src/seqdb/lmdb_read.rs`, `tax_files.rs`, `tax_sqlite.rs`). `lmdb/*.{pdb,ndb,ptf,ntf,pto,nto,pot,not,pos,nos}` are the
lookup files of small databases written by makeblastdb 2.17.0 (`p5`, `n5`, `mv` of four volumes, `dup` with 700
OIDs under one accession, `thr500`/`thr501` at the 500-duplicate sub-page threshold, `mvbig` with 400 volumes,
`tx` with one to three taxids per OID, taxid 0 and a taxid listed twice). `lmdb/expected/` holds what a Python
reference reader (`dump.py`) printed for them and the OID lists NCBI's `blastdbcmd -taxids` and
`blastp -negative_taxids` selected for `tx` (`tx_oracle.txt`) and for the NCBI-distributed `18S_fungal_sequences`
database (`18S_oracle.txt`; the database itself is not in the repository, the `#[ignore]` tests read it from
`LOSAT_TEST_18S_DB`, with the dump files in `LOSAT_TEST_18S_EXPECTED`). `taxsqlite/*.sqlite3` are small databases made
with Python's `sqlite3` from the schema text of NCBI's `taxonomy4blast.sqlite3` (`p512`: 3005 rows in 512-byte pages,
`utf16`, `noindex`, `notable`, and `cycle` with a two-node parent cycle); their `.ref` files list, for every node,
the count and FNV-1a digest of the descendants that NCBI's recursive query gave when run by SQLite (`notdb.bin` and
`empty.bin` are not databases). `LOSAT_TEST_TAX_SQLITE` names the real file for the `#[ignore]` comparison, `LOSAT_TEST_TAX_SQLITE_EXPECTED` the file that
`tax_sqlite_ref.py` writes from it (Python's `sqlite3` running NCBI's query; the test itself starts no program).
`lmdb/tx.{pin,psq,phr}` complete `tx` as a database (made from `tx.fa` with the taxid map `tx.map`; `tx_q.fa` holds
one query per sequence): `src/seqdb/filter.rs` restricts it to each taxid set of `tx_oracle.txt` and compares the
OIDs with NCBI's.

## Mask-column writer fixtures (`mask_writer/`)

`make_mask_writer_fixtures.py` (`--ncbi-bin DIR --work DIR`) builds the inputs and NCBI's outputs
for the tests of the mask-column writer (`src/writedb/mask_data/tests_e2e.rs`, which reads the mask files with `process_mask_data`; the same file runs every file of `mask_parse/` end to end). Each row of
`mask_writer/cases.tsv` is one makeblastdb 2.17.0 run (`-in` files, `-mask_data` files, `-mask_id`,
`-mask_desc`, `-max_file_sz`); its database files are in `mask_writer/db/`, its exit status and
standard error in `mask_writer/expected/<case>.rc` and `.stderr`, its log (standard output, with the
build time, database path and seconds masked as `<T>`, `<OUT>`, `<S>`) in `.stdout`, and what blastdbcmd reads back
in `expected/<case>.algorithms.txt` and `.masks.tsv` (as for the reader fixtures above). The
mask files are ASN.1 text: hand-written ones (`*.asnt`) and `dust.asnt` / `seg.asnt` from
dustmasker / segmasker. The cases cover the combination of several Seq-locs of one sequence, the
ids of algorithms, masks that match no sequence, ranges past the end of a sequence, and volumes
without a masked sequence. A rebuild changes the dates stored in the index and column files.
`src/writedb/mask_data/tests_cli.rs` runs the cases through the command line (`makeblastdb::run_with`) and
compares the exit status, the log and the standard error of each with NCBI's, and every file of the database.
`mask_parse/expected_stdout.tsv` holds the log of each run of `mask_parse/expected.tsv` the same way.

## GI-indexed masks (`gi_mask/`)

`make_gi_mask_fixtures.py` (`--ncbi-bin DIR --work DIR`, 2026-11-20, makeblastdb, dustmasker, segmasker and
convert2blastmask 2.17.0 with `BLAST_USAGE_REPORT=false`) writes the inputs (`gi_mask/in/`: FASTA files from fixed
seeds, their masks, one hand-written mask file) and runs makeblastdb on each row of `gi_mask/cases.tsv` (name, what
LOSAT does, setup, the words after `makeblastdb`) in a directory of its own, as `losat makeblastdb -gi_mask` is
tested: `-gi_mask` and `-gi_mask_name` on single-volume (version 5 and 4) and multi-volume databases, nucleotide and
protein, the argument rules, masks of sequences without a GI, ranges past a sequence's end, `-logfile`, a missing
directory, a rebuild in the same directory, duplicate GIs of version 4, 1 to 1100 GIs (pages of 512), data volumes
cut by `-max_file_sz`, descriptions that pad the index header, two algorithms of one name (one that writes nothing,
data volumes written one after the other, and the refused cases that write one data file at the same time), and
`-input_type blastdb` from `gi_mask/in/src4.*` (a version 4 database with a mask column that makeblastdb wrote; the
"bytes-total" of each metadata file a run writes is checked to be the sum of the sizes of the files it lists). Each
run's exit status, standard output and error
(the run directory, time and seconds masked as `<C>`, `<T>`, `<S>`) and its files with their sizes are in
`gi_mask/expected/<case>.{rc,stdout,stderr,files}`; `gi_mask/db/<case>/` keeps the files that `-gi_mask` writes or
changes (GI-mask files, alias, metadata, mask column, log) and, for six cases, the volume files too.
These are NCBI's bytes but in one place, the script's `single_volume_exception` (PD-LOSAT-NCBI-DEFECTS approved
exception 25, the Owner's decision of 2026-11-20): for a database of one volume NCBI's alias file names the database
itself by its absolute path (`DBLIST <C>/db`), which no program can open, and makeblastdb ends with `DB alias files
are mutually recursive` (exit status 2, no metadata file). LOSAT names the volume by its file name (`DBLIST db`). So
after each such run the script changes that line in the run directory (before the next case of a chain runs) and,
where NCBI ended on that error, takes the run as the successful run it then is: the line leaves the standard error
(or the `-logfile` file), the exit status is 0, and the metadata file is NCBI's `blastdbcmd -metadata` of the
database. `gi_mask/in/ncbi_single.nal` is NCBI's own alias file of that kind (`{C}` for its directory), which the
cases `x10_*` put beside a database (setup `F=>G`) as an earlier run of NCBI's leaves it.
`tests/makeblastdb_gi_mask.rs` runs every case through the `losat` binary and compares them;
`src/writedb/gi_mask_tests.rs` writes every GI mask of `db/` again from the records, GIs, description and dates it
holds and compares the six files byte for byte (but a name that the alias file lists twice, case `x12`). A rerun changes the dates of the `.gmi`/`.gni`, mask column and
alias files.

## Readers of GI-indexed masks (`gi_mask_read/`)

`make_gi_mask_read_fixtures.py` (`--ncbi-bin DIR --work DIR`, 2026-11-20, makeblastdb, blastdbcmd,
blastdb_aliastool and the search programs 2.17.0 with `BLAST_USAGE_REPORT=false`) builds the databases of
`gi_mask_read/dbs.tsv` with NCBI's makeblastdb from the writer's inputs (`gi_mask/in/`, the same functions):
multi-volume version 5 and 4, protein, a database without masks, a description of another program, the GI
4294967295 beside a sequence without a GI, duplicate GIs (version 4), an absolute `-gi_mask_name`, data volumes of
10 and 20 bytes (one of them never written, survey D8), a `-gi_mask_name` with a space, and single-volume databases
(whose alias file it rewrites to `DBLIST db`, the form LOSAT's makeblastdb writes: PD-LOSAT-NCBI-DEFECTS approved
exception 25). It then freezes what NCBI's programs print for each row of `gi_mask_read/cases.tsv` (name, databases,
directory under the root, environment, words): `blastdbcmd -info`, `-outfmt %m`, `-mask_sequence_with`, `blastn
-db_soft_mask` / `-db_hard_mask` by id and name, from another directory, with `BLASTDB`, by absolute name, through an
alias file over the database and as `-db 'db db'`, a plain search of each program family, and (rows added
2026-11-23) `blastp`, `blastx`, `tblastn` and `tblastx` with `-db_soft_mask` / `-db_hard_mask` by id (`-outfmt 6`)
and by name (`-outfmt 0`) on the protein databases `m09`, `svp` and the nucleotide databases `m01`, `sv5`; the rows
`raw_*` copy `m01`, `sv5` and `svp` with `makeblastdb -input_type blastdb` and read the copies back (`-info`, `%m`,
a hard-masked search), from another directory, and with `-in 'db db'` (the rows run in order in one root). `gi_mask_read/db/<db>/`
holds the databases (the root `<R>` in their text files), `gi_mask_read/expected/<case>.{rc,stdout,stderr}` the runs
(the root as `<R>`), `gi_mask_read/in/` the queries. `bdc_m09_msw_0` (a protein in FASTA with masks) is approved
difference 5 of PD-LOSAT-BLASTDBCMD-DEVIATIONS: NCBI's headers with the letters of NCBI's `-outfmt %s` run, NCBI's
own FASTA in `bdc_m09_msw_0.ncbi_stdout`. `tests/gi_mask_read.rs` runs every case through the `losat` binary
(searches also with `-num_threads 4`, whose reference is NCBI's `-num_threads 1` output, Owner decision DB-3) and on
databases that `losat makeblastdb` builds;
`src/seqdb/gi_mask_tests.rs` reads every GI of `gi_mask/db/big*` (1 to 1,100 GIs, one and four data volumes).

## Masks ordered by address, and the errors that come first (`mask_order/`)

`make_mask_order_fixtures.py` (`--ncbi-bin DIR --work DIR`, 2026-11-14, makeblastdb 2.17.0 with
`BLAST_USAGE_REPORT=false`) writes the inputs and NCBI's outputs of five runs on `in.fna`, whose first
sequence has the Seq-ids `gi|5` and `gb|AY000005.1|` (one packed number): masks of both on the plus
strand, which makeblastdb writes in the order of the memory addresses of their handles, and an error of
makeblastdb that does not depend on that order (a range of the same sequence past its end, also of a
second algorithm; a range of the next sequence past its end; a second input file without sequences),
plus `tie`, the masks alone. Each row of `mask_order/cases.tsv` (name, `-in` files, `-mask_data` files)
ran with `-dbtype nucl -parse_seqids -title in.fna -out db -blastdb_version 5` and `4`; its exit status,
standard error, log (masked as above) and the names of the files `db.*` it left are in
`mask_order/expected/<case>.v<N>.{rc,stderr,stdout,files}`. `src/writedb/mask_data/tests_cli.rs` runs
them through the command line: LOSAT gives makeblastdb's outcome for the first four and refuses `tie`
(which NCBI builds) at the close, leaving the files of the failed build of `oob_later_sequence`.
