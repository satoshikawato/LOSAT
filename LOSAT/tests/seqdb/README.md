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
    (`unsupported:TEXT`, `format:TEXT`), empty where it must do as NCBI.
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
`empty.bin` are not databases). `LOSAT_TEST_TAX_SQLITE` names the real file for the `#[ignore]` comparison.
