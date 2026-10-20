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
  A column is the files `.?aa` (index), `.?ab` (big-endian data) and `.?ac` (little-endian data).
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

`make_fixtures.py` rebuilds everything but `expected/search_path/`. NCBI BLAST+ is the comparison oracle only; the tests
read the committed files and run no NCBI program. A rebuild changes the dates stored in the
index files and in `expected/`.
