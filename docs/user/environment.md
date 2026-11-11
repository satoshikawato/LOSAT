# Environment variables

LOSAT reads its environment as NCBI BLAST+ 2.17.0 reads it, with the exceptions named on this
page: a variable here has the same effect as in NCBI BLAST+, or LOSAT stops with an error because
it does not do what the variable asks. The `LOSAT_*` variables are for development and are listed
at the end.

A variable or configuration entry that LOSAT refuses stops the run with a message that names it
(for an entry of a configuration file: the file and the entry) and exit status 1. When the checks
run:

- NCBI's toolkit variables, the `NCBI_CONFIG_*` variables and the configuration files: when the
  program starts, before anything is written (`-h`, `-help` and `-version` still answer).
- `SEQ_ID_PREFER_ACCESSION_OVER_GI`, `SEQ_ID_AVOID_GI` and a `[BLAST] BLASTDB` value that LOSAT
  cannot read: where the `-db` database opens, which in every search program is before an empty
  query file ends the run.
- The other variables of a search: after the options and the input files, before the search. An
  empty query file ends the run first, as in NCBI BLAST+ (`Query is Empty!`, exit status 0).
- The query splitting variables in blastx (at every query batch) and tblastx, the splitting values
  that LOSAT refuses in the other programs, and `ADAPTIVE_CBS` in tblastn: at a query batch, so
  the reports of earlier batches, or a report's opening lines, may already be written.

## Finding a database

| Variable or file | Effect |
| :--- | :--- |
| `BLASTDB` | Directories where `-db` names are looked for, separated by `:` (`;` on Windows). The search order is the current directory, then `BLASTDB`, then `[BLAST] BLASTDB` of NCBI's configuration. `taxdb.bti` is found the same way, with `taxdb.btd` in its directory, and so are the files of `-seqidlist`, `-taxidlist` and their negative forms. |
| `.ncbirc`, `<program>.ini` | NCBI's configuration files (for example `blastn.ini`). LOSAT reads `[BLAST] BLASTDB` and `[BLAST] DATA_LOADERS` from them: the first file of each name on the search path, `<program>.ini` before `.ncbirc`. |
| `NCBI_CONFIG__BLAST__BLASTDB` | `[BLAST] BLASTDB` set in the environment; it comes before the files. |
| `[BLAST] DATA_LOADERS`, `NCBI_CONFIG__BLAST__DATA_LOADERS` | Which of NCBI's data loaders are on. With none on (`none`, or neither `blastdb` nor `genbank`, matched anywhere in the value in any case; or an empty `NCBI_CONFIG__BLAST__DATA_LOADERS`), a query or subject line that is a Seq-id is read as FASTA text instead of being refused. In a configuration file an empty value counts as no entry. |
| `NCBI_CONFIG_PATH`, `NCBI`, `HOME`, `NCBI_DONT_USE_LOCAL_CONFIG` | Where the configuration files are looked for: the directories of `NCBI_CONFIG_PATH` if it is set (an empty element stands for the default path); otherwise the current directory, the home directory, `NCBI`, `/etc` (the Windows directory on Windows) and the directory of the `losat` executable. `NCBI_DONT_USE_LOCAL_CONFIG` leaves out the current directory and the home directory. |
| `NCBI_DONT_USE_NCBIRC`, `NCBI_CONFIG__NCBI__DONT_USE_NCBIRC` | Do not read `.ncbirc`. |

Where LOSAT looks differs in three cases. The home directory: LOSAT reads `HOME`, and
`USERPROFILE` on Windows; NCBI BLAST+ reads `APPDATA` first on Windows, and without `HOME` on Unix
it uses the account's home directory. An empty `NCBI_CONFIG_PATH`: NCBI BLAST+ then reads no
configuration file, LOSAT the default path. A program run through a symbolic link (such as a
`blastn` shim): NCBI BLAST+ also looks in the directory of the link, LOSAT only in that of the
`losat` executable.

A configuration file may hold only these entries: `[BLAST]` `BLASTDB` and `DATA_LOADERS`, and
`[NCBI] DONT_USE_NCBIRC` (in `<program>.ini` it skips `.ncbirc`), which LOSAT reads as NCBI
BLAST+ does; `[BLAST] MAX_SEQID_LENGTH`, which `makeblastdb -parse_seqids` refuses (see below);
and `BLASTMAT`, `BLASTDB_NUCL_DATA_LOADER`, `BLASTDB_PROT_DATA_LOADER`, `IGDATA`,
`BLAST_USAGE_REPORT` (a Boolean) and `[NCBI] DATA`, which change no output of LOSAT's programs. A
file with any other entry, one that LOSAT cannot read, or one that NCBI BLAST+ reports as
malformed, stops the run, and so does an `NCBI_CONFIG__<SECTION>__<NAME>` variable for another
entry. A value with a backslash or a quote inside stops every program for `DATA_LOADERS`, a run
that opens a database for `BLASTDB`, and `makeblastdb -parse_seqids` for `MAX_SEQID_LENGTH`.

## Variables that change a search or a report

These change the output as they change NCBI BLAST+'s output. They are not needed for normal use.

| Variable | Effect |
| :--- | :--- |
| `BATCH_SIZE` | The number of query letters searched together (NCBI's query batch). A value that is not an integer stops the run. `psiblast` searches one query at a time and does not read it, as in NCBI BLAST+. |
| `CHUNK_SIZE`, `OVERLAP_CHUNK_SIZE` | The length of the pieces a long query is split into, and their overlap. An empty value counts as unset; a value that is not an integer stops the run, and a non-integer `OVERLAP_CHUNK_SIZE` also stops `rpsblast`, where NCBI BLAST+ does not read it. With a translated query (`blastx`, `tblastx`) a chunk size that is not a multiple of 3 is NCBI's error `Split query chunk size must be divisible by 3` (exit status 3). Splitting values that NCBI BLAST+ does not search correctly are refused: a negative chunk size above a negative overlap, an overlap so close to the chunk size that NCBI splits a chunk again, and in blastx an overlap that is not a multiple of 3. |
| `ARCHIVE_FORMAT` | The encoding of `-outfmt 11`: unset or empty for text ASN.1, a value starting with `bin` (any case) for binary ASN.1, `xml` (any case) for XML; any other value writes no archive, as in NCBI BLAST+. |
| `CTOOLKIT_COMPATIBLE` | Any value: `-outfmt 0` writes `(bits)` in place of `(Bits)` above the hit list. |
| `OLD_FSC` | Any value, an empty one too: `blastp`, `blastx` and `tblastn` compute E-values and cutoffs without the finite-size correction, as NCBI BLAST+ does for any `-searchsp`. It changes nothing in `blastn`, `tblastx` or an ungapped `blastx` search. `psiblast` and `rpsblast` refuse it (see [Variables that LOSAT refuses](#variables-that-losat-refuses)). |
| `PRE_FETCH_SEQS_LIMIT` | With `-db`, an integer is the largest number of alignments of a query batch for which NCBI BLAST+ reads the report's sequences from the database before the report (0: never); LOSAT applies the same rule. A value that is not an integer stops the run (NCBI BLAST+ stops before its report), except when every format is `-outfmt 11`, which `psiblast` and `rpsblast` do not write. |

One variable is not yet read as NCBI BLAST+ reads it in every program; do not set it:

- `ADAPTIVE_CBS`: with composition-based statistics and a `-max_target_seqs` below 1000, NCBI
  BLAST+ keeps a larger preliminary hit list and then trims it by E-value (except in psiblast).
  `blastp` and `psiblast` keep the larger list, but `blastp` does not trim it; `blastx` does not
  read the variable; `tblastn` with `-comp_based_stats 2` refuses it.

## makeblastdb

| Variable | Effect |
| :--- | :--- |
| `TZ`, `TZDIR` | The time zone of the date stored in the database and of the times in the log (Unix; on Windows LOSAT writes them in UTC, where NCBI BLAST+ uses the local time). |
| `BLASTDB_LMDB_MAP_SIZE` | The map size of the LMDB index of a version 5 database. |
| `MAX_LMDB_TXN_ENTRY`, `LMDB_MIN_SPLIT_SIZE`, `LMDB_SPLIT_CHUNK_SIZE` | How NCBI's makeblastdb splits the writes of the LMDB index. |

## Temporary files

`TMPDIR` (`TEMP` or `TMP` on Windows): where a search program inflates a `.gz`, `.bz2` or `.zst`
`-query` or `-subject` when the input's own directory refuses the file or has no room (see
"Common Command-Line Options" in the [README](../../README.md#common-command-line-options)). The
WASI command has no temporary directory: it inflates beside the input or stops with an error.

## Variables that LOSAT refuses

Each of these makes NCBI BLAST+ do something that LOSAT does not do, so LOSAT stops instead of
giving a different output.

| Variable | Where | What it does in NCBI BLAST+ |
| :--- | :--- | :--- |
| `BL2SEQ_LEGACY` | a search with `-subject` | Searches each subject on its own and writes the legacy bl2seq report. |
| `OLD_FSC` | `psiblast`, `rpsblast` | Computes E-values without the finite-size correction. |
| `SEQ_ID_PREFER_ACCESSION_OVER_GI`, `SEQ_ID_AVOID_GI` | the search programs, with a `-db` built with `-parse_seqids` | Changes which id of a sequence is shown. |
| `NCBI_MAKEBLASTDB_DONT_SCAN_BIOSEQ_FOR_CFASTAREADER_USER_OBJECTS` (a non-empty value) | makeblastdb | Stores other titles for some FASTA definition lines (with a tab, a Ctrl-A, or white space at their ends). |
| `NCBI_CONFIG__BLAST__MAX_SEQID_LENGTH`, `[BLAST] MAX_SEQID_LENGTH` (a value other than 0 that NCBI converts) | `makeblastdb -parse_seqids` | Sets the longest Seq-id that makeblastdb parses. |
| `BLAST_USAGE_REPORT` with a value that is not a Boolean | every program | Aborts. (A Boolean is accepted; LOSAT reports no usage.) |
| `DIAG_*` (any), `ABORT_ON_THROW`, `DEBUG_STACK_TRACE_LEVEL`, `EXCEPTION_STACK_TRACE_LEVEL`, `NCBI_ABORT_ON_COBJECT_THROW`, `NCBI_ABORT_ON_NULL`, `LOG_TRUNCATE`, `LOG_NOCREATE`, `LOG_PERFLOGGING`, `THREAD_STACK_SIZE` | every program | Parameters of NCBI's diagnostics, error handling and threads, which can change its messages or exit status (some only on an error). |
| `NCBI_CONFIG_*` other than `NCBI_CONFIG_PATH` and the entries listed above | every program | Set entries of NCBI's configuration. |

An environment prepared for NCBI BLAST+ that sets one of the variables of the last three rows (for
example `DIAG_POST_LEVEL` in a container) stops every LOSAT program before it runs; the other rows
stop only the programs and runs named in them.

Not read: `BLASTMAT`, because NCBI BLAST+ reads a matrix file from it only for a protein matrix
that it does not have built in, and LOSAT has the same built-in matrices and no others; `IGDATA`,
which only NCBI's igblastn and igblastp read; and `BLAST_MT_QUERY_BATCH_SIZE`, the query batches of
NCBI's searches split by queries (`-mt_mode 1`, or NCBI's own choice at `-num_threads` 2 or
more), which LOSAT does not split that way (its output is NCBI's `-num_threads 1`
output).

## WebAssembly

When the host sets `PWD`, the WASI command first changes to that directory, so relative file
names are read from it.
LOSAT Web gives the engine an empty environment: none of these variables applies there, and NCBI's
configuration files are not read.

## Development variables (unsupported)

Every variable whose name starts with `LOSAT_` (`LOSAT_X_*`, `LOSAT_TIMING`, `LOSAT_X_STATS`,
`LOSAT_X_REACH`, `LOSAT_DEBUG_*`, `LOSAT_TRACE_*` and the others), `BLEMIR_DEBUG` and
`BLEMIR_DEBUG_WINDOW` are development variables, as are the cargo features `diagnostics`,
`xstats`, `blastx-worker-probe` and `tblastx-wasm-scalar`. They are unsupported and may change or
disappear without notice; do not set them in production runs.

- The performance switches (the `LOSAT_X_*` list under "Performance switches" in
  [docs/agents/diagnostics.md](../agents/diagnostics.md)) choose between ways of computing the same
  output and do not change it. The defaults (no variable set) are what the release tests cover.
- The debug and trace variables write to standard error, and some of them turn off a faster path
  (any `LOSAT_DEBUG_*` variable turns off TBLASTX's pipelined batches).
- A performance switch is on unless it is `0`; most other development variables, `LOSAT_X_STATS`
  and `LOSAT_TIMING` among them, are on when set to any value, even `0`.
- A diagnostic variable set for a build without the `diagnostics` feature is ignored with a note
  on standard error.

[docs/agents/diagnostics.md](../agents/diagnostics.md) describes each development variable.
