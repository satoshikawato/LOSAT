# Environment variables

LOSAT reads environment variables and NCBI configuration files (`<program>.ini`, `.ncbirc`) the same way as NCBI BLAST+ 2.17.0, with explicit refusals if a setting requests an unimplemented backend or unsafe behaviour.

When LOSAT refuses an environment variable or configuration entry, it exits with status `1` and names the offending variable or file and key.

## Finding databases and configuration files

| Variable or file | Effect |
| :--- | :--- |
| `BLASTDB` | Directories searched for `-db` names after the current working directory (`:`-separated on Unix, `;`-separated on Windows). Also searched for `taxdb.bti` / `taxdb.btd`, `-seqidlist`, `-taxidlist`, and `taxonomy4blast.sqlite3`. |
| `.ncbirc`, `<program>.ini` | NCBI configuration files (`ncbi.ini` replaces `.ncbirc` on Windows; `<program>.ini` takes precedence). Reads `[BLAST] BLASTDB`, `[BLAST] DATA_LOADERS`, `[NCBI] DONT_USE_NCBIRC`, and `[WINDOW_MASKER] WINDOW_MASKER_PATH`. |
| `NCBI_CONFIG__BLAST__BLASTDB` | Sets `[BLAST] BLASTDB` from the environment ahead of config files. |
| `NCBI_CONFIG__BLAST__DATA_LOADERS` | Controls active data loaders. When set cleanly to `none` (or omitting `blastdb`/`genbank`), Seq-id lines are parsed as FASTA text instead of failing. |
| `NCBI_CONFIG_PATH`, `NCBI`, `HOME` | Override configuration file search directories (`NCBI_DONT_USE_LOCAL_CONFIG` skips `.` and `HOME`). |
| `NCBI_DONT_USE_NCBIRC` | Skips reading `.ncbirc` (`NCBI_CONFIG__NCBI__DONT_USE_NCBIRC` has the same effect). |

Run `losat blastdbcmd -show_blastdb_search_path` to inspect the active database search order in your shell environment. See [BLAST databases](databases.md#where-losat-looks-for-a-database) for resolution order.

## WindowMasker paths (`-window_masker_taxid`)

`blastn -window_masker_taxid <taxid>` locates `<directory>/<taxid>/wmasker.obinary` using the first set source in order:

1. `NCBI_CONFIG__WINDOW_MASKER__WINDOW_MASKER_PATH`
2. `WINDOW_MASKER_PATH`
3. `[WINDOW_MASKER] WINDOW_MASKER_PATH` in `.ncbirc`

An empty value resolves to the current working directory; taxid `0` disables WindowMasker. Refused edge cases (such as 17–251-bit ASCII count units where NCBI behaviour is undefined) are noted in [What LOSAT refuses](refusals.md).

## Search batching and statistics modifiers

These standard NCBI environment variables alter search execution or statistics identically to NCBI BLAST+ 2.17.0:

| Variable | Effect |
| :--- | :--- |
| `BATCH_SIZE` | Number of query letters searched per batch (`psiblast` always runs one query per batch). |
| `CHUNK_SIZE`, `OVERLAP_CHUNK_SIZE` | Long-query splitting window size and overlap (`blastx`/`tblastx` require multiples of 3; `rpsblast`/`rpstblastn` never split queries). |
| `ARCHIVE_FORMAT` | Encoding for `-outfmt 11`: unset/empty for ASN.1 text, `bin*` for binary ASN.1, `xml*` for XML. |
| `ADAPTIVE_CBS` | Enables adaptive composition-based statistics hit list expansion/trimming in `blastp`, `psiblast`, `rpsblast`, `blastx`, and `tblastn`. |
| `CTOOLKIT_COMPATIBLE` | Prints `(bits)` instead of `(Bits)` in `-outfmt 0` headers. |
| `OLD_FSC` | Disables finite-size correction in E-value calculations when set (applies to `blastp`, `blastx`, `tblastn`, `psiblast`, `rpsblast`). |
| `PRE_FETCH_SEQS_LIMIT` | Controls `-db` alignment subject sequence preloading threshold per batch (`0` disables preloading). |

## Database building and temporary files

| Variable | Effect |
| :--- | :--- |
| `TZ`, `TZDIR` | Time zone used by `makeblastdb` for database metadata timestamps (`UTC` if unavailable in WASI). |
| `BLASTDB_LMDB_MAP_SIZE` | Initial LMDB map size when creating BLAST v5 databases. |
| `MAX_LMDB_TXN_ENTRY`, `LMDB_MIN_SPLIT_SIZE`, `LMDB_SPLIT_CHUNK_SIZE` | Controls LMDB index write chunking in `makeblastdb`. |
| `TMPDIR` (`TEMP` / `TMP` on Windows) | Fallback directory when native search commands decompress `.gz`/`.bz2`/`.zst` inputs and the input directory is read-only. WASI builds have no temporary directory and require uncompressed inputs or `-query -`. |

## Variables that LOSAT refuses or ignores

- **Refused (exits with status `1`):**
  - `BL2SEQ_LEGACY` (in `-subject` searches)
  - `SEQ_ID_PREFER_ACCESSION_OVER_GI`, `SEQ_ID_AVOID_GI` (with `-db` parsed IDs)
  - `NCBI_MAKEBLASTDB_DONT_SCAN_BIOSEQ_FOR_CFASTAREADER_USER_OBJECTS` (in `makeblastdb`)
  - Non-Boolean `BLAST_USAGE_REPORT` or `NCBI_CONFIG__BLAST__BLAST_USAGE_REPORT` (Boolean values such as `false` are accepted; LOSAT never transmits usage telemetry)
  - NCBI C++ diagnostics/abortion overrides (`DIAG_*`, `ABORT_ON_THROW`, `DEBUG_STACK_TRACE_LEVEL`, `EXCEPTION_STACK_TRACE_LEVEL`, `NCBI_ABORT_ON_*`, `LOG_TRUNCATE`, `LOG_NOCREATE`, `LOG_PERFLOGGING`, `THREAD_STACK_SIZE`) and unrecognized `NCBI_CONFIG_*` keys
- **Ignored (no effect on LOSAT output):**
  - `IGDATA` (specific to `igblastn`/`igblastp`)
  - `BLAST_MT_QUERY_BATCH_SIZE` (LOSAT always emits deterministic single-thread-equivalent output across any `-num_threads`)
  - `BLASTMAT`, `BLASTDB_NUCL_DATA_LOADER`, `BLASTDB_PROT_DATA_LOADER`, `[NCBI] DATA`

Variables prefixed with `LOSAT_` are internal development options and are unsupported for production use; keep all `LOSAT_*` variables unset for standard runs.
