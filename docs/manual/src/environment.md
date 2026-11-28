# Environment variables

LOSAT reads its environment and NCBI's configuration files (`.ncbirc`, `<program>.ini`) as NCBI BLAST+ 2.17.0 reads them, with the exceptions named on this page. A variable listed here has the same effect as in NCBI BLAST+, or LOSAT stops with an error because it does not do what the variable asks. Most runs need none of them, except `BLASTDB` to find databases.

## When LOSAT refuses a setting

A variable or configuration entry that LOSAT refuses stops the run with a message that names it (for an entry of a configuration file: the file and the entry) and exit status 1. When the checks run:

- NCBI's toolkit variables, the `NCBI_CONFIG_*` variables and the configuration files: when the program starts, before anything is written. `-h`, `-help` and `-version` still answer, unless NCBI BLAST+ would write a message about a configuration file before them that LOSAT does not write: then LOSAT stops there.
- `SEQ_ID_PREFER_ACCESSION_OVER_GI`, `SEQ_ID_AVOID_GI` and a `[BLAST] BLASTDB` value that LOSAT cannot read: where the `-db` database opens, which in every search program is before an empty query file ends the run. A `[BLAST] BLASTDB` value that LOSAT cannot read also stops a `-subject` search whose report looks `taxdb.bti` up (a tabular report with `sscinames`, `scomnames`, `sblastnames` or `sskingdoms`, or the organism report, outfmt 18), at the same point. The singular names (`ssciname`, ...) look nothing up for `-subject` records, which have no taxid, so they do not read it, as in NCBI BLAST+.
- The other variables of a search: after the options and the input files, before the search. An empty query file ends the run first, as in NCBI BLAST+ (`Query is Empty!`, exit status 0).
- The query splitting variables in `blastx` (at every query batch) and `tblastx`, and the splitting values that LOSAT refuses in the other programs: at a query batch, so the reports of earlier batches, or a report's opening lines, may already be written.

## Finding a database

| Variable or file | Effect |
| :--- | :--- |
| `BLASTDB` | Directories where `-db` names are looked for, separated by `:` (`;` on Windows). The search order is the current directory, then `BLASTDB`, then `[BLAST] BLASTDB` of NCBI's configuration. `taxdb.bti` is found the same way, with `taxdb.btd` in its directory, and so are `taxonomy4blast.sqlite3` and the files of `-gilist`, `-ipglist`, `-seqidlist`, `-taxidlist` and their negative forms. |
| `.ncbirc`, `<program>.ini` | NCBI's configuration files (for example `blastn.ini`; on Windows `ncbi.ini` stands for `.ncbirc`). LOSAT reads `[BLAST] BLASTDB` and `[BLAST] DATA_LOADERS` from them: the first file of each name on the search path, `<program>.ini` before `.ncbirc`. A file found that cannot be opened counts as no file, and no later directory is searched for it. |
| `NCBI_CONFIG__BLAST__BLASTDB` | `[BLAST] BLASTDB` set in the environment; it comes before the files. |
| `[BLAST] DATA_LOADERS`, `NCBI_CONFIG__BLAST__DATA_LOADERS` | Which of NCBI's data loaders are on. With none on (`none`, or neither `blastdb` nor `genbank`, matched anywhere in the value in any case; or an empty `NCBI_CONFIG__BLAST__DATA_LOADERS`), a query or subject line that is a Seq-id is read as FASTA text instead of being refused. In a configuration file an empty value counts as no entry. |
| `NCBI_CONFIG_PATH`, `NCBI`, `HOME`, `NCBI_DONT_USE_LOCAL_CONFIG` | Where the configuration files are looked for: the directories of `NCBI_CONFIG_PATH` if it is set (an empty element stands for the default path, and an empty `NCBI_CONFIG_PATH` means no directory); otherwise the current directory, the home directory, `NCBI`, `/etc` (the Windows directory on Windows) and the program's directories. `NCBI_DONT_USE_LOCAL_CONFIG` leaves out the current directory and the home directory. |
| `NCBI_DONT_USE_NCBIRC`, `NCBI_CONFIG__NCBI__DONT_USE_NCBIRC` | Do not read `.ncbirc`. |

`losat blastdbcmd -show_blastdb_search_path` prints the directories LOSAT searches for a database name. [BLAST databases](databases.md#where-losat-looks-for-a-database) shows the order with an example.

The home directory is `HOME` (an empty `HOME` means no home directory), or without `HOME` the account's home directory; on Windows it is `APPDATA`, else `USERPROFILE`. Without `HOME` and without an account entry for the user, NCBI BLAST+ looks the account up by the login name, which LOSAT cannot do, so LOSAT stops instead.

The program's directories are, as NCBI BLAST+ finds them, the directory of the `losat` executable (links resolved) on Linux and none elsewhere. With `BLAST_USAGE_REPORT` set to a false value (`false`, `0`, `no`, `off`, ...), or with `NCBI_DONT_USE_NCBIRC` or `NCBI_CONFIG__NCBI__DONT_USE_NCBIRC` set, NCBI BLAST+ makes the search path later, and the directories are then the directory of the name the program was run by (looked up in the current directory or on `PATH` when it has no directory, so a `blastn` link's directory), then that of the `losat` executable when it is another one (on Windows the directory of the `losat` executable).

## What a configuration file may hold

A configuration file may hold only these entries:

- `[BLAST] BLASTDB` and `DATA_LOADERS`, `[NCBI] DONT_USE_NCBIRC` (in `<program>.ini` it skips `.ncbirc`) and `[WINDOW_MASKER] WINDOW_MASKER_PATH` (see below), which LOSAT reads as NCBI BLAST+ does;
- `[BLAST] MAX_SEQID_LENGTH`, which `makeblastdb -parse_seqids` refuses (see [Variables that LOSAT refuses](#variables-that-losat-refuses));
- `[BLAST] BLASTMAT`, `BLASTDB_NUCL_DATA_LOADER`, `BLASTDB_PROT_DATA_LOADER`, `IGDATA` and `BLAST_USAGE_REPORT` (see the refused variables), and `[NCBI] DATA`, which change no output of LOSAT's programs.

A file with any other entry stops the run, and so does an `NCBI_CONFIG__<SECTION>__<NAME>` variable for another entry.

LOSAT reads the files as NCBI BLAST+ 2.17.0 does: comments, sections, entries before the first section (left out), a value continued on the next line after a backslash, quotes at the ends of a value, and backslash escapes. `BLASTDB` and `WINDOW_MASKER_PATH` are read as their bytes, as NCBI BLAST+ reads them; a non-UTF-8 value for `MAX_SEQID_LENGTH` stops `makeblastdb -parse_seqids`.

### Malformed files and their messages

A file that NCBI BLAST+ reports as malformed gets NCBI BLAST+'s messages, where NCBI BLAST+ writes them.

- NCBI BLAST+ reads `.ncbirc` twice, each time with its `Critical:` message: when the program starts (before `-version` answers; not with `BLAST_USAGE_REPORT` set to a false value or `.ncbirc` turned off) and when it loads its configuration (after `<program>.ini`, before the options are read, `-h` included). The program then goes on with the entries above the line it reports.
- A malformed `<program>.ini` ends the program with exit status 2 and NCBI BLAST+'s messages, before the options are read: its report, then every message so far again.
- The messages go to the standard error stream, and with `-logfile` as in NCBI BLAST+: those of the start to the standard error stream and into the file, the others into the file (`-logfile /dev/null` and `-logfile ""` discard them), except that `blastn`, `blastp`, `blastx`, `tblastn` and `psiblast` write the messages again to the standard error stream. The messages of the run go into the file in `tblastx`, `rpsblast` and `rpstblastn`, and to the standard error stream in the other search programs, as in NCBI BLAST+.
- `-conffile FILE` loads FILE in place of `<program>.ini` (a name without a directory is looked for where `<program>.ini` is looked for; a name with one is that file). A FILE that is not a file ends the program with NCBI BLAST+'s message and exit status 255, before `-h`.
- With `-dryrun` the configuration is loaded and the settings that LOSAT refuses are refused as without it (some of them, such as `DIAG_POST_LEVEL`, would make NCBI BLAST+ write a message); the program then stops after its options.
- `rpsbproc` writes its messages to the `-logfile` file as NCBI's `rpsbproc` does, and its program information, progress and dry-run line to the standard error stream; its `-dryrun` runs its own dry run (its start, then the line of its status). `rpsbproc` reads `.ncbirc` once, when its configuration is loaded (it has no usage report): a `.ncbirc` with a syntax error is reported there as NCBI's `rpsbproc` reports it, naming a source file of that `rpsbproc`'s own build (the Linux x64 binary of NCBI's package of 2025-05-13), and the run goes on.
- LOSAT stops instead where it does not write NCBI BLAST+'s message: a file in UTF-16, or a malformed line that is not UTF-8.

The library runs the same start at the beginning of every search call and every tool call, as the command line does in every run: it reads `.ncbirc`, `<program>.ini` and the `NCBI_*` variables and refuses what the command line refuses. Its messages go where the command line writes them, and a start that stops the command line's program is the call's error.

## WindowMasker files of `blastn -window_masker_taxid`

`blastn -window_masker_taxid <taxid>` masks the queries with the counts file `<directory>/<taxid>/wmasker.obinary`. LOSAT finds the directory as NCBI BLAST+ does, in this order: the first one that is set (an empty value too) wins, and an empty value or none means the current directory. With the taxid 0, no WindowMasker is used and no file is needed.

| Variable or entry | Effect |
| :--- | :--- |
| `NCBI_CONFIG__WINDOW_MASKER__WINDOW_MASKER_PATH` | The directory, before every other source. |
| `WINDOW_MASKER_PATH` | The directory, after the variable above. |
| `[WINDOW_MASKER] WINDOW_MASKER_PATH` in `.ncbirc` | The directory, after both variables: the first `.ncbirc` on the search path, also when `<program>.ini` sets `[NCBI] DONT_USE_NCBIRC`; `NCBI_DONT_USE_NCBIRC` and `NCBI_CONFIG__NCBI__DONT_USE_NCBIRC` skip it. `<program>.ini`'s entry is not read. |

The file may be in any of the four counts formats of `windowmasker` (optimized binary, binary, optimized ascii, ascii), as in NCBI BLAST+. A taxid without that file, or with a file that NCBI BLAST+ cannot read, stops the search with NCBI's error at the first query batch the file cannot mask (exit status 3). NCBI BLAST+ reads `.ncbirc` again for the directory each time it masks a query batch, and each chunk of a batch that it splits; a malformed `.ncbirc` gives its message each time, with `[blastn]` after the error code, and the directory comes from the entries above the line it reports. LOSAT does the same.

LOSAT refuses, before any output:

- an ascii counts file whose unit size is 17 to 251 (NCBI BLAST+ masks with units that do not fit its 32-bit units, an undefined result);
- an optimized ascii counts file whose hash parameters line (its third) leaves unread a value that NCBI BLAST+ then checks (the line ends, or a word before it is not a number: NCBI checks it uninitialized, an undefined result);
- a `.ncbirc` that sets `[NCBI] .Inherits`;
- a current directory it cannot read.

The Web GUI refuses `-window_masker_taxid` other than 0.

## Variables that change a search or a report

These change the output as they change NCBI BLAST+'s output. They are not needed for normal use.

NCBI BLAST+ reads `BATCH_SIZE`, `CHUNK_SIZE`, `OVERLAP_CHUNK_SIZE` and `PRE_FETCH_SEQS_LIMIT` as whole numbers: an optional `+` or `-` and digits, within the range of a 32-bit integer, without spaces. Where it reads a value that is not such a number (`abc`, `1.5`, ` 5`, `5 `, `2147483648`), the program stops as NCBI BLAST+ stops: exit status 255 and an error message that names the value, and in a search with `-outfmt 11` NCBI's error archive. The table says where each one stops (in `blast_formatter`, after reading an archive object's results).

| Variable | Effect |
| :--- | :--- |
| `BATCH_SIZE` | The number of query letters searched together (NCBI's query batch). A value that is not a number (an empty one too) stops the run before anything is written; `blastn` stops after the start of the `-outfmt 0` report. `psiblast` searches one query at a time and does not read it, as in NCBI BLAST+. |
| `CHUNK_SIZE`, `OVERLAP_CHUNK_SIZE` | The length of the pieces a long query is split into, and their overlap. An empty value counts as unset; a value that is not a number stops the run at the first query batch, after the start of the `-outfmt 0` report and the warnings about the batch's query records (`blastn` reads `CHUNK_SIZE` before it reads the first batch). `rpsblast` and `rpstblastn` never split a query and, as NCBI BLAST+ does, do not read `OVERLAP_CHUNK_SIZE` (a non-integer value is ignored there; `CHUNK_SIZE` is still read); neither does `blastx` with several queries in a batch. With a translated query (`blastx`, `tblastx`; not `rpstblastn`) a chunk size that is not a multiple of 3 is NCBI's error `Split query chunk size must be divisible by 3` (exit status 3). Splitting values that NCBI BLAST+ does not search correctly are refused: a negative chunk size above a negative overlap, an overlap so close to the chunk size that NCBI splits a chunk again (except in `blastn`, see [Differences from NCBI BLAST+](differences.md#query-chunks-that-ncbi-would-split-again)), and in `blastx` an overlap that is not a multiple of 3 (refused after the results of the query batches before, where the batch has subjects to search whose sequences can be read; otherwise the batch is reported, or fails, as in NCBI BLAST+). Where NCBI BLAST+ stops at such a chunk (one that it would split again in `blastp`, `blastx`, `psiblast` and `tblastn`, or one left without a query, also by a negative pair whose wrapped range covers or passes the batch's end, and in `blastx` by every negative pair), `blastn`, `blastp`, `blastx`, `psiblast` and `tblastn` write the results of the query batches before (and for the first batch, the prologs that come with them) and then stop with the error, before the errors of the database or of `-searchsp`/average subject length that NCBI BLAST+ would meet later (also when the database is left without a subject: a seqid list that keeps none, an alias over volumes without sequences); a batch whose queries are all unsearchable (all `N` or all `X`, for example) is not split. So does a negative pair whose pieces leave gaps between them, where a hit could be lost (NCBI BLAST+ searches with the hits there lost), at the first piece with searchable query letters (in `blastn` after that piece's setup, where the database's errors come first); with no subject to search, or no such piece, the batch is reported as NCBI BLAST+ reports it. |
| `ARCHIVE_FORMAT` | The encoding of `-outfmt 11`: unset or empty for text ASN.1, a value starting with `bin` (any case) for binary ASN.1, `xml` (any case) for XML; any other value writes no archive, as in NCBI BLAST+. |
| `ADAPTIVE_CBS` | Any value, an empty one too, with composition-based statistics: `blastp`, `psiblast`, `rpsblast`, `blastx` and `tblastn` keep a larger preliminary hit list (1,500 subjects, or twice the hit list size and 50 more from a hit list size of 1,000), and with a hit list size (`-max_target_seqs`) below 1,000 all of them but `psiblast` cut that list by E-value before the final alignments, as NCBI BLAST+ does. The reported hits can grow or shrink. |
| `CTOOLKIT_COMPATIBLE` | Any value: `-outfmt 0` writes `(bits)` in place of `(Bits)` above the hit list. |
| `OLD_FSC` | Any value, an empty one too: `blastp`, `blastx`, `tblastn`, `psiblast` (every round) and `rpsblast` compute E-values and cutoffs without the finite-size correction, as NCBI BLAST+ does for any `-searchsp`. It changes nothing in `blastn`, `tblastx` or an ungapped `blastx` search. |
| `PRE_FETCH_SEQS_LIMIT` | With `-db`, an integer is the largest number of alignments of a query batch for which NCBI BLAST+ reads the report's sequences from the database before the report (0: never); LOSAT applies the same rule. A value that is not a number (an empty one too) stops the run after the search of the first query batch, before its report (after the start of the `-outfmt 0` report), except when every format is `-outfmt 11`, which every program but `psiblast` writes without reading the variable. `blast_formatter` also reads `PRE_FETCH_SEQS_LIMIT` after reading an archive object's results, before formatting it: a value that is not a number stops with exit status 255 and NCBI's message, and `2147483647` is refused by LOSAT for `-outfmt` `0` to `4`, `8`, `9` and `12` (see [blast_formatter](programs/blast_formatter.md#what-losat-refuses)). |

## makeblastdb

| Variable | Effect |
| :--- | :--- |
| `TZ`, `TZDIR` | The time zone of the date stored in the database and of the times in the log, with `/etc/localtime` when `TZ` is unset (Unix, and the WASI command from what its host gives it, see [WebAssembly](#webassembly)). On Windows, as Microsoft's C runtime reads it for NCBI BLAST+: `TZ` in its form `tzn[+\|-]hh[:mm[:ss]][dzn]` (hours west of UTC; a daylight-time name means the United States rules), else the zone of the Windows settings with its daylight-time rules; `TZDIR` is not read. |
| `BLASTDB_LMDB_MAP_SIZE` | The map size of the LMDB index of a version 5 database. |
| `MAX_LMDB_TXN_ENTRY`, `LMDB_MIN_SPLIT_SIZE`, `LMDB_SPLIT_CHUNK_SIZE` | How NCBI's `makeblastdb` splits the writes of the LMDB index. |
| `NCBI_MAKEBLASTDB_SKIP_UNVERIFIED_BIOSEQS` | Any non-empty value with `-input_type asn1_txt` or `asn1_bin`: skips `Bioseq` records marked `Unverified` (or inside an `Unverified` `nuc-prot` set), as in NCBI's `makeblastdb`. |
| `NCBI_MAKEBLASTDB_DONT_SCAN_BIOSEQ_FOR_CFASTAREADER_USER_OBJECTS` | Any non-empty value with `-input_type asn1_txt` or `asn1_bin`: ignores `CFastaReader` `DefLine` user objects on ASN.1 `Bioseq` records and builds definition lines from the `Bioseq` descriptors instead. (Refused with `-input_type fasta`; see below.) |

## Temporary files

`TMPDIR` (`TEMP` or `TMP` on Windows): where a search program inflates a `.gz`, `.bz2` or `.zst` `-query` or `-subject` when the input's own directory refuses the file or has no room (see [Differences from NCBI BLAST+](differences.md#compressed--query--subject-and--in)). The WASI command has no temporary directory: it inflates beside the input or stops with an error.

## Variables that LOSAT refuses

Each of these makes NCBI BLAST+ do something that LOSAT does not do, so LOSAT stops instead of giving a different output.

| Variable | Where | What it does in NCBI BLAST+ |
| :--- | :--- | :--- |
| `BL2SEQ_LEGACY` | a search with `-subject` | Searches each subject on its own and writes the legacy bl2seq report. |
| `SEQ_ID_PREFER_ACCESSION_OVER_GI`, `SEQ_ID_AVOID_GI` | the search programs, with a `-db` built with `-parse_seqids` | Changes which id of a sequence is shown. |
| `NCBI_MAKEBLASTDB_DONT_SCAN_BIOSEQ_FOR_CFASTAREADER_USER_OBJECTS` (a non-empty value) | `makeblastdb -input_type fasta` | Stores other titles for some FASTA definition lines (with a tab, a Ctrl-A, or white space at their ends). With `-input_type asn1_txt` or `asn1_bin`, it is read as in NCBI BLAST+ (see above). |
| `NCBI_CONFIG__BLAST__MAX_SEQID_LENGTH`, `[BLAST] MAX_SEQID_LENGTH` (a value other than 0 that NCBI converts) | `makeblastdb -parse_seqids` | Sets the longest Seq-id that `makeblastdb` parses. |
| `BLAST_USAGE_REPORT` with a value that is not a Boolean | every program, also for `-version` and `-h` | Aborts. (A Boolean is accepted; LOSAT reports no usage.) |
| `NCBI_CONFIG__BLAST__BLAST_USAGE_REPORT`, or else `[BLAST] BLAST_USAGE_REPORT` of `.ncbirc`, with a value that is not a Boolean, unless `BLAST_USAGE_REPORT` is false (for `.ncbirc` also unless `.ncbirc` is turned off in the environment) | every program, also for `-version` and `-h` | Aborts. |
| `OBJECTS_PACK_TEXTID`, `OBJECTS_PACK_GENERAL` (any value) | every program | Changes how Seq-ids are matched, and with it which sequences a mask of `makeblastdb -mask_data` reaches. |
| `OBJECTS_SEQ_DESCR_ALLOW_EMPTY` (any value, with an archive as input); `SERIAL_VERIFY_DATA_READ` or `SERIAL_VERIFY_DATA_GET` set to `NO`, `NEVER`, `DEFVALUE` or `DEFVALUE_ALWAYS` (any case) | `rpsbproc` | Changes how NCBI's object stream reads its input. |
| `SERIAL_VERIFY_DATA_READ` set to `NO`, `NEVER`, `DEFVALUE` or `DEFVALUE_ALWAYS` (any case) | a run that reads an input without a member that NCBI's reader requires: an `-in_pssm` file of `psiblast` or `tblastn`, a profile or the binary definition lines (`ASN1_BlastDefLine`) of `makeprofiledb`, a `-mask_data` file or an `-input_type asn1_txt` / `asn1_bin` input of `makeblastdb`, or the definition lines of a database read by a search program, `blastdbcmd`, `blastdbcheck`, `makeblastdb -input_type blastdb`, `dustmasker`, `segmasker`, `windowmasker` or `makembindex` | Reports the member as missing and reads on, where it otherwise stops with an error (`blastdbcheck`: fails the check of that OID). LOSAT stops where NCBI's error would, after the same output. Other runs with these values are answered as NCBI does. |
| `SERIAL_VERIFY_DATA_GET` set to `NO`, `NEVER`, `DEFVALUE` or `DEFVALUE_ALWAYS` (any case) | `psiblast` and `tblastn` with an `-in_pssm` file whose query has a feature table column without the `field-name` that NCBI reads, or `makeblastdb -input_type asn1_txt` / `asn1_bin` with a `Bioseq` whose `Seq-inst.length` (or `Seq-literal.length`) is unset | Takes a default value and goes on, where it otherwise stops with an error. LOSAT stops where NCBI's error would, after the same output. Other runs with these values are answered as NCBI does. |
| `SERIAL_VERIFY_DATA_WRITE` set to `NO`, `NEVER`, `DEFVALUE` or `DEFVALUE_ALWAYS` (any case) | `psiblast` and `tblastn` with an `-in_pssm` file that holds a list element stated as `NULL`, in a run that writes the PSSM: `-export_search_strategy` with `-db`, `-out_pssm`, `psiblast -outfmt 11` with `-db` | Leaves the element out of the PSSM it writes, where it otherwise stops with an error. Other runs with these values are answered as NCBI does. |
| `SERIAL_VERIFY_DATA_WRITE` set to `ALWAYS` or `DEFVALUE_ALWAYS` (any case) | the search programs, with `-outfmt 5` | Stops with an error in the XML report (`ALWAYS`), or adds an empty `BlastOutput_iterations` element to it (`DEFVALUE_ALWAYS`). Other output formats are answered as NCBI does. |
| `DIAG_*` (any), `ABORT_ON_THROW`, `DEBUG_STACK_TRACE_LEVEL`, `EXCEPTION_STACK_TRACE_LEVEL`, `NCBI_ABORT_ON_COBJECT_THROW`, `NCBI_ABORT_ON_NULL`, `LOG_TRUNCATE`, `LOG_NOCREATE`, `LOG_PERFLOGGING`, `THREAD_STACK_SIZE` | every program | Parameters of NCBI's diagnostics, error handling and threads, which can change its messages or exit status (some only on an error). |
| `NCBI_CONFIG_*` other than `NCBI_CONFIG_PATH` and the entries listed above | every program | Set entries of NCBI's configuration. |

An environment prepared for NCBI BLAST+ that sets one of the variables marked "every program" (for example `DIAG_POST_LEVEL` in a container) stops every LOSAT program before it runs; the other rows stop only the programs and runs named in them.

Matrix files: for a `-matrix` that NCBI BLAST+ does not have built in (which only an ungapped `tblastx`, `blastx -ungapped -comp_based_stats 0`, `blastp -ungapped -comp_based_stats 0` or `tblastn -ungapped -comp_based_stats 0` search takes), NCBI BLAST+ looks for a file of that name in upper case in the directories of `NCBI_DATA_PATH` and `[NCBI] Data`, then in `BLASTMAT` and its `nt` subdirectory, then in `data` under the current directory. When it finds one that it can open, NCBI BLAST+ 2.17.0 reads it as a matrix and then stops abnormally (or with another error code when the file is not a matrix), so LOSAT stops with a message before the search. When it finds none, or cannot open it, both stop with NCBI's "Unknown error code -1".

## Variables that LOSAT does not read

- `IGDATA`, which only NCBI's `igblastn` and `igblastp` read.
- `BLAST_MT_QUERY_BATCH_SIZE`, the query batches of NCBI's searches split by queries (`-mt_mode 1`, or NCBI's own choice at `-num_threads` 2 or more). LOSAT does not split them that way: its output is NCBI's `-num_threads 1` output.

## WebAssembly

When the host sets `PWD`, the WASI command first changes to that directory, so relative file names are read from it, and `makeblastdb`'s `New DB name:` is the absolute path in it, as natively (without `PWD` the directory is `/`). The WASI command reads `TZ`, `TZDIR` and the zone files (`/etc/localtime`, `/usr/share/zoneinfo`) as on Unix, from the environment and the directories that its host gives it; without them `makeblastdb`'s times and the database's date are UTC. See [Platforms and WebAssembly](platforms.md) for what a host must pass.

The Web GUI gives the engine an environment of its own: it has only the directory variable of `-window_masker_taxid` (see above), and none of the other variables of this page applies there. The engine reads the configuration files that the page gives it, the file of `-conffile` and a `<program>.ini` or `.ncbirc` given under that name, as files of the command line's working directory. The Web GUI reads a `-logfile` name as the command line on Linux does, but carries no list of the system's files: a name that climbs with `..` out of a file or out of a missing name in the system's directories, such as `-logfile /etc/passwd/../../tmp/x`, gives a file among the run's files in the Web GUI, where the command line refuses the path and writes the messages to the standard error stream. The Web GUI reads an `-export_search_strategy` name the same way, against the run's files, and gives the search strategy file among the run's files under the command line's name of the file (`s.asn` for `./s.asn`), with the command line's bytes. A name in a folder that the run's files do not have (`out/s.asn` without a folder `out`) or that names a folder of them is NCBI's error, as on the command line run in that folder. Where the command line could not keep that file apart from the run's other files and streams, or would empty or make a file of the run, the Web GUI refuses the name: the standard output (`/dev/stdout`, `/dev/fd/1`), which is the report there, the standard error stream, the system's devices and files, the file of `-query`, `-subject`, `-in_pssm`, `-in_msa`, `-out_pssm`, `-out_ascii_pssm`, `-logfile` or `-window_masker_db`, a file that the search reads or checks in the run's folder, also where it is missing (a file of the database of `-db` or of BLASTN's `-filtering_db`, the list of `-gilist`, `-negative_gilist`, `-ipglist`, `-negative_ipglist`, `-seqidlist`, `-negative_seqidlist`, `-taxidlist` or `-negative_taxidlist`), and any other file of the run's files (such as `.ncbirc` or the file of `-conffile`). `/dev/null` gives no file, as on the command line.

## Development variables (unsupported)

Every variable whose name starts with `LOSAT_` (`LOSAT_X_*`, `LOSAT_TIMING`, `LOSAT_X_STATS`, `LOSAT_X_REACH`, `LOSAT_DEBUG_*`, `LOSAT_TRACE_*` and the others), `BLEMIR_DEBUG` and `BLEMIR_DEBUG_WINDOW` are development variables, as are the cargo features `diagnostics`, `xstats`, `blastx-worker-probe` and `tblastx-wasm-scalar`. They are unsupported, may change or disappear without notice, and are outside what the release tests cover: the release's results hold with none of them set and a normal build. Do not set them in production runs.

- The performance switches (`LOSAT_X_*`) choose between ways of computing the same output. The defaults (no variable set) are what the release tests cover.
- The debug and trace variables write to standard error, and some of them turn off a faster path (any `LOSAT_DEBUG_*` variable turns off `tblastx`'s pipelined batches).
- A performance switch is on unless it is `0`; most other development variables, `LOSAT_X_STATS` and `LOSAT_TIMING` among them, are on when set to any value, even `0`.
- A diagnostic variable set for a build without the `diagnostics` feature is ignored with a note on standard error.

For developers, [docs/agents/diagnostics.md](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/agents/diagnostics.md) describes each development variable.
