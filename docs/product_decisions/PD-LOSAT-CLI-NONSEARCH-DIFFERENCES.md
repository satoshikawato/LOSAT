# Product Decision: CLI behaviour outside the search results

- Decision ID: `PD-LOSAT-CLI-NONSEARCH-DIFFERENCES`
- Version: 1.21
- Date: 2026-10-02 (1.0); 1.1 the same day (exception 5, Session S07+++b); 1.2 the same day
  (the timing of a closed pipe, exceptions 3 and 5, Session S07+++b); 1.3 2026-10-03
  (exception 6, a standard output closed at the start, Session S08b); 1.4 2026-10-05 (NCBI
  C++ Toolkit words in an option's value, Session S08+b, plan DW-19); 1.5 2026-10-25
  (exception 7, a failed write to a side file of outfmt 13 or 14, recorded for the v0.2.0 release;
  Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07); 1.6 2026-11-03
  (exception 8, a `-query` or `-subject` named `.gz`; the Owner asked for `.gz` inputs on
  2026-11-03, the form is Owner-delegated); 1.7 2026-11-04 (exception 8 also for `.bz2` and
  `.zst`, the Owner's ask of 2026-11-04); 1.8 2026-11-04 (exception 8: the inflated copy is
  made beside the input, the Owner's ask of 2026-11-04; the fallback to `TMPDIR` is
  Owner-delegated); 1.9 2026-11-08 (exception 8: a `.gz` file is read as makeblastdb's
  decompressor reads it, its trailer and header CRC not checked; line GZ); 1.10 2026-11-11
  (`losat rpsbproc` in the scope, exception 9, a failed write of rpsbproc's report; line PS,
  item PS16a, DECISIONS 27 of its task record; Owner-delegated, standing instruction
  2026-09-29, restated 2026-10-07); 1.11 2026-11-12 (`-logfile` and `-conffile` run on the command
  line as in NCBI, line RG item E1, after the Owner's ask of 2026-11-12 to add the features NCBI
  has and LOSAT lacked); 1.12 2026-11-13 (`-dryrun` runs on the command line as in NCBI, line RG
  item E1b, the same ask); 1.13 2026-11-13 (`losat blastdbcmd` and `losat dustmasker` in the
  scope, exception 10, a closed pipe; line RG, item K1, RG DECISIONS K1-8; Owner-delegated,
  standing instruction 2026-09-29, restated 2026-10-07); 1.14 2026-11-13 (exception 10:
  blastdbcmd's retrieval goes on after a failed write and its `-metadata` refuses one; line RG,
  item K1b, RG DECISIONS K1b-3; Owner-delegated, the same standing instruction); 1.15 2026-11-13
  (`-logfile`, `-conffile` and `-dryrun` run on the command line of `losat rpsbproc` as in NCBI's
  rpsbproc 0.5.1, line RG item E1c, the ask of 2026-11-12); 1.16 2026-11-13 (`losat segmasker`
  in the scope and in exception 10; line RG, item K2, RG DECISIONS K2-5; Owner-delegated, the
  same standing instruction); 1.17 2026-11-13 (`-help-full` runs in every program as `-help`, and
  the reader stops at a help flag as NCBI's; line RG item E2, the ask of 2026-11-12);
  1.18 2026-11-13 (`losat convert2blastmask` in the scope and in exception 10; line RG, item K3;
  Owner-delegated, the same standing instruction); 1.19 2026-11-13 (the full versions run in every
  program; a failed write of a version or help is dropped as NCBI's, and exception 10 covers a
  closed pipe there; line RG item E3, RG DECISIONS E3-1 to E3-4, Owner-delegated, the same
  standing instruction); 1.20 2026-11-14 (`-xmlhelp` runs on the command line of every program:
  NCBI's XML without the arguments that LOSAT refuses by name, exception 1; line RG item E2x,
  the Owner's choice of 2026-11-13, RG DECISIONS E2x-1); 1.21 2026-11-14 (the library API runs
  `-logfile`, `-conffile` and `-dryrun` as the command line does; line RG item E1w-2, RG DECISIONS
  E1w2-1 to E1w2-10, Owner-delegated, the same standing instruction)
- Status: Accepted by the maintainer on 2026-10-02, in Session S07+++ (E2g), on the items
  that the BLASTN inventory (`docs/evidence/losat_web_e2g/INVENTORY.tsv`, moved to a release
  asset on 2026-10-23: [ARCHIVED.md](../evidence/losat_web_e2g/ARCHIVED.md); actions `S08` and
  `OPEN`) left for a maintainer decision. Plan decision DW-13 in
  [`docs/losat_web_gui_plan.md`](../losat_web_gui_plan.md).

## Scope

The command-line behaviour of every LOSAT program (BLASTN, BLASTP, TBLASTN, TBLASTX,
BLASTX, since version 1.10 `losat rpsbproc`, NCBI's rpsbproc 0.5.1, whose report takes the
place of a search result, since version 1.13 `losat blastdbcmd` and `losat dustmasker`, and
since version 1.16 `losat segmasker` and since version 1.18 `losat convert2blastmask`, whose
outputs do) that is not part of a search result: argument parsing, help text, thread
warnings, output-write failures and memory exhaustion. Search results, their outfmt
0/6/7 bytes, warnings and errors raised after argument parsing, and exit codes of
those errors remain under the root [`AGENTS.md`](../../AGENTS.md) bit-perfect rule.

## Approved exceptions

Each exception is narrow. Nothing outside the listed behaviour may differ from NCBI
BLAST+.

1. **Argument syntax errors and help.** An error that the argument parser detects
   (unknown option, missing value, a value that does not parse, a value outside a
   declared range such as `-word_size 3`) prints LOSAT's parser message and exits with
   code 2; NCBI prints its USAGE text and the `CArgException` and exits with code 1.
   `-help`, `-help-full` (since version 1.17) and `-h` print LOSAT's help text. `-xmlhelp`
   (since version 1.20) prints NCBI BLAST+ 2.17.0's XML of the program's arguments without the
   arguments that LOSAT refuses by name and the dependencies that name them: the definitions of
   the arguments that LOSAT runs are NCBI's (the Owner's choice of 2026-11-13, RG DECISIONS
   E2x-1). Errors that NCBI raises after parsing the
   arguments (`BLAST query/options error: ...`, input errors such as a missing
   `-subject`, engine errors) are not covered and must match NCBI.
2. **Thread count with `-subject`.** LOSAT searches with the requested `-num_threads`.
   NCBI reduces a local-subject search to one thread and prints
   `'num_threads' is currently ignored when 'subject' is specified.`, and above the CPU
   count `Number of threads was reduced to N to match the number of available CPUs`
   (`blast_args.cpp:3203-3236`). LOSAT prints neither warning. Every other byte of stdout
   and stderr must equal NCBI's.
3. **Write failure in tabular output.** When writing outfmt 6 or 7 fails (for example
   `-out /dev/full` or a closed pipe), NCBI aborts; LOSAT reports the write error and
   exits with a non-zero code. For outfmt 0 LOSAT must follow NCBI
   (`BLAST failed to write output: <msg>`, exit code 6 `BLAST_OUTPUT_ERROR`,
   `blast_app_util.hpp:242-254`); this is ported program by program (BLASTN in S07+++b).
4. **Memory exhaustion.** LOSAT aborts when an allocation fails. NCBI catches the
   allocation failure and prints `BLAST ran out of memory`, exit code 4
   (`BLAST_OUT_OF_MEMORY`, `blast_app_util.hpp:216-235,255-258`).
5. **Closed pipe in outfmt 0** (version 1.1, accepted by the maintainer on 2026-10-02 in
   Session S07+++b). NCBI keeps the C runtime's default action for SIGPIPE, so a write to
   a pipe whose reader has closed ends it by the signal (exit status 141 from the shell,
   no message), in every format. The Rust runtime ignores SIGPIPE; restoring it needs a
   native `signal` call that the project's pure-Rust runtime boundary check rejects. LOSAT
   reports the failed outfmt 0 write as NCBI reports other outfmt 0 write failures:
   `BLAST failed to write output`, exit code 6 (outfmt 6/7: exception 3).
   Version 1.2: whether a write to a pipe fails depends on when the reader closes it. When
   LOSAT's writes have completed before the reader closes (for example a short outfmt 6 or
   7 report read by `head -c 100`), LOSAT exits with code 0, where NCBI, which writes later
   (at its flushes), is ended by SIGPIPE (141). This timing is part of exceptions 3 and 5
   (round-3 audit of Session S07+++b).
6. **Standard output closed at the start** (version 1.3, accepted by the maintainer on
   2026-10-03 in Session S08b, plan DW-17). When a program starts with its standard output
   closed (`>&-`), NCBI's first write to `cout` fails: outfmt 0 reports `BLAST failed to
   write output`, exit code 6; outfmt 6 and 7 end on an uncaught `std::ios_base::failure`
   (abort, exit status 134). Before `main`, the Rust runtime opens `/dev/null` on a closed
   standard descriptor, so LOSAT cannot tell it from a `/dev/null` that the caller opened
   for reading and writing (Python's `subprocess.DEVNULL`, Node's `stdio: 'ignore'`); a
   check before `main` needs a native function that the pure-Rust runtime boundary check
   rejects, and a Linux check of `/proc/self/fdinfo/1` failed such callers (Session S08,
   audit (c) round 2, N2). LOSAT writes the report to that `/dev/null` and exits as the
   search ends (0 when it succeeds), in every program and format; the report is discarded
   either way. Evidence `docs/evidence/losat_web_e2b/closed_stdout/`: BLASTN, TBLASTX,
   BLASTP, TBLASTN and BLASTX in outfmt 0, 6 and 7 (LOSAT exit 0 without stderr in all
   15; NCBI exit 6 in the 5 outfmt 0 runs and 134 in the 10 others).
7. **Write failure in a side file of outfmt 13 or 14** (version 1.5, recorded 2026-10-25;
   the behaviour is older: it is the decision numbered 35 S2 in the outfmt line's decision
   record, which is not in this repository, and is quoted in the code comment at
   `LOSAT/src/report/blast_xml2/framing.rs:104-111`). `-outfmt 13` and `-outfmt 14` write one
   side file per report (`<base>_<n>.json`, `<base>_<n>.xml`) and a master file at `-out`
   that lists them (BLASTN, BLASTP, TBLASTN, TBLASTX and BLASTX). When a write to a side file
   fails after the file was opened (a full disk; `-out` next to a link to `/dev/full`), LOSAT
   stops the run with the message and the exit status of a failed outfmt 0 write
   (`BLAST failed to write output`, exit code 6): the side files written before stay and the
   master file stays empty. NCBI's `BlastXML2_FormatReport(data, file_name)`
   (`algo/blast/format/blastxml2_format.cpp:767-778`) opens the file and writes the object
   through `CObjectOStreamXml` without checking the stream afterwards; the NCBI behaviour
   recorded in the code comment is that it reports the failed flush on the error stream and
   goes on, writes the remaining side files and the master file, and exits with 0 (TBLASTX:
   an empty master file, also exit 0). That NCBI behaviour was observed, not derived from the
   source above, and its observation record is outside the repository. A side file that
   cannot be opened is not part of this exception: LOSAT gives NCBI's
   `Command line argument error: Cannot open output file`, exit code 1
   (`framing.rs:94`, from `blastxml2_format.cpp:773-774`). All five programs make the same
   error from the same function (`side_file_write_error`, `framing.rs:111`; call sites
   `framing.rs:353`, `LOSAT/src/algorithm/blastx/blast_xml.rs:428`). The only test in the
   repository is BLASTX's `a_side_file_write_error_at_the_release_is_the_write_error`
   (`LOSAT/src/algorithm/blastx/blast_xml.rs:1029`, Linux, `/dev/full`), which checks exit 6,
   the message and the empty master file. A successful run is identical to NCBI's; this is
   an error path, so no search result differs. Release-notes sentence: If a side file of
   `-outfmt 13` or `14` cannot be written, LOSAT stops with `BLAST failed to write output`
   (exit status 6), where NCBI BLAST+ 2.17.0 reports it on standard error and continues with
   exit status 0.

8. **A `-query` or `-subject` whose name ends in `.gz`, `.bz2` or `.zst`** (version 1.6,
   2026-11-03, `.gz`; version 1.7, 2026-11-04, `.bz2` and `.zst`; version 1.8, 2026-11-04, the
   copy beside the input; line GZ, DECISIONS 1, 3, 5, 8, 9 and 12). LOSAT inflates the file
   through gzip, bzip2 or zstd before reading it (into an unnamed temporary file beside the
   input, in its directory, or under `TMPDIR` when that directory refuses the file or has no
   room for it; the copy is then read exactly as the plain file would be: the same records,
   messages and line numbers), concatenated members, streams and
   frames (bgzip, pbzip2, pzstd) included; the name alone decides, in any case, as it does for
   NCBI's `makeblastdb -in` (`makeblastdb.cpp:1037-1039`), and `-` is read as it comes. A `.gz`
   file is read as makeblastdb's decompressor reads it (version 1.9): neither a member's trailer
   (CRC-32 and length) nor its header CRC is checked, a trailer cut short reads as the whole
   member, and a file that does not start with a gzip header is read as a zlib stream (the rules
   are listed under `-input_type` in `PD-LOSAT-DATABASE-SEARCH`). NCBI's
   search programs open the file as text (`CNcbiIfstream`) and stop at its first line:
   `BLAST query error: CFastaReader: Near line 1, there's a line that doesn't look like
   plausible data, but it's not marked as defline or comment.`, exit code 1. A file that is
   not what its name says, or is cut short (for `.gz`: where makeblastdb's decompressor fails
   too), stops LOSAT before the search, with
   `Error: -query: 'x.fa.gz' is named .gz but is not gzip data (...); LOSAT reads a file named
   .gz through gzip` (`.bz2` / bzip2, `.zst` / zstd alike; exit code 1; outfmt 11 writes the
   error archive it writes for a file that does not open). A plain file, or a compressed file
   read through `-query -`, is unchanged, so no search result that NCBI can produce differs.
   The decoders are pure Rust (flate2's zlib-rs, a port of zlib's inflate, bzip2's
   libbz2-rs-sys, ruzstd). Tests:
   `LOSAT/tests/cli_smoke.rs` (`compressed_query_and_subject_read_as_the_plain_files` and the
   tests after it), `LOSAT/src/blastinput/compressed.rs` and
   `LOSAT/src/blastinput/input_files.rs`.

9. **Write failure in rpsbproc's report** (version 1.10, 2026-11-11, Owner-delegated; line PS,
   item PS16a). NCBI's rpsbproc 0.5.1 writes its report to `cout` or the `-o` file through
   iostreams (`RpsbProc/src/rpsbproc_main.cpp:727-917`) and never checks the stream: when a
   write fails (a full disk, `-o /dev/full`), the stream sets its bad bit, the later writes
   are dropped, the progress lines go on, and the run exits with 0; a write to a pipe whose
   reader has closed ends it by SIGPIPE (exit status 141 from the shell, no message;
   exception 5). LOSAT stops at the first failed write of the report (a closed pipe is such a
   write) with `rpsbproc failed to write output: <error>` on standard error and exit code 4,
   rpsbproc's own code for its output file (`e_OutputFileOpenError`, main:91-101), as the
   search programs use their output-error code 6 for a failed write (exception 3). As in
   exception 5 (version 1.2), a closed pipe that LOSAT's completed writes never reach ends
   with exit code 0. Nothing is written to standard error after the failure (no progress
   line). A failed write to standard error (its reader has closed: NCBI's rpsbproc ends by
   SIGPIPE, exit status 141 from the shell, PS16a review oracle M1) ends LOSAT at once with exit
   code 4 and no message; so does a failed write of rpsbproc's help (its text, or since version
   1.20 its `-xmlhelp` XML) or version, after the message above. A successful run is identical to NCBI's. The function is `rpsbproc::run`
   (`LOSAT/src/rpsbproc/mod.rs`, `Report`). Release-notes sentence: If `losat rpsbproc` cannot
   write its report, it stops with `rpsbproc failed to write output` (exit status 4), where
   NCBI's rpsbproc 0.5.1 drops the writes and continues with exit status 0.
10. **Closed pipe in blastdbcmd's, dustmasker's, segmasker's and convert2blastmask's output,
   and in every program's version and help**
   (version 1.13, 2026-11-13, Owner-delegated; line RG, item K1, RG DECISIONS K1-8; version 1.14
   the same day, item K1b, K1b-3 and K1b-5; version 1.16 the same day, segmasker, item K2, K2-5;
   version 1.18 the same day, convert2blastmask, item K3). A write of NCBI's blastdbcmd 2.17.0,
   dustmasker 1.0.0, segmasker 1.0.0 or convert2blastmask 2.17.0 to a pipe whose
   reader has closed ends the program by SIGPIPE outside Windows (exit status 141 from the
   shell, no message; exception 5). LOSAT's runtime ignores SIGPIPE (exception 5), so LOSAT
   sees the failed write and ends the program at once, with exit status 0 and no message, in
   every mode and format: nothing more is read, written or reported, as in NCBI's ended process.
   As NCBI's `cerr` is tied to `cout`, the programs flush their output before each line of
   diagnostics, so the closed pipe is met before the line (and with `2>&1` the lines come after
   the output written before them, as NCBI's do).
   Since version 1.19 (line RG item E3, RG DECISIONS E3-3) the same holds for the version, help
   and argument-error texts of every other program but rpsbproc (exception 9): a closed pipe
   drops the text and the program ends with its exit status (0 for a version or a help), where
   NCBI's is ended by SIGPIPE (`cli::main_entry::write_dropping_errors`); a full disk is dropped
   as NCBI drops it.
   On Windows, where NCBI's write fails without a signal, a closed pipe is a failed write like
   a full disk. A failed write that is not a closed pipe (`-out /dev/full`) is NCBI's: where
   NCBI writes to the `ostream` itself (dustmasker's `interval`, `acclist` and `fasta`,
   segmasker's `interval` and `fasta`, convert2blastmask's `interval`,
   blastdbcmd's retrieval, `-info`, `-list`, `-tax_info`) the later writes are dropped and the
   run goes on to its end with its exit status (0, or the error of a later input or entry, such
   as blastdbcmd's `Skipped` of an id that is not found, exit status 1); where NCBI writes
   through its serial stream, which reports a failed write (dustmasker's `seqloc_*` and
   `maskinfo_*`: exit status 3, 0 or 134 by the case; segmasker's: 1 or 134, the local official
   release with `-out /dev/full`; convert2blastmask's `maskinfo_*`; blastdbcmd's `-metadata`:
   exit status 6), LOSAT refuses it (exit status 1). A run whose writes all succeed is identical to NCBI's; so
   is a closed pipe that NCBI's writes never reach (a short output read by `head -c 100`).
   Functions: `ncbi_app::files::Ostream` (one stream of `-out` for each program;
   `ends_like_sigpipe`, `flush_before_diagnostic`), `dustmasker::work::run`
   (`Failure::ClosedPipe`), `segmasker::work::run_masks` (`Failure::ClosedPipe`), `convert2blastmask::work::run_masks`
   (`Failure::ClosedPipe`), `blastdbcmd::retrieve::process` (`Stop::ClosedPipe`, `Run::post`),
   `blastdbcmd/metadata.rs`,
   `taxonomy.rs`, `list.rs` and `mod.rs`. Release-notes sentence: When the reader of `losat
   blastdbcmd`'s, `losat dustmasker`'s, `losat segmasker`'s or `losat convert2blastmask`'s output closes the pipe early, LOSAT stops at once with
   exit status 0, where NCBI's programs are ended by SIGPIPE (141).

## Decided handling that is not an exception

- **File names that are not UTF-8** (`-query`, `-subject`, `-out`): LOSAT rejects them
  with an explicit "not supported by LOSAT" error instead of printing them differently
  (NCBI prints the raw bytes). BLASTN in S07+++b; the other programs in their sessions;
  `losat rpsbproc`'s `-i`, `-o` and `-d` before it writes anything (line PS, item PS16a).
- **`.ncbirc`**: LOSAT rejects, with an explicit error, an NCBI configuration file that
  NCBI would read and that sets a key changing the program's output (for BLASTN, for
  example `[BLAST] LONG_SEQID`; the list comes from the NCBI source). A `.ncbirc` that
  only sets keys without an output effect (such as `BLASTDB`) is accepted. BLASTN in
  S07+++b.
- **NCBI C++ Toolkit words** (`-version`, `-version-full*`, `-dryrun`, `-logfile`,
  `-conffile`, `-xmlhelp`, `-help-full`; decisions D9 and D14 of
  `docs/evidence/losat_web_e2e/AUTHORITY.md` §M): LOSAT rejects them with the explicit
  "not supported by LOSAT's <PROGRAM>" error, in an option's place and, since version 1.4
  (accepted by the maintainer on 2026-10-05 in Session S08+b, plan DW-19), also in an
  option's value (`-out -version`), where NCBI's argv pre-pass (`ncbiapp.cpp:926-1001`)
  acts on them up to `--` (it prints the version and exits 0, drops `-dryrun`, or takes
  the next word as the log or configuration file). BLASTN, BLASTP, TBLASTN, TBLASTX and the
  adapter's `validate` (S08+a). A last `--` changes nothing, as in NCBI
  (`ncbiargs.cpp:2866-2872`); a word after it is an argument syntax error (exception 1).
  Since version 1.11 `-logfile` and `-conffile` are no longer refused on the command line of
  BLASTN, BLASTP, BLASTX, TBLASTN, TBLASTX, PSIBLAST and RPSBLAST (and `-conffile` of
  makeblastdb and blastdbcmd): they run as in NCBI (`README.md`, NCBI C++ Toolkit options).
  The library API (until version 1.21) and the Web engine still refuse them, and BLASTX refuses
  the other words by name, as the other programs do. Since version 1.12 `-dryrun` is no longer refused on
  the command line of those seven programs, makeblastdb and blastdbcmd: the arguments are read
  and checked, then nothing runs, exit status 0, as NCBI's default `DryRun`
  (`ncbiapp.cpp:705-707`, 292-296). An argument error is reported as without `-dryrun`
  (exception 1), and the parser's refusals of an NCBI option that LOSAT does not run, or of a
  value form it does not read, stay (LOSAT has no definition to check them as NCBI does; line RG
  DECISIONS E1b-2), as do the refusals of registry and environment settings (some, such as
  `DIAG_POST_LEVEL`, make NCBI write `DryRun`'s message; DECISIONS E1b-6). The Web engine refuses
  it, and so did the library API until version 1.21. Since version 1.21 the library API (and
  pylosat) runs `-logfile`, `-conffile` and `-dryrun` as the command line does (line RG item
  E1w-2): a call whose arguments have one of them, or that reads NCBI's settings, runs the
  application start; the `-logfile` file gets what the command line writes there, the registry
  and environment settings that the command line refuses are the call's error, and a `-dryrun`
  call reads no record and searches nothing (an empty report per format, no `-out` file). The
  library's tools print the versions to their output. Since version 1.15 `losat rpsbproc` runs the three on the command
  line too: `-logfile` takes rpsbproc's diagnostics (its `LOG_POST_X` messages and the reports
  of an exception; its own `cerr` lines stay on the standard error), `-conffile` is loaded in
  place of `rpsbproc.ini`, and `-dryrun` runs rpsbproc's own `DryRun` after its `Init`
  (`rpsbproc_main.cpp:924-932`: the input, output and data files are opened and loaded, then
  `Dryrun passed, ...` or `Dryrun failed, error subcode = N`, exit status N). LOSAT's refusals in
  that start (a file name, a registry or environment setting) stay with `-dryrun`; those of what
  only the run reads (the XML input, the archive) do not apply, as NCBI's `DryRun` reads no
  archive (line RG DECISIONS E1c-2). Since version 1.17 `-help-full` is no longer refused: every
  program writes its `-help` text (exception 1), as NCBI's programs write the bytes of `-help`
  (`ncbiapp.cpp:685-697`; NCBI's makeblastdb adds its hidden `-randomize`, which LOSAT refuses, so
  LOSAT's text has no hidden option to add). The reader stops at `-h`, `-help` and `-help-full` as
  NCBI's (`ncbiargs.cpp:2840-2852`): an argument after one is not read. `-xmlhelp` (in an option's
  place) and the full versions stay refused, in makeblastdb too (line RG item E2); since version
  1.20 `-xmlhelp` is not refused on the command line (exception 1): it stops the reader as `-help`
  does and prints NCBI's XML without the arguments refused by name (`ncbiapp.cpp:685-690`,
  `ncbiargs.cpp:3870-3950`), and the library API and the Web engine refuse it, in every program
  (line RG item E2x). The library API
  answers `-help-full` as `-help` (the help text in its error), and the Web app and ABI v1's BLASTX
  and TBLASTN refuse it with their other help words. Since version 1.19 the full versions are no
  longer refused: every program prints `-version-full`, `-version-full-xml` and
  `-version-full-json` where NCBI's application reads them (`ncbiapp.cpp:980-997`, before `--`, in
  an option's value too, before the configuration and the arguments), in NCBI's forms
  (`version.cpp:728-884`) with NCBI's first line and the application's version, and LOSAT's
  package as in `-version`: the text's `Package:` line names LOSAT, the XML and JSON package is
  `LOSAT` with its version, and the values that NCBI's build writes about itself
  (`Package-Config`'s configure line, `Build-Signature`, the build dates and revisions,
  `Stable-Components-Version`, `Production-Version`, `Built-As`, rpsbproc's toolkit component)
  are not written; the keys NCBI writes whenever there is a package stay, empty (line RG
  DECISIONS E3-1). The library API and the Web engine refuse the versions, `-version` too, as
  they refuse `-logfile`, `-conffile` and `-dryrun`, but in blastdbcmd, dustmasker, segmasker and
  convert2blastmask, whose library runs print them as their `-version` (E3-2). A failed write of a program's
  version, help or argument error (`-version > /dev/full`) is dropped, as NCBI's `cout` drops it
  (exit status 0 for a version or a help), where LOSAT's earlier write panicked (exit status
  134); a closed pipe there is exception 10; rpsbproc keeps exception 9 (E3-3).
- **FASTA input that NCBI's `CFastaReader` reads differently**: the explicit rejections
  stay (TD-12) until a dedicated session before S17 ports the `CFastaReader` functions
  on the path, together with the adapter's index scan (plan §10).

## Closed earlier items (accepted as recorded)

- The V-PERF judgements of E1a, E1b and E1c (no regression: the alternating 12-sample
  measurements gave ×0.989 to ×1.008 where the 3-sample runs had exceeded +5%).
- The two CLI behaviour differences of E1c: write failures are reported instead of
  being ignored with exit code 0 (item 3 above governs the exit code); the order of
  LOSAT's own diagnostic stderr lines (not NCBI output).
