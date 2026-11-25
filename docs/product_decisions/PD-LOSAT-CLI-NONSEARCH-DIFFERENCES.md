# Product Decision: CLI behaviour outside the search results

- Decision ID: `PD-LOSAT-CLI-NONSEARCH-DIFFERENCES`
- Version: 1.33 (pending: the Owner's review of exception 14 and of exception 12's prefixes, line NW2's pull request)
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
  E1w2-1 to E1w2-11, Owner-delegated, the same standing instruction); 1.22 2026-11-20 (exception
  11, BLASTN's progress bar on a terminal, recorded for the v0.2.0 release on the Owner's choice of
  2026-11-20; the Owner approved its text the same day, with no words: at 14:24 JST the Owner
  pressed "Ready for review" on pull request #1355, which its description named as the Owner's
  decision on the text; the bar is cleared at the end of each batch's search since this version,
  Owner-delegated, the same standing instruction); 1.23 2026-11-20 (exception 12:
  messages in LOSAT's words; the Owner's instruction of 2026-11-20, line NW, NW DECISIONS 1, 5, 6;
  exception 13: the SAM `@PG CL:` tag, recording `docs/decisions/OWNER_DECISIONS.md` entry 7 of
  2026-10-19, Owner-delegated, which said the difference goes into a product decision; found by
  line MAN's review); 1.24 2026-11-21 (`losat windowmasker` in the scope and in exception 10; line
  NP, item M5a; Owner-delegated, the same standing instruction); 1.25 2026-11-21 (exception 14:
  LOSAT's name and version where LOSAT states the running program's version, in the search outputs
  too, and LOSAT's error prefixes; the Owner's 「Q1=B Q2=B」 of 2026-11-21, line NW2, NW2
  DECISIONS 1; the forms, NW2 DECISIONS 2-4, Owner-delegated, the same standing instruction);
  1.26 2026-11-21 (every search call of the library API runs the application start at its
  beginning, as its tool calls and the command line do; line LIB item E1; Owner-delegated,
  standing instruction 2026-09-29, restated 2026-10-07); 1.27 2026-11-22 (exception 2: the memory
  LOSAT uses with more than 64 threads where NCBI's does not grow with `-num_threads`, recorded
  as a known difference for the v0.2.0 release on the Owner's choice of 2026-11-22, line EB,
  EB DECISIONS 55; the Owner approved its text the same day, with no words: at 22:55 JST the
  Owner pressed "Ready for review" on pull request #1484, which its description named as the
  Owner's decision on the text); 1.28 2026-11-23 (exception 10: blastdbcmd's `-metadata` reports
  a failed write as NCBI's serial stream does, no longer refused; line RG item KC-11b, RG
  DECISIONS KC-11b-1; Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07);
  1.29 2026-11-24 (exception 1 keeps only `-xmlhelp`: argument errors and `-help`/`-h` are
  NCBI's, line RG item K4i-2b, RG DECISIONS K4i-2, K4i-2b-1 and K4i-2b-2; Owner-delegated,
  standing instruction 2026-09-29, restated 2026-10-07); 1.30 2026-11-24 (exception 3
  withdrawn: a failed write of outfmt 6, 7 or 10 aborts as NCBI's; exception 7
  withdrawn: a failed write of a side file of outfmt 13 or 14 ends as NCBI's; a failed write to
  standard error is NCBI's, the run goes on; exception 5 for every format of the search programs
  and for their standard error, its timing of 1.2 restated, the Owner's answer A to line RG's T6-Q1 of 2026-11-17; line RG item
  T6a, RG DECISIONS T6-1, T6-2 and T6-3 of 2026-11-17, Owner-delegated, the same standing
  instruction);
  1.31 2026-11-24 (`losat makeblastdb`'s
  log in the scope and in exception 10; the Owner's choice of 2026-11-24, `Q-KD-1=A`, line RG
  item KD; numbered after line RG item T6a's 1.30, first 1.30);
  1.32 2026-11-25 (`losat makembindex` in the scope and in exception 10; line RG item KD, RG
  DECISIONS KD-mbi-pipe-1; Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07);
  1.33 2026-11-26 (`losat makeprofiledb` in the scope and in exception 10; line RG item KD,
  Owner-delegated, standing instruction 2026-09-29, restated 2026-10-07)
- Status: Accepted by the maintainer on 2026-10-02, in Session S07+++ (E2g), on the items
  that the BLASTN inventory (`docs/evidence/losat_web_e2g/INVENTORY.tsv`, moved to a release
  asset on 2026-10-23: [ARCHIVED.md](../evidence/losat_web_e2g/ARCHIVED.md); actions `S08` and
  `OPEN`) left for a maintainer decision. Plan decision DW-13 in
  [`docs/losat_web_gui_plan.md`](../losat_web_gui_plan.md).

## Scope

The command-line behaviour of every LOSAT program (BLASTN, BLASTP, TBLASTN, TBLASTX,
BLASTX, since version 1.10 `losat rpsbproc`, NCBI's rpsbproc 0.5.1, whose report takes the
place of a search result, since version 1.13 `losat blastdbcmd` and `losat dustmasker`, and
since version 1.16 `losat segmasker`, since version 1.18 `losat convert2blastmask` and since
version 1.24 `losat windowmasker`, whose outputs do, since version 1.31 the log of `losat
makeblastdb`, since version 1.32 `losat makembindex`, and since version 1.33 `losat makeprofiledb`) that is not part of a search result: argument parsing, help text, thread
warnings, since version 1.27 the memory used with more than 64 threads, output-write failures, memory exhaustion and, since version 1.22, BLASTN's progress
bar on a terminal, and since version 1.25 the statement of the running program's version in every
output, search results included (exception 14). Search results, their outfmt
0/6/7 bytes, warnings and errors raised after argument parsing, and exit codes of
those errors remain under the root [`AGENTS.md`](../../AGENTS.md) bit-perfect rule, apart from that
version statement and the words of exceptions 12 and 14.

## Approved exceptions

Each exception is narrow. Nothing outside the listed behaviour may differ from NCBI
BLAST+.

1. **`-xmlhelp`.** `-xmlhelp` (since version 1.20) prints NCBI BLAST+ 2.17.0's XML of the
   program's arguments without the arguments that LOSAT refuses by name and the dependencies that
   name them: the definitions of the arguments that LOSAT runs are NCBI's (the Owner's choice of
   2026-11-13, RG DECISIONS E2x-1). LOSAT Web's ABI v1, which is frozen, keeps the argument
   texts it had ([`docs/web/abi_v1.md`](../web/abi_v1.md) §3). Up to version 1.28 this exception
   also covered LOSAT's own argument errors (exit status 2) and help texts. Since version 1.29
   (line RG item K4i-2b) the programs read their arguments with NCBI's argument reader: an
   argument error prints NCBI's usage and posts (the second post without NCBI's
   `(CArgException::eX) ` part, exception 12) and exits with status 1, and `-help`, `-help-full`
   and `-h` print NCBI's texts in LOSAT's words (exceptions 12 and 14). Errors that NCBI raises
   after parsing the arguments (`BLAST query/options error: ...`, input errors such as a missing
   `-subject`, engine errors) are not covered and must match NCBI.
2. **Thread count with `-subject`.** LOSAT searches with the requested `-num_threads`.
   NCBI reduces a local-subject search to one thread and prints
   `'num_threads' is currently ignored when 'subject' is specified.`, and above the CPU
   count `Number of threads was reduced to N to match the number of available CPUs`
   (`blast_args.cpp:3203-3236`). LOSAT prints neither warning. Every other byte of stdout
   and stderr must equal NCBI's (BLASTN's progress bar on a terminal: exception 11).
   Since version 1.27: with more than 64 threads, LOSAT can use more memory than NCBI in the
   searches where NCBI's memory does not grow with `-num_threads`, such as those with
   `-subject`, which NCBI runs on one thread: measured for v0.2.0 (development build
   `ff07c6a7d` against NCBI BLAST+ 2.17.0 on 64 logical CPUs), LOSAT's median peak resident
   memory was up to 1.21 times NCBI's at 96 threads and up to 1.44 times at 128 threads
   (TBLASTN with `-subject`; also BLASTX, BLASTN and BLASTP with `-subject` and TBLASTN with a
   database of one sequence, but not TBLASTX or `blastn-short` with `-subject`), while LOSAT
   finished 8 to 69 times sooner; the Owner accepted this difference for v0.2.0 on 2026-11-22.
3. **Write failure in tabular output.** When writing outfmt 6 or 7 fails (for example
   `-out /dev/full` or a closed pipe), NCBI aborts; LOSAT reports the write error and
   exits with a non-zero code. For outfmt 0 LOSAT must follow NCBI
   (`BLAST failed to write output: <msg>`, exit code 6 `BLAST_OUTPUT_ERROR`,
   `blast_app_util.hpp:242-254`); this is ported program by program (BLASTN in S07+++b).
   Version 1.30: withdrawn. Every program's report stream fails where NCBI's does
   (`report::ostream::ReportStream`, line RG item T6a): a failed write of outfmt 6, 7 or 10
   aborts with libstdc++'s `terminate called after throwing an instance of
   'std::__ios_failure'` text, as NCBI's (the destructor of `CBlastTabularInfo` flushes the failed
   stream, `tabular.cpp:160-163`). A closed pipe is exception 5.
4. **Memory exhaustion.** LOSAT aborts when an allocation fails. NCBI catches the
   allocation failure and prints `BLAST ran out of memory`, exit code 4
   (`BLAST_OUT_OF_MEMORY`, `blast_app_util.hpp:216-235,255-258`).
5. **Closed pipe** (outfmt 0 in version 1.1, accepted by the maintainer on 2026-10-02 in
   Session S07+++b; every format in version 1.30). NCBI keeps the C runtime's default action for SIGPIPE, so a write to
   a pipe whose reader has closed ends it by the signal (exit status 141 from the shell,
   no message), in every format. The Rust runtime ignores SIGPIPE; restoring it needs a
   native `signal` call that the project's pure-Rust runtime boundary check rejects. LOSAT
   reports the failed outfmt 0 write as NCBI reports other outfmt 0 write failures:
   `BLAST failed to write output`, exit code 6 (outfmt 6/7: exception 3, withdrawn in 1.30).
   Version 1.2: whether a write to a pipe fails depends on when the reader closes it. When
   LOSAT's writes have completed before the reader closes (for example a short outfmt 6 or
   7 report read by `head -c 100`), LOSAT exits with code 0, where NCBI, which writes later
   (at its flushes), is ended by SIGPIPE (141). This timing is part of exceptions 3 and 5
   (round-3 audit of Session S07+++b; exception 3 withdrawn in 1.30, the timing restated below).
   Version 1.30: every format of the search programs (the Owner's answer A to line RG's
   question T6-Q1, 2026-11-17: no crate and no native call; RG DECISIONS T6-3). A write to a
   closed pipe (`EPIPE`; the Rust runtime ignores SIGPIPE) is a failed write of the report's
   stream like any other (`report::ostream::ReportStream`, line RG item T6a), and the run ends as
   NCBI's run ends when that format's write to its standard output fails, where NCBI is ended by
   SIGPIPE (141, no text): outfmt 0 to 4 and 18 with `BLAST failed to write output`, exit code 6
   (outfmt 0 to 4 abort, exit status 134, where the flush at the end of the alignments fails:
   its 802.4 post's `cerr` flushes the failed `cout`); outfmt 5, 8, 9, 11, 12, 15 and 16 (whose
   failures post, and each post's `cerr` flushes the failed `cout`) and outfmt 6, 7 and 10 with
   libstdc++'s `terminate` text and an abort (exit status 134); SAM (17) with exit code 0. A pipe
   that only the writes after the formatter reach (`~CBlastFormat`'s flush and the exit's) is not
   reported (exit code 0), as NCBI does not report a failed device there.
   The timing of version 1.2, as it happens since version 1.30: a write to the pipe fails only
   when it comes after the reader has closed, in LOSAT as in NCBI. LOSAT keeps up to 8,192 bytes
   of NCBI's writes before it makes them together (where one fails, the run ends where NCBI's
   failed write ends it), so a short report reaches the pipe in one write at the end of the run;
   when that write lands before the reader closes, the run sees no failure and exits with code
   0 (`head -c 0` on a one-query BLASTN outfmt 0 report: LOSAT 0 in 2 of 40 runs, code 6 in the
   others; NCBI 2.17.0 141 in 40 of 40; line RG item T6a, review 2 N1, 2026-11-18).
   A closed pipe on standard error (version 1.30) is a failed write of standard error, which
   NCBI's run survives: NCBI's `cerr` that cannot write (`2>/dev/full`) takes no later message,
   and the run goes on to the end of its report and its exit code (0 for a search that ends
   well; NCBI 2.17.0 observed, line RG item T6a, review 2 D2). LOSAT does the same where NCBI
   is ended by SIGPIPE at its first message (exit status 141 from the shell): the report is
   whole and the messages are lost.
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
   Version 1.30: withdrawn. A side file is written as NCBI's object stream and `CNcbiOfstream`
   write it (`report::ostream::write_side_file`, line RG item T6a; NCBI 2.17.0 probed on Linux,
   RG survey T6 D8): a write that fails inside the document is the buffer's write fault, posted
   as `Error: (802.5) [<program>] CObjectOStream: error at byte <N>: stack is empty: write
   fault` (`N`: the bytes the object stream had put into its buffer), and the run stops there
   with exit code 0 (`CATCH_ALL` swallows that `CIOException`; no later report and no master
   file); a failure of the file's last flush alone posts NCBI's chain of `cannot close output
   stream` and the run goes on. The earlier text above recorded an observation that the probe
   does not repeat (NCBI does not write the remaining side files).

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
   rpsbproc's own code for its output file (`e_OutputFileOpenError`, main:91-101), as NCBI's
   search programs use their output-error code 6 for a failed report write (`BLAST failed to
   write output`; exception 3, which had LOSAT's search programs do so in every format, is
   withdrawn in 1.30). As in
   exception 5 (version 1.2), a closed pipe that LOSAT's completed writes never reach ends
   with exit code 0. Nothing is written to standard error after the failure (no progress
   line). A failed write to standard error (its reader has closed: NCBI's rpsbproc ends by
   SIGPIPE, exit status 141 from the shell, PS16a review oracle M1) ends LOSAT at once with exit
   code 4 and no message; so does a failed write of rpsbproc's help (its text, or since version
   1.20 its `-xmlhelp` XML) or version, after the message above. A successful run is identical to NCBI's. The function is `rpsbproc::run`
   (`LOSAT/src/rpsbproc/mod.rs`, `Report`). Release-notes sentence: If `losat rpsbproc` cannot
   write its report, it stops with `rpsbproc failed to write output` (exit status 4), where
   NCBI's rpsbproc 0.5.1 drops the writes and continues with exit status 0.
10. **Closed pipe in blastdbcmd's, dustmasker's, segmasker's, windowmasker's and convert2blastmask's output,
   in makeblastdb's log, in makembindex's `-show_filters`, `-stat` and index outputs that are named pipes, and in every program's version and help**
   (version 1.13, 2026-11-13, Owner-delegated; line RG, item K1, RG DECISIONS K1-8; version 1.14
   the same day, item K1b, K1b-3 and K1b-5; version 1.16 the same day, segmasker, item K2, K2-5;
   version 1.18 the same day, convert2blastmask, item K3; version 1.24, 2026-11-21, windowmasker,
   line NP item M5a; version 1.31, 2026-11-24, makeblastdb's log, the Owner's choice, line
   RG item KD; version 1.32, 2026-11-25, makembindex, Owner-delegated, line RG item KD). A write of NCBI's blastdbcmd 2.17.0,
   dustmasker 1.0.0, segmasker 1.0.0, windowmasker 1.0.0 or convert2blastmask 2.17.0 to a pipe whose
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
   segmasker's `interval` and `fasta`, windowmasker's `interval` and `fasta`, convert2blastmask's `interval`,
   blastdbcmd's retrieval, `-info`, `-list`, `-tax_info`) the later writes are dropped and the
   run goes on to its end with its exit status (0, or the error of a later input or entry, such
   as blastdbcmd's `Skipped` of an id that is not found, exit status 1); where NCBI writes
   through its serial stream, which reports a failed write (dustmasker's `seqloc_*` and
   `maskinfo_*`: exit status 3, 0 or 134 by the case; segmasker's: 1 or 134, the local official
   release with `-out /dev/full`; windowmasker's: 3 for `seqloc_*`, 134 for `maskinfo_*`, NCBI
   2.17.0 with `-out /dev/full` on ws, 2026-11-21; convert2blastmask's `maskinfo_*`), LOSAT refuses it (exit
   status 1). blastdbcmd's `-metadata` is not refused since version 1.28: its serial stream, its buffer
   and the stream under it (the `-out` file's `ofstream`, or `cout`) are NCBI's, so a failed write gives
   NCBI's lines, file bytes and exit status (6 when the flush of the part of the document that the buffers
   still hold fails, else 0; 300 runs of NCBI 2.17.0 on a full device and under file size limits,
   `tests/blastdbcmd_metadata_failures.rs`). A run whose writes all succeed is identical to NCBI's; so
   is a closed pipe that NCBI's writes never reach (a short output read by `head -c 100`).
   Since version 1.31 the same holds for makeblastdb's log, on the standard output or on a
   `-logfile` that is a pipe. The Owner, 2026-11-24 12:44 JST, answered line RG item KD's question
   Q-KD-1 with 「Q-AGVI-CAP=A Q-FM-R2=B Q-D23=A Q-KD-1=A Q-MIXIDX=C Q-GBHELP=A Q-DL-TIE=B
   Q-DL-QLOC=B Q-PS7c-SPOUGE=B Q-PS10e=A」, choosing the option 「閉じたパイプへの最初のログ書き込み失敗で即停止し、ファイルを作らず
   exit 0（PD 例外 10 に `makeblastdb` を追記）」 (stop at the first failed write of the log to the
   closed pipe, make no files, exit status 0; name makeblastdb in exception 10). Every line of
   NCBI's log ends with a flush (`endl`), and the first one that meets the closed pipe ends
   makeblastdb 2.17.0 by SIGPIPE; LOSAT ends the build at the end of that line with exit status 0
   and no message, and nothing more is read, written, removed or reported (no `EndBuild`, no
   close of the writer, no message of `~CBuildDatabase`, no metadata file). A reader that is gone
   when the program starts meets the first line of the log, which the constructor of the build
   writes after it has created the output directory and before `DeleteBlastDb`: the output
   directory is made and nothing else, and an existing database of the name stays as it was (NCBI
   2.17.0, nucleotide and protein, versions 4 and 5, probes of 2026-11-24). A reader that closes
   during the build ends it at a later line, and the files written before that line stay, as
   NCBI's ended process leaves its own. At a line of `EndBuild`, after the close (the reader
   closed at the line "Adding sequences from FASTA"), both leave the volume files without the
   metadata file, the same bytes. At the line "Adding sequences from FASTA" itself, before the
   close, both leave the files of the volume that is still open and LMDB's data files, but not
   the same: LOSAT's volume files hold what its writer has passed on (every sequence but the
   last), where NCBI's file streams still hold part of it (NCBI 2.17.0: empty files for three
   short sequences; for 400 sequences of 100,000 bases the same sequence file and 24,656 of
   LOSAT's 29,976 header bytes), and NCBI's ended process also leaves LMDB's lock files
   (`.ndb-lock`, `.ntf-lock`), which LOSAT's writer removes as it does after an error
   (PD-LOSAT-DATABASE-SEARCH). A build that LOSAT
   refuses at its end (mask data whose order NCBI takes from memory addresses) removes its files
   when the closed pipe meets the first line of `EndBuild`, so that LOSAT leaves no database it
   refuses. With `-logfile` naming a file nothing goes to the standard output, and the build runs
   to its end in both. On Windows the log's write fails without a signal and the build goes on, as
   NCBI's does.
   Since version 1.32 the same holds for `losat makembindex` (probes against NCBI BLAST+ 2.17.0
   makembindex, 2026-11-25). A normal `makembindex` build writes progress (`creating ...`, `done`,
   `removed (empty)`) only to standard error and nothing to standard output, so a closed standard-output
   pipe finishes the full index with exit status 0 in both NCBI and LOSAT. `-show_filters` writes
   filter descriptions to standard output (`mkindex_app.cpp:236-240`): on a closed pipe outside
   Windows, NCBI 2.17.0 is ended by SIGPIPE (141) and LOSAT exits with status 0 and no extra
   message (earlier warnings on standard error stay; a missing database fails with exit status 3
   before writing; `/dev/full` exits with status 0 in both). `-stat` flushes each non-empty word
   count line (`dbindex_factory.cpp:1408-1436`), and `.NN.idx` writes 8 KiB blocks after the
   44-byte header (`dbindex_factory.cpp:1604-1634`). Index outputs (`.NN.idx`, `.map`, `.shd`)
   also write through `Ostream`: when `-stat` meets a closed pipe outside Windows (such as `-stat
   /dev/stdout` with standard output closed, or a pipe whose reader closes mid-run), NCBI 2.17.0 is
   ended by SIGPIPE at that line, leaving the current `.NN.idx` incomplete and writing no `done`,
   no `.map`, no `.shd`, and no later index or database volumes; LOSAT flushes `-stat` at each
   8 KiB `.NN.idx` block boundary, so it stops in the same 8 KiB block with exit status 0 and no
   message, leaving the same incomplete `.NN.idx` bytes (the 44-byte header followed by the
   complete 8 KiB blocks, for both regular files and named pipes) and the same files as NCBI. When
   an index output (`.NN.idx`, `.map`, or `.shd`) is a named pipe whose reader closes early outside
   Windows, NCBI 2.17.0 is ended by SIGPIPE and LOSAT stops at once with exit status 0 and no
   message, leaving the same files and `-stat` progress as NCBI. When a volume has no
   indexed words (so `-stat` receives no write) or `-stat` hits a non-pipe write failure
   (`/dev/full`, or a closed pipe on Windows), the index completes with exit status 0 in both.
   Since version 1.33 the same holds for `losat makeprofiledb` (probes against official NCBI
   BLAST+ 2.17.0 makeprofiledb, 2026-11-26). A fresh `-dbtype rps` build writes nothing to
   standard output and finishes with exit status 0 on a closed standard-output pipe in both NCBI
   and LOSAT. Overwriting an existing database of the same name deletes the old database
   (`s_DeleteMakeprofileDb`) and then writes `Deleted existing BLAST database with identical name.`
   with `endl` to `m_LogFile` (standard output unless `-logfile` names a file, `makeprofiledb.cpp:1331-1333`):
   when that stream is a closed pipe outside Windows (stdout, `-logfile -`, `-logfile /dev/stdout`,
   or a closed named pipe), NCBI 2.17.0 is ended by SIGPIPE (141) before building any new volume,
   and LOSAT stops at once with exit status 0 and no message. When `CWriteDB` writes an
   `ib=...,n4sz=...` line to `cout` (`writedb_convert.cpp:464-466`) on a closed standard-output
   pipe (even with `-logfile`), NCBI 2.17.0 is ended by SIGPIPE inside `AddSequence` without
   closing the writer or running `~CMakeProfileDBApp`, and LOSAT stops there with exit status 0
   and no message. Non-pipe write failures (`/dev/full`) and `-logfile` naming a regular file
   complete with exit status 0 in both.
   Functions: `ncbi_app::files::Ostream` (one stream of `-out` for each program;
   `ends_like_sigpipe`, `flush_before_diagnostic`), `dustmasker::work::run`
   (`Failure::ClosedPipe`), `segmasker::work::run_masks` (`Failure::ClosedPipe`), `windowmasker::work`
   (`Failure::ClosedPipe`), `convert2blastmask::work::run_masks`
   (`Failure::ClosedPipe`), `blastdbcmd::retrieve::process` (`Stop::ClosedPipe`, `Run::post`),
   `blastdbcmd/metadata.rs` (`writedb::metadata::write_metadata`, `writedb::cout`),
   `taxonomy.rs`, `list.rs` and `mod.rs`; `writedb::build_db::endl` (`WriteDbError::ClosedPipe`)
   and `writedb::makeblastdb::CoutLog` for makeblastdb's log; `makembindex::work::run`
   (`Failure::ClosedPipe`), `dbindex::factory::make_index` (`Stop::ClosedPipe`) and
   `dbindex::format::SuperHeader::save` for makembindex; `makeprofiledb::work::run`
   (`RunError::ClosedPipe`, `write_line_or_closed_pipe`, `post_writer`) for makeprofiledb. Tests:
   `tests/makeblastdb_closed_pipe.rs` (`a_closed_pipe_at_the_first_line_makes_no_file`,
   `an_old_database_of_the_name_is_left_as_it_was`, `logfile_builds_with_the_pipe_closed`,
   `a_logfile_that_is_the_closed_pipe_makes_no_file`), `writedb/tests.rs`
   (`a_closed_pipe_ends_the_build_at_its_line`), `tests/makembindex_closed_pipe.rs`
   (`normal_index_builds_complete_with_a_closed_stdout_pipe`,
   `show_filters_help_and_version_on_a_closed_pipe_end_silently`,
   `stat_on_a_closed_pipe_stops_at_once`, `stat_pipe_closed_after_first_line_stops_at_once`,
   `volume_and_idmap_named_fifos_stop_at_once`,
   `volume_fifo_matches_ncbi_partial_bytes_and_stat_progress`), `makembindex/tests.rs`
   (`show_filters_on_a_closed_pipe_ends_silently`), `dbindex/format_tests.rs`
   (`super_header_save_reports_closed_pipe_and_ignores_other_write_errors`),
   `tests/makeprofiledb_closed_pipe.rs` (`fresh_build_completes_with_a_closed_stdout_pipe`,
   `overwrite_on_closed_stdout_pipe_stops_after_deleting_existing_db`,
   `overwrite_with_regular_logfile_or_dev_full_completes`,
   `writer_cout_on_closed_stdout_pipe_stops_without_cleanup_even_with_logfile`,
   `help_version_and_dryrun_on_closed_stdout_pipe`) and `makeprofiledb/work_tests.rs`
   (`write_line_or_closed_pipe_maps_broken_pipe_only_off_windows`). Release-notes sentence: When the reader of `losat
   blastdbcmd`'s, `losat dustmasker`'s, `losat segmasker`'s, `losat windowmasker`'s or `losat convert2blastmask`'s output, of `losat makeblastdb`'s log, of `losat makembindex`'s `-show_filters`, `-stat`, or index outputs that are named pipes, or of `losat makeprofiledb`'s log or writer output closes the pipe early, LOSAT stops at once with
   exit status 0, where NCBI's programs are ended by SIGPIPE (141).
11. **BLASTN's progress bar on a terminal** (version 1.22, 2026-11-20, recorded for the v0.2.0
   release on the Owner's choice of 2026-11-20; the Owner approved this text the same day, with no
   words: at 14:24 JST the Owner pressed "Ready for review" on pull request #1355, which its
   description named as the Owner's decision on the text). `losat blastn` on the command line draws
   a progress bar on its standard error while it searches a query batch with `-num_threads` 2 or
   more against two subjects or more (the records of `-subject`, or the sequences of `-db`); the
   serial WASI command runs on one thread and draws none. The bar is one line, redrawn in place: a
   spinner, the elapsed time, a bar of 40 columns and the number of subjects traced over the number
   of subjects (`19/33`). When the traceback of the batch ends, or fails, the bar is cleared before
   anything else is written (the batch's report, or an error): its line is blanked and the cursor
   goes back to the start of that line, so nothing of the bar stays on the terminal and what
   follows appears where it would appear without the bar. NCBI's blastn writes nothing to its
   standard error while it searches its query batches (`blastn_app.cpp:277-279`). The bar is drawn
   only when the standard error is a terminal (a pseudo-terminal too, such as the one of `script`)
   and, outside Windows, `TERM` is set and is not `dumb` (the check of the `indicatif` crate:
   `console::Term::is_term` and `console::is_dumb`); a standard error sent to a file or a pipe
   receives no byte of it. A failed write of the bar is ignored. The bytes of the standard output,
   the `-out` file and the exit status are the same with and without the bar. The library API (and
   pylosat), the Web engine and the other programs draw no bar: only the command line's
   `blastn::run` sets `BlastnArgs::progress_bar` (`LOSAT/src/algorithm/blastn/args.rs:798`,
   `blast_engine/run.rs:8018`); the library and the Web engine search through `run_local` and
   `run_local_database`, which leave it off, and ABI v1 sets it off (`LOSAT/src/web_api.rs:364`).
   The bar is made by `batch_progress_bar` (`blast_engine/subject_env.rs:49-57`: the standard error
   as its target, `ProgressFinish::AndClear` as its end, which also clears a bar dropped
   unfinished), moves on by one for each subject traced
   (`blast_engine/traceback_subject.rs:292-294` and `418-420`) and is finished at
   `blast_engine/run.rs:13533-13535`. Tests: `the_end_of_a_batch_clears_the_bar` and
   `a_search_that_stops_with_an_error_clears_the_bar` (`blast_engine/subject_env.rs`) draw the bar
   on a terminal kept in memory and check that nothing of it stays; no test runs it on a real
   terminal. Up to version 1.21 (LOSAT v0.1.0 and dev before this version) the finished bar stayed
   on the terminal, one line per query batch. Release-notes sentence: On a terminal, `losat blastn`
   with `-num_threads` 2 or more and two subjects or more shows a progress bar on the standard
   error while it searches, where NCBI BLAST+ 2.17.0 shows none; the bar is cleared when the search
   of each query batch ends, and a standard error sent to a file or a pipe gets none of it.

12. **Messages in LOSAT's words** (version 1.23, 2026-11-20). The Owner, 2026-11-20 (verbatim):
   「てかメッセージとかで「BLAST+のマニュアルを参照してください」とかNCBIの文言がそのまま出てきてるやつ、全部LOSATに直しといて」;
   the prefixes stayed (the Owner's 「すべて推奨通り」 to keeping `BLAST query/options error:` and its
   siblings, the same day) until the Owner's 「Q1=B Q2=B」 of 2026-11-21 (exception 14). Where NCBI's messages name NCBI's product or show the internals of
   NCBI's C++ Toolkit, LOSAT's say it in LOSAT's words, in standard error, in the library's and the
   Web engine's errors and in an archive's messages (`Blast4-error.message`):
   - the footer of the query, input, filtering, scoring, formatting, strategy, archive and
     `blastdbcmd` taxonomy errors, `Please refer to the BLAST+ user manual.`, reads
     `Please refer to the LOSAT user manual.`;
   - an exception report (`CException::what()`, NCBI's `NCBI C++ Exception:` and one
     `T0 "<NCBI's build tree>/<file>", line <n>: Error: (<class>::<code>) <function>() - <message>`
     line for each exception of the chain) keeps the messages only: the first after the text
     before the report, each further one on its own line with NCBI's indentation; a post of an
     exception (`Error: (<class>::<code>) <message>`) drops the class and code, and inside a line
     (`Application's initialization failed (<class>::<code>) <message>`) they become `: `, as does a
     report glued to the text before it (NCBI's `<message>NCBI C++ Exception:` posts, which
     read `<message>: <message>`);
     `(m_Pos = <n>)` reads `(position <n>)`;
   - LOSAT's refusals say `the option -logfile` where they said `the NCBI C++ Toolkit option
     -logfile`, and LOSAT Web's notices no longer call the running program BLAST+.
   NCBI's diagnostic codes such as `(106.15)`, the message bodies (also where they name an NCBI
   function, such as `CSeq_align::GetSeqRange(): Row is empty`), every error condition, exit status
   and search result stay NCBI's; since version 1.25 the error category prefixes, the version
   statement of the search results, `-version`, `-help` and `-xmlhelp` follow exception 14. LOSAT's refusals that name
   `NCBI BLAST+` speak of NCBI's program and stay. One producer: `LOSAT/src/wording.rs`
   (`USER_MANUAL_FOOTER`; `losat_text`, which `ncbi_exception::NcbiException` applies to the
   report it builds as NCBI's, as do its posts, `post_line_in`). The comparisons with NCBI's
   freezes apply the same function to the freeze (`losat_output_bytes`; AGENTS.md rule 8).
   Release-notes sentence: Error messages refer to the LOSAT user manual and no longer show NCBI's
   C++ exception reports (source paths, class names); the `BLAST ... error:` prefixes and exit
   statuses are unchanged.

14. **LOSAT's name and version, and LOSAT's error prefixes** (version 1.25, 2026-11-21). The Owner,
   2026-11-21, on a LOSAT Web screenshot whose failed run kept an outfmt 0 output that starts
   `BLASTP 2.17.0+` under `Run 1 failed: BLAST engine error: ...` (verbatim, as session
   losat-baselines-58 recorded it): 「エラーメッセージにBLASTP 2.17.0+ とか書いてあるけどなんで？LOSATでしょ？…LOSAT全体から駆逐すべきでは？」,
   then 「Q1=B Q2=B」: Q1 B, the first line of the outputs names LOSAT (no longer the same bytes as
   NCBI's), and Q2 B, the `BLAST ... error:` prefixes in LOSAT's words. The forms
   (Owner-delegated, NW2 DECISIONS 2-4):
   - Wherever LOSAT states the running program's version it states LOSAT and LOSAT's version
     (the crate's, `0.2.0` in v0.2.0). After a program's name: `BLASTP (LOSAT 0.2.0)` where NCBI
     writes `BLASTP 2.17.0+`: outfmt 0's first line, outfmt 7's `# BLASTP (LOSAT 0.2.0)`, outfmt 5's
     `BlastOutput_version`, outfmt 13 and 15's JSON `version`, outfmt 14 and 16's `<version>`, the
     XML that rpsbproc reads and its `SESSION` line, and the descriptions of NCBI's frozen `-help`
     and `-xmlhelp` texts (`Protein-Protein BLAST (LOSAT 0.2.0)`); the program's name stays first,
     so a reader that takes the first word as the program (Biopython's outfmt 7 reader) still reads
     it. Alone: `LOSAT 0.2.0` where NCBI writes `2.17.0+` alone: an archive's (outfmt 11) and a
     search strategy's `ident`, SAM's `@PG VN:`, the first line of `-version`
     (`blastn: LOSAT 0.2.0`; dustmasker, segmasker and rpsbproc too, where NCBI's state their own
     `1.0.0` and `0.5.1`) and the program's version in `-version-full-xml` and `-version-full-json`,
     blastdbcmd's `version LOSAT 0.2.0`. rpsbproc's banner reads
     `Post-RPSBLAST Processing Utility (LOSAT 0.2.0)` (NCBI's `v0.5.1`). The program's `<version>`
     in `-xmlhelp` holds LOSAT's version number, `0.2.0`, as `-version-full-xml`'s numbers do
     (NCBI's `2.17.0`; dustmasker's and segmasker's `1.0.0`, rpsbproc's `0.5.1`).
   - The structure of the XML and JSON outputs stays NCBI's: element and member names, their order
     and nesting, the DOCTYPE (`-//NCBI//NCBI BlastOutput/EN`, `NCBI_BlastOutput.dtd`),
     `xmlns="http://www.ncbi.nlm.nih.gov"` and `xs:schemaLocation` (`NCBI_BlastOutput2.xsd`), which
     name the format's definition, not the program that ran (Biopython's `Bio.Blast` parser checks
     them).
   - The prefixes: `LOSAT query/options error:`, `LOSAT engine error:`, `LOSAT database error:`,
     `LOSAT database creation error:` (makeblastdb), `LOSAT options error:` and
     `LOSAT query error:`, where NCBI writes `BLAST query/options error:`, `BLAST engine error:`,
     `BLAST Database error:`, `BLAST Database creation error:`, `BLAST options error:` and
     `BLAST query error:`; the same length, so an archive's wrapped message keeps NCBI's wraps.
   - Kept: the `Reference:` citations, `BLAST` as the name of the method or the format
     (`# BLAST processed 2 queries`, `BLAST database`), sentences about NCBI's program
     (`output compatible with blast 2.17.0`, `NCBI BLAST+ 2.17.0 aborts`), NCBI's diagnostic codes,
     message bodies (also `BLAST failed to write output`, `BLAST ran out of memory`) and exit
     statuses, and every other byte of the outputs.
   One producer: `LOSAT/src/wording.rs` (`LOSAT_VERSION`, `REPORT_VERSION`, `ERROR_PREFIXES`); NCBI's
   frozen `-help` and `-xmlhelp` texts are printed through its mapping. The comparisons with
   NCBI's freezes apply the same mapping to the freeze (`losat_output_bytes`; AGENTS.md rule 8;
   `python3 -I LOSAT/tests/losat_wording.py map-output <file>` for an ad-hoc comparison).
   Release-notes sentence: Outputs name LOSAT and its version where NCBI BLAST+ names its own
   (`BLASTP (LOSAT 0.2.0)` in place of `BLASTP 2.17.0+`, `-version` prints `blastn: LOSAT 0.2.0`),
   and error messages start with `LOSAT engine error:` and the other `LOSAT ... error:` prefixes in
   place of NCBI's `BLAST ... error:`; every other byte of the search results is unchanged.

13. **SAM `@PG CL:`** (version 1.23, 2026-11-20, Owner-delegated, recording the decision of
   2026-10-19, `docs/decisions/OWNER_DECISIONS.md` entry 7). BLASTN `-outfmt 17` writes in the
   `CL:` tag of the `@PG` line the command line LOSAT was run with, each argument (the program
   name first) followed by one space, as NCBI's `GetCmdlineArgs` writes its own
   (`LOSAT/src/cli.rs` `command_line_bytes`); a library call, which has no command line, writes no
   `CL:` tag (`LOSAT/src/report/sam.rs` `SamDocument::finish`). The rest of the document is NCBI's
   byte for byte; the fixture comparisons leave out this tag only (`strip_sam_cl` of
   `LOSAT/tests/db_search_fixtures.py`, `LOSAT/tests/outfmt17_blastn.rs`).

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
  (`ncbiargs.cpp:2866-2872`); a word after it is an argument error, NCBI's.
  Since version 1.11 `-logfile` and `-conffile` are no longer refused on the command line of
  BLASTN, BLASTP, BLASTX, TBLASTN, TBLASTX, PSIBLAST and RPSBLAST (and `-conffile` of
  makeblastdb and blastdbcmd): they run as in NCBI (`README.md`, NCBI C++ Toolkit options).
  The library API (until version 1.21) and the Web engine (until line RG item E1w-3a) refused them.
  BLASTX refuses the other words by name, as the other programs do. Since version 1.12 `-dryrun` is no longer refused on
  the command line of those seven programs, makeblastdb and blastdbcmd: the arguments are read
  and checked, then nothing runs, exit status 0, as NCBI's default `DryRun`
  (`ncbiapp.cpp:705-707`, 292-296). An argument error is reported as without `-dryrun`
  (NCBI's). Until version 1.28 the parser's refusals of an NCBI option that LOSAT does not run, or
  of a value form it does not read, stayed with `-dryrun` (LOSAT had no definition to check them as
  NCBI does; line RG DECISIONS E1b-2); since version 1.29 LOSAT's refusals of such an option, and
  of a value that is not UTF-8, come after NCBI's reader and after the dry run's end, so a dry run
  of such a command line ends with exit status 0 as NCBI's does (line RG DECISIONS K4i-2b-2 (1));
  `losat rpsbproc`, whose dry run reads its arguments, refuses them first. The refusals of registry
  and environment settings stay (some, such as
  `DIAG_POST_LEVEL`, make NCBI write `DryRun`'s message; DECISIONS E1b-6). The library API
  refused it until version 1.21, and the Web engine until line RG item E1w-3a. Since line RG item
  E1w-3a the Web engine (ABI v2: BLASTN, BLASTP, BLASTX, TBLASTN and TBLASTX) runs `-logfile`,
  `-conffile` and `-dryrun` as the command line does: every `validate` and `run` makes NCBI's
  application start, the file of `-conffile` is a file that the page gives the engine, the `-logfile`
  file comes back as a file of the run, and a `-dryrun` run reads no record and writes no output (a
  `-logfile` name with a directory, refused until line WQA's sweep D11 of 2026-11-20, is read as the command line on Linux reads it where the name's directories exist). Since version 1.21 the library API (and
  pylosat) runs `-logfile`, `-conffile` and `-dryrun` as the command line does (line RG item
  E1w-2): every search or tool call runs the application start at its beginning (since version
  1.26, line LIB item E1; until then a search call ran it when its arguments had one of them or
  where it read NCBI's settings); the `-logfile` file gets what the command line writes there, the registry
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
  program writes its `-help` text, as NCBI's programs write the bytes of `-help`
  (`ncbiapp.cpp:685-697`; NCBI's makeblastdb adds its hidden `-randomize`). Since version 1.29 the
  texts of `-help`, `-help-full` and `-h` are NCBI BLAST+ 2.17.0's frozen texts in LOSAT's words
  (exceptions 12 and 14), and makeblastdb's `-help-full` lists the hidden `-randomize` as NCBI's
  does. The reader stops at `-h`, `-help` and `-help-full` as
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
  DECISIONS E3-1). The library API and the Web engine refuse the versions, `-version` too, but
  in blastdbcmd, dustmasker, segmasker and convert2blastmask, whose library runs print them as
  their `-version` (E3-2), and, since version 1.21, makeblastdb, whose library run prints them
  to its output (line RG DECISIONS E1w2-6). A failed write of a program's
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
