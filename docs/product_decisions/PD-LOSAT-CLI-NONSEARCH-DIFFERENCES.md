# Product Decision: CLI behaviour outside the search results

- Decision ID: `PD-LOSAT-CLI-NONSEARCH-DIFFERENCES`
- Version: 1.4
- Date: 2026-10-02 (1.0); 1.1 the same day (exception 5, Session S07+++b); 1.2 the same day
  (the timing of a closed pipe, exceptions 3 and 5, Session S07+++b); 1.3 2026-10-03
  (exception 6, a standard output closed at the start, Session S08b); 1.4 2026-10-05 (NCBI
  C++ Toolkit words in an option's value, Session S08+b, plan DW-19)
- Status: Accepted by the maintainer on 2026-10-02, in Session S07+++ (E2g), on the items
  that the BLASTN inventory (`docs/evidence/losat_web_e2g/INVENTORY.tsv`, actions `S08` and
  `OPEN`) left for a maintainer decision. Plan decision DW-13 in
  [`docs/losat_web_gui_plan.md`](../losat_web_gui_plan.md).

## Scope

The command-line behaviour of every LOSAT program (BLASTN, BLASTP, TBLASTN, TBLASTX,
BLASTX) that is not part of a search result: argument parsing, help text, thread
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
   `-help` and `-h` print LOSAT's help text. Errors that NCBI raises after parsing the
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

## Decided handling that is not an exception

- **File names that are not UTF-8** (`-query`, `-subject`, `-out`): LOSAT rejects them
  with an explicit "not supported by LOSAT" error instead of printing them differently
  (NCBI prints the raw bytes). BLASTN in S07+++b; the other programs in their sessions.
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
- **FASTA input that NCBI's `CFastaReader` reads differently**: the explicit rejections
  stay (TD-12) until a dedicated session before S17 ports the `CFastaReader` functions
  on the path, together with the adapter's index scan (plan §10).

## Closed earlier items (accepted as recorded)

- The V-PERF judgements of E1a, E1b and E1c (no regression: the alternating 12-sample
  measurements gave ×0.989 to ×1.008 where the 3-sample runs had exceeded +5%).
- The two CLI behaviour differences of E1c: write failures are reported instead of
  being ignored with exit code 0 (item 3 above governs the exit code); the order of
  LOSAT's own diagnostic stderr lines (not NCBI output).
