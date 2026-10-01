# Product Decision: CLI behaviour outside the search results

- Decision ID: `PD-LOSAT-CLI-NONSEARCH-DIFFERENCES`
- Version: 1.0
- Date: 2026-10-02
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

## Decided handling that is not an exception

- **File names that are not UTF-8** (`-query`, `-subject`, `-out`): LOSAT rejects them
  with an explicit "not supported by LOSAT" error instead of printing them differently
  (NCBI prints the raw bytes). BLASTN in S07+++b; the other programs in their sessions.
- **`.ncbirc`**: LOSAT rejects, with an explicit error, an NCBI configuration file that
  NCBI would read and that sets a key changing the program's output (for BLASTN, for
  example `[BLAST] LONG_SEQID`; the list comes from the NCBI source). A `.ncbirc` that
  only sets keys without an output effect (such as `BLASTDB`) is accepted. BLASTN in
  S07+++b.
- **FASTA input that NCBI's `CFastaReader` reads differently**: the explicit rejections
  stay (TD-12) until a dedicated session before S17 ports the `CFastaReader` functions
  on the path, together with the adapter's index scan (plan §10).

## Closed earlier items (accepted as recorded)

- The V-PERF judgements of E1a, E1b and E1c (no regression: the alternating 12-sample
  measurements gave ×0.989 to ×1.008 where the 3-sample runs had exceeded +5%).
- The two CLI behaviour differences of E1c: write failures are reported instead of
  being ignored with exit code 0 (item 3 above governs the exit code); the order of
  LOSAT's own diagnostic stderr lines (not NCBI output).
