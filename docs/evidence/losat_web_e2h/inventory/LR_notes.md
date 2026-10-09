# Range LR notes: line reading and the input stream

Table: `LR.tsv` (30 rows). Scratch: `scratch_LR/` (inputs `in/`, runs `out/<label>.cmd/.out/.err/.rc`; `n_*` = NCBI BLAST+ 2.17.0, `L_*` = LOSAT base binary; `plain/` = working directory without `.ncbirc`; `lrmodel.py` = Python model of CStreamLineReader + NcbiGetline + CPushback_Streambuf; `fuzz*.py`, `fixcheck.py`, `show.py` (prints a run), `rows/` = the Python sources of the table). The ~380 runs in `out/` dated before this session were made by the agent that stopped on 2026-10-06 (their scripts `run.sh`, `runs*.sh`, `matrix.sh` still name the old scratch path); this session checked them, added `n_batch_*`, `n_lost_*`, `n_near_*`, `n_cmt_probe*`, `n_task_*`, `n_blastn_noquery_arg`, `n_q_nul_data` and re-ran the model checks. No input starts with a letter or digit before a `>` line or blank line, no `.ncbirc` is needed; no NCBI run could reach the network.

## Call path (all four programs, bl2seq)
1. Arguments: `-query` is `AddDefaultKey ... eInputFile` default `"-"` (cmdline_flags.cpp:47). `CArg_InputFile::x_Open` (ncbiargs.cpp:695-737): `-` -> `std::cin`, else `ifstream.open`; failure -> `CArgException eNoFile` (`Command line argument error: Argument "query". File is not accessible:  `NAME'`, rc 1). The subject stream is opened and fully read in `SetOptions` (blast_args.cpp:2525-2557, `ReadSequencesToBlast` -> `CBlastFastaInputSource` -> `CStreamLineReader` -> `GetAllSeqs`), so subject warnings/errors come first and a bad subject stops before the query is opened. The query stream is `m_InputStream` (blast_args.cpp:3458-3468; no gzip, that is magicblast only).
2. `InitializeSubject`, then `IsIStreamEmpty(query stream)` (blastn_app.cpp:209, tblastx_app.cpp:132, tblastn_app.cpp:213 (no PSSM), blastp_app.cpp:211): true -> `Warning: [prog] Query is Empty!` (stderr), rc 0, before `CBlastFastaInputSource` is built and before `PrintProlog`, so stdout is empty even for outfmt 0/7.
3. `CBlastFastaInputSource(stream)`: `new CStreamLineReader(stream)` (auto EOL). `CStreamLineReaderConverter` (magicblast only) and `CMemoryLineReader` (string constructor, not on the command line) are unreachable.
4. `PrintProlog` (outfmt 0 header), then `for (; !input.End(); ...) GetNextSeqBatch`: `End()` = `AtEOF()`; per record `ReadOneSeq` reads lines through the reader; eEOF ('Expected defline') ends the batch silently; a batch without queries throws `Empty CBlastQueryVector` (rc 3).

## State that persists
Per stream: line number, ungot-line flag, EOL style (decided from the FIRST line only, then possibly switched), the pushback streambufs and the stream's eof/fail bits. One reader per file; the subject and query readers are independent except when both are `-` (same `cin`). Nothing is reset between records or query batches.

## Line reading rules (line_reader.cpp:77-297, ncbistre.cpp:55-176, stream_utils.cpp:379-443)
- First line: `NcbiGetline(.., "\r\n")`; the delimiter read back by `unget(); get()`: CR not followed by LF -> style CR; LF (also the LF of a CRLF) -> style "CRLF" (deferred LF/CRLF decision); no terminator (one-line file) -> unknown.
- Style CR: split at CR; an LF directly after the CR is eaten (CRLF = one break). Style LF/CRLF (`x_AdvanceEOLCRLF`): split at LF; a CR right before the LF is dropped; a lone CR elsewhere in the physical line splits it there, switches the reader to CR style and mixed, and the break at the end of that physical line is LOST (the tail is pushed back without its terminator and is glued to the next line). Mixed style: CR, LF and CRLF each one break, LFCR two (so `\n\r` files have doubled line numbers but correct records).
- Observed: LF, CRLF, CR-only, and each without the final EOL give byte-identical results (blastn and blastp, outfmt 6 and 0); the mixed file LF/CRLF/CR gives q2 with its sequence glued into the title (`Query_2 q2 secondTGGC..`, no data, 20-nt title warning); `\r\r\n` = a break plus an empty line; hyphen-warning line numbers for LF, CRLF, CR, LFCR files: 2,3,5,7,8 / 2,3,5,7,8 / 2,3,5,7,8 / 3,4,7,11,13. Blank/comment/white-space lines before the first defline count (`around line 5` for 3 skipped lines + `>p` + data). `Near line N` for a bad data line is the same in LF, CR, CRLF files (line 4 of the probe).
- No limit on line length (3 MB line, 5058-byte line tested). BOM, NUL, 0x1A, 0xA0, 0xFF are ordinary bytes of a line; a BOM before `>` fails as `Near line 1` (CheckDataLine, range RD).
- Pushback quirk (data loss): when a pushed-back tail is recycled into an existing pushback buffer (<= 256 bytes, no new streambuf) the stream's eofbit is not cleared; if the last, unterminated physical line was read to EOF, `AtEOF()` is true and the tail is never returned (`n_lost_g55`, `n_lost_g55_nl`: the last line `-ACGT..` is never read; with only ONE embedded EOL (`n_lost_a..d`) nothing is lost). The model `lrmodel.py` reproduces these and the fuzz comparison (hyphen-warning lines against model lines) found 0 mismatches in this session's re-runs (`fuzz.py 901/904/905`, `fuzz2.py 902`), and `fixcheck.py` reproduces all 1413 `file_` rows of `LOSAT/tests/unit/blastx_stage_e_io_stream_expected.tsv`.

## `IsIStreamEmpty` and pipes (blast_app_util.cpp:845-875, non-Windows branch 855-874)
`tellg() < 0` -> false at once (pipe, FIFO, tty, or a stream whose fail/eof bit is set). Else `in >> c` with skipws (C isspace: space, \t, \n, \v, \f, \r): none -> true. NUL, 0x1A, 0xA0 are not space. Results (all four programs):
| input as `-query` | NCBI |
| --- | --- |
| empty file, /dev/null, white space/blank lines file, directory, `< empty`, `<&-` | `Warning: [prog] Query is Empty!`, rc 0, stdout empty |
| empty pipe (`: \|`) or empty FIFO | silent, rc 0; outfmt 6/10: nothing; outfmt 7: `# BLAST processed 0 queries`; outfmt 0: header + epilog (`Database: ... Matrix: ... Gap Penalties: ...`) |
| pipe with only blank/white-space lines, comment-only file | `BLAST engine error: Empty CBlastQueryVector`, rc 3 (outfmt 0: header already printed) |
| `-query - -subject -` (file or pipe) | silent, rc 0 (subject takes all of cin; no 'Query is Empty!'), outfmt 0 prints header/epilog |
| subject empty / white space / blank lines / directory / empty pipe | `BLAST engine error: Empty CBlastQueryVector`, rc 3, stdout empty (subject is never tested by IsIStreamEmpty) |
| empty query + empty subject | rc 3 `Empty CBlastQueryVector` (subject first) |
LOSAT-before: identical for everything seekable; rejects the empty-pipe and blank-pipe cases (`an empty query from a stream without a position (such as a pipe) is not supported by LOSAT's <PROGRAM>`, rc 1), both-`-` blank streams, comment-only files; BLASTX rejects every `-`.

## Surprises / for the session
1. A lone CR in an LF file (or a lone LF in a CR file) is a break AND eats the next break (tail glue). A port that treats CR/LF/CRLF uniformly will differ on mixed files (records, titles, line numbers).
2. The pushback eofbit quirk loses the last line in rare mixed-EOL files; keep BLASTX's/stream.rs `push_back`/`failed` handling exactly (checked against a C++ oracle in the BLASTX table, not rewritten).
3. Same bytes, different exit: white-space-only seekable file = rc 0 'Query is Empty!'; as a pipe = rc 3 `Empty CBlastQueryVector`.
4. Empty pipe is a normal, successful run with an empty report (outfmt 0 prints a full header and epilog) - the report writers (AD/RP) must support 'zero queries'.
5. Trailing blank/comment lines never create an extra empty batch (consumed by the last ReadOneSeq), also not with `BATCH_SIZE` below the record length.
6. `-num_threads`/`-mt_mode` with `-subject` do not change reading (bl2seq is single-threaded; warning text is the approved exception).

## Decisions / open questions for the session (recommendations)
- Make the shared reader the only reader for all four programs and both roles; open `-` as a File duplicate of fd 0 (as `open_input` does) but share ONE failed/eof state between the subject and query `FastaStream` when both are `-` (row 1); `stream_is_empty` must return false when `stream_position` fails (done in stream.rs).
- Windows builds: NCBI's IsIStreamEmpty differs (empty pipe is empty). Recommend the Linux semantics everywhere (the oracle) unless the Owner wants Windows parity (row 30).
- Produce a 'no queries' report (empty pipe) and the `Empty CBlastQueryVector` rc 3 path in the callers (rows 7, 8).
- Deliver queries batch by batch from the stream (BLASTX's `emit_fasta_batch` pattern) so that later-batch errors follow earlier output (row 24).
- stream.rs `getline` has an extra fast path over the buffered window; it looks equivalent by source but this inventory did not run the Rust tests: run the BLASTX fixture (`blastx_stage_e_io_stream_expected.tsv`, 12 180 rows) against `stream.rs` when it is wired.

## BLASTX fixture `LOSAT/tests/unit/blastx_stage_e_io_stream_expected.tsv`
Independent C++ oracle output, 9 columns: name, input bytes (hex or `~byte:count` runs), seekable flag, check-empty flag, fault offset (-1 = none), (column 6, not used by the test), `IsIStreamEmpty` result, `|`-joined lines read, number of lines. Row groups: 40 named samples (`lf/cr/crlf/mixed/noeol/space/vtab/all_space/comments/empty` x seekable x check_empty); 20 `line_failure_*` rows (a read fault at offset N, with and without the empty check); 8 201 `exhaustive_*` (short byte strings, 7 108 of them with a fault offset; the alphabet was not inspected); 2 506 `memory_*` (CMemoryLineReader); 1 413 `file_*` (real ifstream, seekable, windows around the 8191-byte buffer and pushback/EOF boundaries). It covers the line splitting, EOL switching, pushback, eofbit, `IsIStreamEmpty` and `CMemoryLineReader`; it does NOT cover line numbers, PeekChar/UngetLine, End()/batching, or stdin sharing. My Python model reproduces all 1 413 `file_` rows.

## Where I stopped
All rows of the brief are written (30). Nothing remains except re-running the Rust tests of `stream.rs` (not allowed in this range).
