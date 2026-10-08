# AD notes: Web adapter (scan, register, run), ABI v2 contract, ABI v1

Sources: base-src `web/adapter/src/{scan,store,lib,run,json,emit}.rs`, `docs/web/abi_v2.md`, `LOSAT/src/web_api.rs`, `web/app/src/{ports/engine.ts,domain/{programs,dataset}.ts,application/{draft,coordinator}.ts,infra/data/data-service.ts,infra/fake/fake-fasta.ts,infra/reactor/{scanner,checker}.ts}`. NCBI: `fasta.cpp`, `line_reader.cpp`, `ncbistre.cpp`, `blast_fasta_input.cpp`. 32 rows in AD.tsv. All oracle runs used `scratch_AD/.ncbirc` (`DATA_LOADERS=none`), outfmt 6 `qseqid sseqid qstart qend length qlen`, subject `sub.fa` (one 240 nt record). No build, no /mnt/c. Files: `scratch_AD/run.sh` (NCBI), `runl.sh` (LOSAT-native), `mk.py` and inline Python for the inputs, `o_*` NCBI outputs, `l_*` LOSAT-before outputs (`.cmd .out .err .rc`).

## Call path (adapter)
- `register(program, role, bytes)` (store.rs:39) runs the program's checks, then bio parse, then `scan::scan` + `check_scan`, then the residue/defline/no-residue checks, stores `Vec<bio::Record>` under a handle and emits `{handle, records:[{index,id,length}]}` on stream 2. It never sees argv.
- `scan_begin(0)/scan_chunk/scan_end` (lib.rs:445-489): a `Scanner` per handle, errors are only reported by `scan_end`; the map entry is removed by `scan_end` either way.
- `run(argv, qh, sh)` (lib.rs:494, run.rs:203): `with_inputs` checks role and program of the handles, `parse(argv)` with the CLI parser, one `run_local_*` call writes outfmt 0/6/7 to stream writers, diagnostics to stream 3, HSP records to stream 1. Diagnostics are sent only on success.
- The app builds the run input from record ranges of the sources (data-service.ts runInputBytes), calls register, compares `{id,length}` with its scan table (`recordMismatch`).

## State that persists
Scan today: UTF-8 partial char, pending whitespace run, header bytes, record builder. NCBI reader state that a new scan must carry: EOL style (unknown, crlf, lf, cr, mixed), 1 byte of lookahead (CR vs CRLF), the "dropped LF" of the pushback quirk, line number, need_defline, "record has residues yet" (CheckDataLine), the 70-byte window of the current truncated line, ';' seen, mask state, title first word and last 50 bytes. All are functions of the byte stream, so chunked reproduction is possible; no look-behind beyond one byte.

## Oracle observations (NCBI 2.17.0 blastn, all exit 0 unless noted)
- NBSP at the end of a data line (`a1_nbsp`, `a1b_nbsp_first`): stderr `FASTA-Reader: Ignoring invalid residues at position(s): On line 3: 141-142`, query length 240 (the two bytes removed; also on the first data line, `On line 2: 101-102`, no error). LOSAT-before: error `query record 1 has a non-ASCII byte in a sequence line ...`.
- `>?100` / `>?unk100` (`a2_gap100`, `a2_gapunk`): 100 N inserted (qlen 340, two HSPs). `>?` alone (`a2_gapq`): stderr `CFastaReader: Bad gap size at line 3`, qlen 241. LOSAT-before: `has a defline that starts with '?'`.
- `>?_x desc` (`b9`): a defline with ID `x`. LOSAT-before rejects it as a '?' defline.
- Empty defline (`a3_emptydef`, `b13`): IDs `Query_1`, `Query_2` by position. Empty record: `Warning: [blastn] Query_2 q2: Sequence contains no data ` (trailing space), the other records are searched (`a6_empty_first`, `b12`).
- EOL: LF/CRLF/CR-only/mixed files all give the same 240-residue record (`a4_*`); LOSAT-before accepts CRLF and the two mixed files, rejects CR-only ("defline has the control character 0x0d"). Line numbers in warnings count `\r\n`, `\r`, `\n` as one line each and count comment/blank lines (`c1`-`c3`, `c5`: `On line 9`).
- Lost-LF quirk (`d1_lfcr_defline`): LF file with `...\n\r>q2 x\nACGT...`: stderr title warning and `Warning: [blastn] Query_2 q2 xCGTCAATACGGTT.. : Sequence contains no data`: the defline absorbed the next sequence line. In `d2` (CRLF first) and `d3` (mixed first) q2 is read correctly (length 140). `c4` (LF then `\n\r` inside sequence): reported `On line 4: 71`, the LF was lost. Source: `x_AdvanceEOLSimple` pushes back the rest of the line without the consumed delimiter (line_reader.cpp:245-268). BLASTX's `FastaStream` reproduces this; kind 0 cannot.
- `  >q2` with leading spaces (`b5`) is a data line (warning `On line 3: 1-3`); `>q2 x` at a line start ends the record (`b6`).
- CheckDataLine (`b1`, `b3`, exit 1): `BLAST query error: CFastaReader: Near line 2, there's a line that doesn't look like plausible data, but it's not marked as defline or comment.` A digits-only first line fails at once (`b3`). Not on later lines of a record with residues (`b4`: warning only). `b2` (100 N first line) prints nothing: the >40% ambiguity warning is ignored by `IgnoreProblem(eProblem_TooManyAmbiguousResidues)` (blast_fasta_input.cpp:358-360). `b7`: `CFastaReader: Hyphens are invalid and will be ignored around line 2`; `b8`: `;` comment to end of line is silent; `b11`/`b12`: trailing blank/comment lines add no record.
- Text before the first defline (`a5_blank`, `a5_comment`, `a5_nodef`): blank and comment lines skipped; sequence without defline is a record named `Query_1`. LOSAT-before: `FASTA that bio cannot read ... Expected > at record start.`
- Titles over 1000 bytes only warn (TooLong ignored): no effect on records.

## ABI v1 (TD-1) summary
BLASTN v1: no accepted input found where bio and NCBI differ (checks listed in row 23). TBLASTX v1 (outfmt 6) and BLASTP v1 run no defline or sequence-line byte check, so inputs accepted today with empty/tab/non-ASCII/leading-space deflines or `>?` lines would be read differently by the new reader (IDs, gap lines); an NBSP at a line end gives the same residues plus a warning. TBLASTN is not in v1; BLASTX v1 uses its own NCBI-style reader. The maintainer's rule (current checks first, then the new reader) makes these intended changes; list them in the V-ABI/v1 evidence.

## Surprises
1. HSP `q_idx` is the index in the searched query list: with `-query_loc` and a skipped record it differs from the registered index (§8 says registered index). `s_idx` is fine.
2. The kind numbering is per program in the app (`fastaParser`), but TBLASTN/BLASTX need one kind per role.
3. `checkRecordTable` forbids a record without header (`header_offset < sequence_offset`), and run-input concatenation would merge a headerless record into the previous record.
4. The lost-LF pushback quirk changes record structure; the scanner must emulate the exact line splitter, not a "\r\n or \n" rule.
5. CheckDataLine fires on every data line until the record has a residue and needs a 70-byte window.
6. A failing run drops stream 3, so warnings the CLI prints before an error are not visible in ABI v2; errors in a later query batch (CLI prints batch 1 first) fail the whole run.

## Decisions for the session
- Header-less record and text before the first defline: `header_offset`/`sequence_offset` convention, app prefixing `>` on concatenation, or reject.
- `residue_counts`: stored upper-case residues (proposed) vs original bytes; mask kept as case or intervals in the registered record.
- Message texts for the scan/register errors (CheckDataLine, `>?`), with or without the `BLAST query error: ` prefix; query reader errors rejected at register (proposed) although NCBI reports them when the batch is read.
- Whether to keep scan kind 0 after the move (the app would no longer use it).
- `>?_x` is a legitimate defline (port); other `>?` lines stay rejected (maintainer).
- q_idx mapping for skipped `-query_loc` records.
- UNSURE: BLASTN q_idx across query batches (global?) not verified; the HSP record of a skipped record cannot exist.

## Test inputs to reuse
`scratch_AD/*.fa` (a*, b*, c*, d*) and the 12180-row C++ oracle `LOSAT/tests/unit/blastx_stage_e_io_stream_expected.tsv` for the line splitter.
