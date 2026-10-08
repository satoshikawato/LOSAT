# SF (E2h) port plan: instruction steps 3 (engine) and 4 (adapter)

Draft written 2026-10-08 by a read-only planner. Inputs: the SF instruction (decisions 1-5), plan §4.2 and §7 SF, TD-2, TD-8 and TD-12, the AD/RD/RP/BI tables and notes (LR is still open), and `$WT` at eba45f75 (line numbers are from that commit).

## 1. Record type and how it reaches each consumer
**Decision:** `blastinput::fasta_reader::FastaRecord` is the record of BLASTN, TBLASTX, TBLASTN and BLASTP everywhere: the CLI, `run_local`, `seq_range.rs` and ABI v2 `register`.
- `local_id` is `Query_N` or `Subject_N`, where N is the record's position and every record counts.
- `title: Vec<u8>` holds any bytes; it is the title after `x_ApplyMods`.
- `sequence` holds IUPAC letters. U becomes T for nucleotides, `>?` gaps become N or X, and letters are lower case exactly on the `x_OpenMask`/`x_CloseMask` intervals.
- `warnings` holds the reader's lines in NCBI order: line warnings first, then the title warning.

**Additions in S1:**
- `FastaRecord::cut(from,to)` for ranges; it keeps the ID and title and drops the warnings.
- `FastaRecord::from_bio(rec, n, prefix, protein)` bridge: title = id + ' ' + desc, plus today's title warning. It is used only where LOSAT's current checks guarantee that bio reads like NCBI: ABI v1 (Q1), the adapter until S10, and TBLASTN's query until S8.
- `read_subjects(source, range, warn)`, a port of `ReadSequencesToBlast` + `GetAllSeqs` + the range check of `x_FastaToSeqLoc` (blast_input_aux.cpp:221-247, blast_input.cpp:198-220, blast_fasta_input.cpp:376-466):
  - warnings are written as each record is read;
  - a `ReadError::Parse` becomes `BLAST query error: <msg>` with exit 1 (`NativeError`);
  - a range starting past a record's end becomes today's `options_error`.
- Queries:
  - `FastaInputSource::stream_is_empty()` ports IsIStreamEmpty: a regular file with only white space is empty, a pipe never is.
  - `read_queries` returns `QueryRecords{records, end}`.
  - Each program's shared `search` takes `(&[FastaRecord], QueryEnd)`.
  - `run_local` passes `QueryEnd::Input`, so its TD-2 signature only changes `fasta::Record` to `FastaRecord`; Q2 explains why that is enough.
- `ncbi_environment::check_ncbi_application_settings` returns `ApplicationSettings{data_loaders}`:
  - precedence: environment `NCBI_CONFIG__BLAST__DATA_LOADERS`, then `<prog>.ini`, then `.ncbirc` (ncbireg.cpp:1567-1583);
  - the `blastdb`/`genbank`/`none` substring rules follow blast_scope_src.cpp:76-94;
  - `main.rs` passes the result to each `run`, which builds `ReaderConfig`;
  - the DATA_LOADERS entries stay accepted (they now only choose a reader class that LOSAT reproduces or rejects).
- Temporary trait `InputRecord {seq, cut, title_bytes}` for `bio::Record` and `FastaRecord`. It lets `seq_range.rs`, `report/query_warnings.rs`, `coordination.rs`, `seed/na_word_finder.rs` and `blastn/lookup.rs` build for both until S8 deletes the bio impl.

**Consumers:**
- **Masks:** `collect_lowercase_masks(record.seq())` (coordination.rs:938) does not change. Programs without `-lcase_masking` ignore case as today. `with_u_as_t` is deleted (the reader does U→T, reader.rs:1648-1660).
- **seq_range.rs:** `cut`, `cut_subjects`, `subjects_read`, `cut_queries` and `QueryInput::whole` take FastaRecord. `ordinals` stays (it equals N-1). `check_no_empty_interval` goes in S4, because an empty interval is a record without data.
- **Titles as bytes** (`Arc<[u8]>`, `&[u8]` writers):
  - tabular ID = `shown_id()` (s_ReplaceLocalId, tabular.cpp:474-504);
  - sseqid/sacc go through GetSeqIdList: `lcl|Subject_` gives the first token of GenerateDefline, and a BLASTP subject with an empty title gives `unnamed` (RP 17-19);
  - `Query= ` is raw and wrapped by NStr::Wrap at 68 bytes; `# Query: ` is raw;
  - the description list and the heading use GenerateDefline on bytes (RP 30-39);
  - warnings are `Warning: [prog] Query_N|id[ title]: …` with the 35/25-byte cut (query_warnings.rs:95-120);
  - outfmt 0/6/7 and stderr are raw bytes; only JSON is lossy.
- **BLASTX freeze:** nothing under `algorithm/blastx/` changes.
  - These keep their fields and signatures: `PairwiseHit` (subject_title: Option<String>), `BlastpPairwiseQuery`, `BlastpPairwiseReport`, `BlastxPairwiseOptions`, `PairwiseConfig`, `write_blastx_*`, `outfmt6::format_*`, `common::Hit`, `cli::{inaccessible,try_parse_from}`.
  - The four programs pass titles separately (`&[Arc<[u8]>]` by q_idx/s_idx) instead of `PairwiseHit.subject_title`.
  - Shared helpers (`write_flatfile_wrapped`, `write_blastp_query_header`, …) take bytes; BLASTX's entry points in pairwise.rs pass `as_bytes()`.
- **Adapter:** `Registered.records: Vec<FastaRecord>`. The JSON `id` is the title up to its first ' ' (empty with no title), as lossy UTF-8. `q_idx`/`s_idx` count records without residues.

## 2. Ordered implementer steps (one agent each; Q = quick tier, F = focused checks)
**Q:**
- `cargo fmt --check`;
- `cargo build --release --locked --target-dir $BUILD_ROOT/sf-native`;
- filtered `cargo test --all-features --target-dir $BUILD_ROOT/sf-test`;
- `check_losat.py --threads 1 --programs P` (the outfmt 0 manifest rows).

**F:**
- `LOSAT/tests/{blastn,tblastx,range}_regression_fixtures.py check --jobs 4`;
- `ci_fast_regressions.py --programs P --jobs 4` (S02, Gate A and TLOSAN Stage G hashes);
- E2g `check_inputs.py`, diffed against a run of the S11 binary: only rows that LOSAT rejects today may change;
- replay of the inventory cases (`scratch_{RD,AD,BI,RP}` `.cmd` and `runboth.py`, with paths moved to `$BUILD_ROOT/sf-e2h/inventory`) against the stored NCBI outputs (`*.n.*`, `o_*`). This needs no new NCBI run.

Run the standard tier at the end of each program step, before a push. Steps that touch `web_api`, `web/adapter` or `run_local` add V-ABI quick, the v1 WASI matrix (`check_wasm_threading.py`) and `v1_requests.js`.

**S0 Reader prerequisites (no output change; the reader-tests agent is on it now).**
- Reader unit tests green.
- A test of `fasta_reader::stream` against `LOSAT/tests/unit/blastx_stage_e_io_stream_expected.tsv` (12180 rows).
- A property test that `from_bytes(b)` equals `from_file(temp file with b)` for records, warnings and errors (the in_avail and lost-LF quirk).
- Reader-only timing against bio: a 100 Mb genome in 80-column lines, the same as one line, and 100k × 300 nt queries. Add the bulk path from R1 if the reader is more than 2× slower than bio.

**S1 Shared plumbing (no output change).**
- Files: `fasta_reader/mod.rs`, `ncbi_environment.rs`, `main.rs` (settings passed to the four `run` functions), `seq_range.rs`, `report/query_warnings.rs`, `coordination.rs`, `seed/na_word_finder.rs`, `blastn/lookup.rs` (made generic over `InputRecord`).
- NCBI: as in §1, plus blast_app_util.cpp:845-875.
- Removes nothing.
- Checks: Q and F for all programs, all byte-identical; unit tests for `read_subjects` (warnings stop at a range past the end) and for DATA_LOADERS precedence.

**S2 Report byte layer (no output change on accepted inputs).**
- `report/defline.rs`: a byte GenerateDefline:
  - strip trailing `.,;~ `, then the TPA/MAG prefixes;
  - HtmlDecode with its entity table, CUtf8::GuessEncoding and the Windows-1252 conversion (`x_AppendChar`);
  - strip trailing `,;~ `, then `x_CleanAndCompress` (signed-char last byte, approved exception 2).
  - Also: the `Unknown` / `Sequence with id … skipped` path (RP 36, 41) and a GetSeqIdList helper.
- Byte forms of the wrap helpers in `report/pairwise.rs`, the outfmt 7 `# Query:` writer and the warning ID/title formatting.
- NCBI: create_defline.cpp:3431-3446, 3952-3960, 4050-4095; ncbistr.cpp HtmlDecode; ncbi utf8 GuessEncoding and CharToSymbol; showdefline.cpp; showalign.cpp; align_format_util.cpp:729-746; tabular.cpp:474-504.
- Checks:
  - unit tests from the RP oracle (`scratch_RP/out*`, u05-u07, v01, v09);
  - Q for all programs; F with `ci_fast --programs blastx` plus the BLASTX cargo tests, because shared helpers changed;
  - the 204 outfmt 0 fixtures unchanged.

**S3 BLASTN reading and reports.**
- `blastn/blast_engine/run.rs`:
  - `run` reads subjects with `read_subjects`;
  - `search_cli` reads the query with `FastaInputSource` and `stream_is_empty`. An empty pipe gives zero batches (outfmt 0 prolog and epilog, `# BLAST processed 0 queries`). A white-space pipe gives the prolog, then `Empty CBlastQueryVector` with exit 3;
  - `run_local` takes FastaRecord; IDs and titles are bytes (`write_batch_reports`, `fasta_defline`).
- `coordination.rs` switches record type.
- v1 `run_web_pair` keeps its checks, then goes through `from_bio`.
- Adapter `run.rs` converts the registered bio records with `from_bio` (temporary).
- NCBI: blastn_app.cpp:199-214, 277-282; blast_args.cpp:2525-2557; blast_fasta_input.cpp:316-503; fasta.cpp:312-1679.
- Removed from BLASTN:
  - `check_deflines*`, `check_sequence_lines*`, `check_residues*`, `unreadable_fasta`, `read_records`, `parse_fasta`, `bio_records_of`;
  - `with_u_as_t`, the HTML-decoded-title rejection (run.rs:5720), and `an empty query from a stream without a position`;
  - `read_blastn_fasta_records`, and the dead `coordination::read_sequences` / `read_fasta_records`.
- Checks: Q and F with the 183 BLASTN fixtures, the BLASTN range fixtures, `check_inputs.py`, the RD/AD/BI/RP BLASTN cases, and the v1 web_api tests.
- Expected changes:
  - the input kinds rejected today now give NCBI's bytes: deflines (empty, white space, control bytes, tab, CR, non-ASCII, non-UTF-8, `>?100`, `>?unk100`, `>?`, `>?_x`, 20 nt followed by a space), sequence lines (non-IUPAC, digits, `*`, `-`, `;`, inner and leading spaces, NBSP) with `FASTA-Reader:` and `CFastaReader:` warnings, text before the first defline, BOM, CR-only and mixed line ends, CheckDataLine errors with exit 1, entity and non-UTF-8 subject titles in outfmt 0;
  - on accepted input, only the outfmt 6/7 sseqid of `Subject_…` first tokens changes (RP v01, v09).
  - Records without residues and empty ranges stay rejected until S4.

**S4 BLASTN records without residues and batch-time errors.**
- An empty query:
  - `Warning: [blastn] Query_N title: Sequence contains no data ` before its report;
  - the outfmt 0 `Length=0` / `No hits found` body; outfmt 7 counts it;
  - a batch of only empty queries gives `BLAST engine error: Warning: … ` with exit 3.
- An empty subject:
  - its warning comes before the prolog; it counts in the database statistics and effective search space and never hits;
  - if every subject is empty: `The average subject length is too short` with exit 3, after the prolog.
- An empty range interval follows the same path.
- `QueryEnd` at batch time, with `BatchSizeMixer` sizes:
  - an error in batch k comes after batches before k are flushed, with no epilog;
  - an eEOF on an empty first batch gives `Empty CBlastQueryVector`.
- NCBI: blast_setup_cxx.cpp:485-650, 733-800; blast_input.cpp:134-177; blast_format.cpp:1450-1452; blast_message.c:216.
- Removes `check_records_have_residues*` and `check_no_empty_interval` for BLASTN.
- Checks: Q and F, with the BI cases `e*`, `bb1`, `b_q_b2_*`, `ord*`, the RD `e_*` cases and the range fixtures.
- Expected change: those rejected cases become `same` or `same-error`.

**S5 TBLASTX, both roles nucleotide** (split at empty records if the context runs short).
- `tblastx/blast_engine/run_impl.rs` (`run`, `search_cli`, `run_local`, the warnings in `run_in_pool`, `check_subjects_not_empty`, `check_shown_subject_titles`); `tblastx/report.rs` (`fasta_defline`, the outfmt 7 header) as bytes; the S4 behaviours; v1 through `from_bio`.
- Removes the `fasta_input` checks, `read_nucleotide_subjects` for TBLASTX, `with_u_as_t`, and the HtmlDecode and CleanAndCompress title rejections.
- Checks: Q and F with the 84 TBLASTX fixtures, the TBLASTX range fixtures, `ci_fast --programs tblastx` (Gate A and TLOSAN hashes), the `x_*` cases; Gate A itself comes with step 7.
- Expected changes: the same kinds as S3 and S4.

**S6 TBLASTN subject (nucleotide).**
- `tblastn/args.rs` (`run`, `search_cli`, `run_local` take FastaRecord for both roles; the query stays bio plus protein checks through `from_bio`) and `tblastn/stage_e_report.rs` (IDs and titles as bytes).
- Delete `read_nucleotide_subjects` and `NucleotideSubjects`.
- Checks: Q and F with `check_losat --programs tblastn`, the TBLASTN range fixtures, `ci_fast --programs tblastn`, the `t_both` cases.

**S7 BLASTP, both roles protein.**
- `blastp/blast_engine.rs`: `run`, `search_cli`, `run_local`, `run_resolved_with_records`, `fasta_id`, `fasta_defline` and the outfmt 7 header.
- An empty subject title gives `unnamed protein product` in the description list and heading and `unnamed` as sseqid.
- Empty queries and batch errors as in S4; empty subjects stay in the database (the warning exists already).
- v1 (`run_web_pair`, the `run_web_pair_records` handles) goes through `from_bio`.
- Removes the BLASTP uses of `check_protein_input_of`, `check_protein_sequence_lines_of`, `check_protein_residues_of`, `write_protein_title_warnings`, `protein_title_warning` and `check_records_have_residues_of`.
- Checks: Q and F with `check_losat --programs blastp`, the BLASTP range fixtures, `ci_fast --programs blastp`, the `pp_*` cases, and `v1_requests.js` in the standard tier.
- Expected changes: protein kinds are read as NCBI reads them, with warnings for digits, `-`, `.` and NBSP; empty deflines; headerless records.

**S8 TBLASTN query (protein) and cleanup.**
- The TBLASTN query is read with the protein reader.
- Delete the rejection functions of `blastn/input.rs`, the bio impl of `InputRecord` and the `fasta_input` aliases. Move the helpers that survive (`open_input`, `standard_input`, `check_utf8_file_name`, directory-as-empty) to `blastinput/`.
- Proof by grep: `fasta::Record` remains only in `web_api.rs` (v1), adapter scan kind 0 and `algorithm/blastx/`.
- Check `BLASTINPUT_GEN_DELTA_SEQ` again now that gap lines are ported (RD's 6 runs).
- Checks: the standard tier for all programs.

**S9 Adapter scan kinds 1 (nucleotide flags) and 2 (protein flags).**
- `web/adapter/src/scan.rs` gets `enum Scanner {Bio, Ncbi{protein}}`, a push-based port of the CStreamLineReader end-of-line rules (lost-LF quirk, a CR at a chunk end) and of the ReadOneSeq record structure:
  - no warnings;
  - `>?` lines are rejected;
  - the read-forward rules for `line_layout`: skip white space, `-`, everything after `;`, invalid residues and comment lines;
  - `residue_counts` counts the stored residues upper-cased.
- `lib.rs` `scan_begin` accepts 1 and 2; kind 0 stays.
- New property test `web/adapter/tests/scan_ncbi_properties.rs` against `FastaInputSource::from_bytes`, using the generator and edge cases of the instruction's step 4.
- Checks: adapter `cargo fmt` and `cargo test`.

**S10 Adapter `register`, store, `run`, and the docs.**
- `store.rs`:
  - `register` reads with the engine reader for the program and role (`ReaderConfig` with data_loaders on, Q3);
  - it rejects `>?` (decision 3), Seq-id lines, reader errors and a comment-only query (Q2);
  - `check_scan` uses kind 1 (BLASTN, TBLASTX, TBLASTN subject) or kind 2 (BLASTP, TBLASTN query);
  - the JSON `id` and indexes follow §1.
- `run.rs` passes FastaRecord (the bridge goes); `q_idx` is mapped through `ordinals` (Q8).
- `docs/web/abi_v2.md`: the `register`, `run` and `scan` rows, §8 and §9 (kind 1 is no longer "BLASTX, SX").
- `v_abi.js` lines 239-241: a BLASTP record without a defline is now read.
- Checks: adapter tests, the reactors, `check_build_identity.py`, V-ABI quick, the v1 matrix; V-ABI full in step 7.
- Expected changes: `register` accepts the inputs it rejects today; a record with an empty title gets `id: ""`.

**Explicit rejections after S10:**
- a first line that NCBI may read as a Seq-id while data loaders are on (`seq_id.rs`, decision 2);
- `-parse_deflines`, `-html` and the other options already rejected;
- `>?` lines in `register` and `scan` only;
- records longer than 2^31 letters;
- RP row 20 (a `Subject_` token with an undecodable title, NCBI exit 255), pending the PD-LOSAT-NCBI-DEFECTS batch;
- the v1 checks (TD-1), file names that are not UTF-8, and the `-lcase_masking` rejection of BLASTP and TBLASTX in `cli.rs`.

## 3. Risks
- **R1 Reader speed.**
  - `FastaStream` takes one byte at a time (`raw_take`, pushback checks, `consume(1)`).
  - `parse_data_line` pushes each residue and calls `close_mask`.
  - bio uses `read_line` (memchr) and `extend`.
  - Fix: a bulk path when no pushback is active (memchr2 on the BufReader buffer, then table-driven runs of residues with `extend_from_slice`), equal to the byte path by property test; measure in S0 and V-PERF.
  - Gain: the CLI no longer holds the file bytes as well as the records.
- **R2 Timing.**
  - Subjects are read during argument processing: warnings go to stderr at once and errors come before option errors and before any output.
  - Query warnings are precomputed per record and written before a batch's first report, so threads cannot reorder them.
  - The failing batch is never searched; earlier output is flushed first; the outfmt 0 prolog comes before the first query read; `BatchSizeMixer` decides batch ends one at a time.
  - ABI v2 drops stream 3 on a failed run, so parity holds only for runs that succeed.
- **R3 ABI v1 freeze.** Reading v1 input with the new reader changes accepted TBLASTX and BLASTP v1 inputs (tab, control or empty deflines, `>?` lines become gaps; AD rows 24-25). Decision 4 assumed no change, so see Q1.
- **R4 Byte-frozen outputs:** S02 capture (236), Gate A, TLOSAN Stage G, outfmt 0 fixtures 204, BLASTN 183, TBLASTX 84, range 35, `check_inputs.py` 300, the LOSATX gate (BLASTX shares the report helpers), V-ABI expectations. Only the `Subject_` token change is intended on accepted inputs.
- **R5 `from_bytes` versus `from_file`:** the in_avail difference could split ABI v2 from the CLI on the lost-LF quirk; S0's property test covers it.
- **R6 Network:** the Seq-id path cannot be byte-compared. NCBI runs must keep `DATA_LOADERS=none` in their directory.
- **R7 Size:** `run.rs` has 13.9k lines, hence the S3/S4 split; TBLASTX may need S5a/S5b.

## 4. Open questions, with my recommendation
- **Q1 v1:** keep bio plus today's checks and convert with `from_bio`; this is a strict reading of TD-1. Record it as a new judgment because decision 4's premise fails for TBLASTX and BLASTP v1. The other option is the literal decision 4 with the changes listed.
- **Q2 `register` and query endings:** reader errors and comment-only queries end every run in failure, so `register` fails early with NCBI's text. `run_local` then needs only `QueryEnd::Input`.
- **Q3 Web data loaders:** on, as in the CLI default with no `.ncbirc`. Seq-id lines are then rejected in Web the same way as in the CLI.
- **Q4** Port HtmlDecode and GuessEncoding in S2: decision 5 needs them for non-UTF-8 titles in outfmt 0 anyway. RP row 20 stays a rejection until the batch question is answered.
- **Q5** A record without a defline in kinds 1/2 gets `header_offset == sequence_offset`, documented in §9; the app session handles it.
- **Q6** `residue_counts` counts the stored residues upper-cased (after U→T), on both sides.
- **Q7 Seq-id breadth:** keep `seq_id.rs` as it is (a line of letters only is FASTA; any other alphanumeric first line is rejected). Port the table-free shapes later if wanted (BI rows 15, 16, 18, 19).
- **Q8** `q_idx` refers to the registered index, mapped through `ordinals`; this matches §8.
- **Q9** Start freezing the step 5 NCBI fixtures with a sonnet agent during S1-S2. Until then, the inventory replay is the focused oracle.
- **Q10** Gap-line modifiers stay ported (reader.rs `gap_type_info`), with no rejection.
- **Q11** BI's open question on `NCBI_CONFIG__BLAST__DATA_LOADERS`: the environment registry has the highest priority in `app->GetConfig()` (ncbireg.cpp:1567-1583), so treat the variable like a registry entry.

### Critical Files for Implementation
- /home/kawato/losat-work/.worktrees/web-gui/LOSAT/src/blastinput/fasta_reader/mod.rs
- /home/kawato/losat-work/.worktrees/web-gui/LOSAT/src/algorithm/blastn/blast_engine/run.rs
- /home/kawato/losat-work/.worktrees/web-gui/LOSAT/src/report/pairwise.rs (and report/defline.rs)
- /home/kawato/losat-work/.worktrees/web-gui/LOSAT/src/blastinput/seq_range.rs
- /home/kawato/losat-work/.worktrees/web-gui/web/adapter/src/store.rs (and scan.rs)