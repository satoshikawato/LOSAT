# Audit: protein/nucleotide FASTA input handling (BLASTP, TBLASTN, TBLASTX) and web entry points

Harness: `$D/work/inputs/c.sh NAME PROG QUERY SUBJECT [opts]` runs NCBI 2.17.0 and `$D/LOSAT` with -outfmt 0, 6 and 7
(stdout, stderr, exit status, and `2>&1`) and stores outputs in `$D/work/inputs/r/NAME/{n,l}.FMT.{out,err,rc,both}`.
All inputs are under `$D/work/inputs/{in,le,ws,nt,pc,fz}`. About 2100 cases were run; no panics were seen.
"D" below is `/home/kawato/.cache/losat-web-gui-target/s08p/audit`.

## Summary
- High (5): IN-2, IN-4, IN-8, IN-9 (executed, CONFIRMED); IN-10 (web adapter, code reading, not executed)
- Medium (7): IN-1, IN-3, IN-5, IN-6, IN-7, IN-11 (CONFIRMED); IN-14 (SUSPECTED)
- Low (2): IN-12 (SUSPECTED), IN-13 (CONFIRMED)

## Findings

### IN-8 (high, CONFIRMED) BLASTP identity/mismatch wrong when O meets X (O->X applied to the subject too and to the identity count)
- LOSAT: `LOSAT/src/algorithm/blastp/encoding.rs:54-62` (`encode_protein_sequence` maps O to X for every sequence), used for queries and subjects at `blastp/blast_engine.rs:660` and `:706`; identities counted from the encoded residues (`blast_engine.rs:2777`).
- NCBI: `src/algo/blast/api/blast_setup_cxx.cpp:894-899` replaces O by X in the QUERY search buffer only. NCBI reports identity from the original residues: query O vs subject O is an identity, query O vs subject X (or query X vs subject O) is a mismatch.
- Evidence (min_qO.faa = real protein with O at position 101, min_sX.faa = same with X there):
  `blastp -query in/min_qO.faa -subject in/min_sX.faa -outfmt 6`
  NCBI `qO sX 99.760 417 1 0 ...`   LOSAT `qO sX 100.000 417 0 0 ...`
  outfmt 0: NCBI `Identities = 416/417 (99%), Positives = 416/417 (99%)`; LOSAT `Identities = 417/417 (100%), Positives = 416/417 (99%)` (identities above positives).
  qX vs sO gives the same difference. qO vs sO is identical (both 100%). With composition-based stats default (2) and `-seg yes`, `-matrix` default.
  Fuzz: 24 random query/subject sets with O/X/U/B/Z/J/* injected; 12 differ, all only in pident/mismatch columns; after replacing every O by K all 24 are identical, so O is the only cause.
  TBLASTN (query O vs subject translated X from NNN) is identical, so only the BLASTP path is wrong.

### IN-2 (high, CONFIRMED) BLASTP outfmt 0, all subjects empty: "Effective search space used" is the query length, NCBI prints 0
- LOSAT: `blastp/blast_engine.rs:4998-5009` (`SearchSpace::for_database_search` with `total_db_len == 0`); footer printed from `search_spaces[q_idx].effective_space`.
- NCBI: `src/algo/blast/core/blast_setup.c:729-732` (`if (db_length == 0 && !SearchSpaceSet) return 0;` leaves eff_searchsp and length adjustment at 0).
- Commands: `blastp -query $F/e2e_protein_query.faa -subject in/s_empty_only.faa -outfmt 0` (s_empty_only = `>emptyonly` with no residues), also in/s_empty_two.faa.
  NCBI: `Effective search space used: 0` (all 8 stat blocks), LOSAT: `Effective search space used: 480`, `386`, `408`, ... (the query lengths). Everything else identical. outfmt 6/7 identical.

### IN-4 (high, CONFIRMED) Subject file whose first record is a bare ">" (no ID, no description, no residues) is read as having no records: generic "Empty CBlastQueryVector", exit 3
- LOSAT: `blastp/blast_engine.rs:4337-4338` and `blastn/input.rs:778-779` (`read_nucleotide_subjects`) return `empty_subjects_error` when bio reads no records, before the deferred defline check. bio's reader stops at the first empty record, so every record after the bare `>` is dropped too.
- NCBI: subject set has the empty-defline record (warning "Subject sequence contains no data") and all following records are searched.
- Commands (all three programs): `printf '>\n' > in/gt_only.faa; blastp -query $F/e2e_protein_query.faa -subject in/gt_only.faa -outfmt 6`
  NCBI rc=0 stderr `Warning: [blastp] Subject_1 : Subject sequence contains no data`; LOSAT rc=3 `BLAST engine error: Empty CBlastQueryVector`.
  Also le/p_gtonly_first.fa (`>` then a real record): NCBI reports the real record's hits, LOSAT rc=3. Same for le/n_gtonly_first.fa and le/n_gtonly_file*.fa with tblastn and tblastx (tblastx with a `>`-only file: NCBI rc=3 "The average subject length is too short" after the warning; LOSAT rc=3 "Empty CBlastQueryVector").
  Also with an empty query: NCBI `Query is Empty!` rc=0, LOSAT rc=3 (eq_bp_gtfirst_empty).
  The same bare ">" in the middle, at the end, or in a query is rejected explicitly (OK). Fix idea: reject (explicit "not supported by LOSAT's X") when the subject bytes are non-blank and bio returned no records.

### IN-9 (high, CONFIRMED, edge) TBLASTX `-query -` with `-subject -`: LOSAT says "Query is Empty!", NCBI runs with no queries
- LOSAT: `tblastx/blast_engine/run_impl.rs:907` (`seekable` ignores that stdin was already read as the subject). BLASTP and TBLASTN have the `query == "-" && subject == "-"` clause (`blastp/blast_engine.rs:4447-4452`) and reject explicitly.
- NCBI: `src/app/blast/blast_app_util.cpp:856-860` (tellg < 0 on the consumed stream = "piped input", not empty).
- Command: `tblastx -query - -subject - -outfmt 6 < in/n_base.fna`  NCBI: no output, rc 0 (outfmt 0 prints the prolog and statistics for 0 queries); LOSAT: `Warning: [tblastx] Query is Empty!`, rc 0. Also `-subject -` with the query omitted (sstdin3_tblastx).

### IN-1 (medium, CONFIRMED) Message order for a query with O that is also invalid (Karlin-Altschul message)
- LOSAT: `report/query_warnings.rs:107-135`, `blastp/blast_engine.rs:6940-6948`, `tblastn/args.rs:965-973` put the O message first. The comment at `query_warnings.rs:~142` ("O messages come before its Karlin-Altschul message", citing blast_setup_cxx.cpp:608-621) is not what NCBI prints.
- NCBI order (observed): `Could not calculate ungapped Karlin-Altschul ... filtering options One or more O characters replaced by X ... positions 0, 1, 2 `.
- Command: `blastp -query in/q_O3.faa -subject in/s_base.faa -outfmt 6` (query `OOO`; also `X*40+O+X*9`, 30 lowercase o, multi-query in/q_allOmulti.faa, and TBLASTN with the same files)
  NCBI: `Warning: [blastp] Query_1 O3: Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options One or more O characters replaced by X for alignment score calculations at positions 0, 1, 2 `
  LOSAT: `Warning: [blastp] Query_1 O3: One or more O characters replaced by X for alignment score calculations at positions 0, 1, 2 Could not calculate ungapped Karlin-Altschul ... filtering options `
  Stdout identical; stderr differs (all outfmts). Affects queries that are all X/O/invalid after replacement.

### IN-3 (medium, CONFIRMED) Protein 50-letter title warning fires for a defline with trailing white space; NCBI does not warn
- LOSAT: `blastn/input.rs:225-243` (`write_protein_title_warnings` uses `id + " " + desc` as bio trimmed it). The nucleotide path rejects this case explicitly (`input.rs:103`, "ends with white space after 20 nucleotide letters"); the protein path has no such check.
- NCBI: `src/objtools/readers/fasta.cpp:1651-1673` checks the raw title, so a trailing space ends the letter run.
- Command: query `>abc <50 letters> ` (one trailing space; in/q_tt_sp1.faa; also two spaces, space+tab) `blastp -query in/q_tt_sp1.faa -subject in/s_base.faa -outfmt 6`
  NCBI: no stderr output; LOSAT: `FASTA-Reader: Title ends with at least 50 valid amino acid characters.  Was the sequence accidentally put in the title line?`. Same for BLASTP subject (in/s_ttl_trailing_space.faa) and TBLASTN query. NCBI's rule for trailing control characters is odd (a trailing TAB, VT, FF or CRLF still warns; ' \t' does not), so a faithful port is "reject any defline ending in white space after 50 letters".

### IN-5 (medium, CONFIRMED) Protein sequence line ending in (or made of) Unicode white space: NCBI warns, LOSAT is silent
- LOSAT: `blastn/input.rs:340-351` (`check_protein_input_of` only checks deflines and letters; bio's `trim_end` drops U+00A0, U+0085, U+2003, U+2028, U+3000 at the end of a line). The nucleotide path rejects them (`check_sequence_lines_of`, input.rs:140).
- NCBI: `fasta.cpp:966-979` eCharType_Bad: `FASTA-Reader: Ignoring invalid residues at position(s): On line 3: 358-359`.
- Command: `blastp -query ws/p_seqlast_nbsp.fa -subject in/s_base.faa -outfmt 6` (last sequence line ends with NBSP, UTF-8 C2 A0). NCBI stderr has the warning, LOSAT none; stdout and rc identical. Also ws/p_seqtrail_*, ws/p_lineonly_* for emsp, ideo, lsep, nbsp, nel; BLASTP query, BLASTP subject, TBLASTN query. (zwsp, VT, FF, bell etc. are handled: VT/FF identical, others rejected.)

### IN-6 (medium, CONFIRMED) Protein query title warning placement with `2>&1` in outfmt 0
- LOSAT: `blastp/blast_engine.rs:4480` and `tblastn/args.rs:721` write the query title warnings before the prolog.
- NCBI prints the prolog (through `N sequences; M total letters` and its blank line) first, then the query warning, then the rest.
- Command: `blastp -query in/q_ttl_multi.faa -subject in/s_base.faa -outfmt 0 2>&1` (queries a and c have 50-letter titles). NCBI: warnings at lines 22-23 (after "5 sequences; 2,467 total letters"); LOSAT: lines 1-2. stdout and stderr separately identical; outfmt 6/7 identical. TBLASTN protein queries identical situation (tnq_ttl_multi). Subject title warnings and nucleotide (TBLASTX query) title warnings are identical.

### IN-7 (medium, CONFIRMED) BLASTP with an empty query: subject warnings differ
- (a) LOSAT writes the "Subject sequence contains no data" warnings before `Query is Empty!` (`blastp/blast_engine.rs:4436` runs before the `is_blank(query)` check at 4462); NCBI returns at `blastp_app.cpp:211-214` before `InitializeSubject`, so only `Query is Empty!` appears.
  Command: `blastp -query le/p_ws.fa -subject in/s_empty_mid2.faa -outfmt 6`  NCBI: `Warning: [blastp] Query is Empty!`; LOSAT: three `Subject_N ...: Subject sequence contains no data` lines first.
- (b) The subject residue/defline checks are deferred (`subject_checks`, 4344) until after `Query is Empty!`, so NCBI's reader warnings are lost: `blastp -query le/p_ws.fa -subject pc/p_mid_31.fa` (subject has a digit) NCBI: `FASTA-Reader: Ignoring invalid residues at position(s): On line 3: 41` then `Query is Empty!`; LOSAT: only `Query is Empty!`. Same for NBSP at the end of a subject line (eq_bp_nbsp_empty).

### IN-10 (high by code reading, NOT EXECUTED) Web adapter `register` and `run_local` do not apply the CLI's input checks for BLASTP and TBLASTN
- `web/adapter/src/store.rs:44-95`: `nucleotide` is `Some` only for BLASTN and TBLASTX; for BLASTP/TBLASTN the bytes go straight to `bio::io::fasta`. Missing versus the CLI: `check_deflines_of` (tabs, control characters, empty deflines), `check_sequence_lines_of`, `is_blank` handling, `check_protein_residues_of`, and for TBLASTN subjects the nucleotide checks.
  Consequences (from code and the CLI results above): `;` or `#` comment lines inside a file and `>` with an empty ID are accepted (bio appends comment text to the sequence, so BLASTP/TBLASTN would search different residues than NCBI, silently); a whitespace-only (non-empty) file fails with the generic bio error "failed to read query FASTA: Expected > at record start." instead of NCBI's `Query is Empty!`; NBSP/tab deflines are accepted.
- `blastp/blast_engine.rs:4537-4573` (`run_local`, used by the adapter and web ABI v1) never calls `check_protein_residues_of`, `check_protein_input_of`; digits, `-`, `.` in protein records are accepted silently while the CLI rejects them (TBLASTN's shared `search()` does call `check_protein_residues_of` and `check_residues_of`, so only its deflines/sequence-line checks are missing).
- Not built or run (rule: build nothing), so status is by code reading; the CLI-side rejections are verified (see "Checked OK").

### IN-14 (medium, SUSPECTED, code reading) BLASTP `run_local` has no empty-query / empty-subject branch
- `blastp/blast_engine.rs:4537-4573` versus `tblastn/args.rs:775-796` (TBLASTN `run_local` returns `Empty subjects` and `Query is Empty!`) and `tblastx/blast_engine/run_impl.rs` (`check_subjects_not_empty`). A 0-byte query or subject registered through the adapter for BLASTP reaches the search with zero records; NCBI prints `Query is Empty!` or exits 3 with `Empty CBlastQueryVector`.

### IN-11 (medium, CONFIRMED by code reading) Adapter `validate` does not check `-num_threads` for BLASTP and TBLASTX
- `web/adapter/src/run.rs:84-96`: BLASTP calls `blastp::blast_engine::check_options` (`blast_engine.rs:4419`, `resolve()` + `validate_requested_blastp_support`, no `validate_threads`); TBLASTX `check_options` (`run_impl.rs:874`) is `check_ncbi_options` + `check_losat_limits`, no `validate_threads`. TBLASTN is fine (its `search_settings`, `tblastn/args.rs:490`, calls it).
- CLI: `blastp|tblastx -num_threads 100000 ...` -> `Error: requested 100000 threads exceeds Rayon maximum 65535, which is not supported by LOSAT`; also `LOSAT_WASI_THREAD_CAP` and serial wasm builds ("unsupported num_threads=N") are checked by `validate_threads` only at run. So `validate` accepts values that `run` rejects (a validate/run mismatch for BLASTP and TBLASTX). Other option values (word size, threshold, evalue, seg, matrix, gap costs, comp_based_stats, window, db_gencode, culling, TBLASTN's max_intron_length etc.) go through the same `check_options`/`search_settings`/`check_ncbi_options`+`check_losat_limits` functions as the CLI, plus clap, so are consistent (CLI rejections were enumerated for tblastx and blastp; all come from clap or those functions).

### IN-12 (low, SUSPECTED, code reading) Web ABI v1 `parse_blastp_args` bypasses clap value constraints
- `LOSAT/src/web_api.rs:549-640`: `-max_target_seqs` and `-max_hsps` are parsed as plain numbers, so `0` is accepted; the CLI rejects both (`invalid value '0' ... expected an integer >= 1`, clap `blastn_count`/value parser). `-word_size 0/1` is still caught later by `validate_requested_blastp_support`. v1 has no TBLASTN. Not run.

### IN-13 (low, CONFIRMED) Reference/comment accuracy in session-added comments
- `report/query_warnings.rs:~142`: comment claims O messages come before the Karlin-Altschul message; NCBI shows the opposite (IN-1).
- `blastn/input.rs:210-224`: cites `fasta.cpp:1650-1673`; the amino-acid block is 1651-1673 (line 1650 is blank). `fasta.cpp:966-979` starts one line early (the quoted `HyphenToIgnoreAndWarn` case is at 967). `fasta.cpp:375-384` omits the `TruncateSpaces_Unsafe` line (376) that gives "lines containing only whitespace". These are off by one and do not change meaning.

## Checked OK (identical stdout, stderr, exit status and `2>&1` for outfmt 0, 6, 7 unless noted)
- Protein query residues, BLASTP and TBLASTN: single O, 5 O, exactly 20 O, 21 O, 25 O (",... (only first 20 shown)"), lowercase o, O first position, O in some queries only (q_mixO), U, B, Z, J, '*' (middle, start, end), lowercase and mixed-case sequences, all-X, all-U, all-'*', all-O (except the ordering in IN-1), X+O invalid-query combo (ordering IN-1 only). Every ASCII letter and '*' at start/middle/end of a sequence (pc/): identical for BLASTP query, BLASTP subject, TBLASTN query. Digits and all other punctuation, space and tab inside a sequence line, leading space: LOSAT explicit rejection "not supported by LOSAT's BLASTP/TBLASTN". Trailing space/tab on a sequence line: identical.
- O in the BLASTP subject (3 of N, all K replaced by O): identical; `-comp_based_stats 2`, `-seg yes`, `-matrix` variants identical (except IN-8 identity). `-comp_based_stats 0/1/3` explicit rejection.
- Records without residues: query empty (first, middle, last, only): LOSAT explicit rejection ("not supported by LOSAT's BLASTP/TBLASTN"; NCBI warns / exit 3); empty BLASTP subjects at any position (first/middle/last/only/two, with blank line, whitespace line, no final newline, long title, pipe IDs, `gi|`, `lcl|`): identical warning text and placement (Subject_N numbering, ID + title), except all-empty outfmt 0 (IN-2) and the bare `>` record (IN-4); `>   leading` empty-ID rejected explicitly. Empty TBLASTN/TBLASTX nucleotide subjects: LOSAT explicit rejection (`subject record N (...) has no residues ... not supported by LOSAT's TBLASTN/TBLASTX`).
- Empty/whitespace-only/newline-only subject file for all three programs: `BLAST engine error: Empty CBlastQueryVector`, exit 3, identical. Blank/empty/whitespace query file or stdin (redirected): `Warning: [prog] Query is Empty!` identical for all three programs; missing query/subject files and a directory as query/subject identical.
- 50-letter title warning, BLASTP/TBLASTN protein query and BLASTP subject: 49/50/51-letter boundaries, ID-only 51/50, lowercase, mixed case, X-only, `*`-only, ID with pipes, `gi|`, digit within the last 50 (no warning), non-letter at the 51st, CRLF, double space, trailing TAB; nucleotide inputs (TBLASTN subject, TBLASTX query and subject): 19/20/21 trailing ACGT, N, U, lowercase, IUPAC, protein letters (no warning), dash inside: identical; trailing-space nucleotide deflines are rejected explicitly. Multi-record title warnings and subject-warning order (titles, empty subjects, query titles, O messages) identical on stderr.
- Line endings and layout (BLASTP q/s, TBLASTN q/s, TBLASTX q/s): LF, CRLF, CRLF without final newline, `\r\r\n`, mixed LF/CRLF, no final newline, one-line sequences (100 kb), 5000-wide lines, width 1, 1.5 Mb single-line nucleotide subject, 200 kb protein subject, 20000-word 100 kb defline (query and subject, all three programs), blank line between/inside records, trailing blank lines/spaces, trailing whitespace on sequence lines, lowercase files: identical. Lone CR line endings, leading blank line, leading spaces line, `;` or `#` comment lines (first or middle), tab in defline, NUL in sequence, UTF-8 BOM, invalid UTF-8, non-ASCII bytes, `>` followed by spaces, bare `>` in the middle/last/no newline, missing defline: explicit "not supported by LOSAT's X" rejections (every format).
- Control/odd characters in deflines and sequence lines (bell, esc, FS, US, DEL, VT, FF, zwsp, NEL, U+2003/2028/3000, NBSP, é): identical or explicit rejection except the Unicode white space cases in IN-5.
- Nucleotide subject/query (TBLASTN subject, TBLASTX both): U/u, all IUPAC ambiguity codes (RYKMSWBDHVN, lowercase n), lowercase, mixed case, all-N, all-U, U-for-T, runs of N, long IUPAC mix: identical. `-`, `*`, digits, `.`, X, E, F, J, O, Q, Z in nucleotide records: explicit rejection.
- Query from stdin: `-query -` and omitted `-query`, valid/empty/blank (redirected file) for BLASTP, TBLASTN, TBLASTX: identical. Through a pipe: valid, O query, title warning identical; empty and blank piped query: LOSAT explicit rejection ("an empty query from a stream without a position ... not supported by LOSAT's X"; NCBI exits 0 or 3). `-query - -subject -`: BLASTP/TBLASTN explicit rejection; TBLASTX differs (IN-9). `-subject -` with a query file for all three programs: identical.
- Repository fixtures: all e2e_* query x subject pairs for BLASTP and TBLASTN, and the edge_*, iupac, rna, lcase_minus, mask, multi, tblastx_ambig pairs for TBLASTX and TBLASTN: identical. TBLASTN `-lcase_masking`, `-soft_masking true/false` with mixed-case and O queries identical.
- Web layer: `web_api.rs` session changes (struct field changes `subject: Some(subject)`, `threshold: 13.0`, `seg: "yes"`, `comp_based_stats: "2"`, `Option` max_target_seqs, word/window clamped to i32) are consistent with the CLI structs. v1 supports BLASTX, BLASTP, BLASTN, TBLASTX only (no TBLASTN).
- Options through adapter `validate`: BLASTP, TBLASTN, TBLASTX share the CLI's clap parser and `check_options` functions; `-out/-outfmt` rejected up front.

## NCBI references spot-checked (all read against /home/kawato/.cache/losat-web-gui-target/s08p/ncbi/c++)
1. fasta.cpp:1650-1673 (1651-1673 correct; off by one at start)  2. blast_setup_cxx.cpp:773-788 OK  3. fasta.cpp:966-979 (off by one at start)
4. blast_args.cpp:2553-2557 OK  5. fasta.cpp:375-384 OK (misses line 376)  6. blast_setup_cxx.cpp:894-899 OK  7. blast_setup_cxx.cpp:920-932 OK
8. blastp_app.cpp:211-214 OK  9. blast_app_util.cpp:856-860 OK  10. blast_format.cpp:1450-1452 OK  11. blast_args.cpp:3624-3627 OK
12. tblastn_app.cpp:213-216 OK  13. blast_app_util.hpp:252-255 OK  14. objmgr_query_data.cpp:375-380 OK  15. blast_args.cpp:3456-3481 OK (starts 3457)
16. seqsrc_multiseq.cpp:226-240 OK  17. blast_seqsrc.h:205 OK  18. local_db_adapter.cpp:133-135 OK  19. local_blast.cpp:177-180 OK  20. blast_aux.cpp:936-951 OK
21. blast_setup_cxx.cpp:608-621 (actual 608-622) OK in substance, but the claim built on it (O message before KA message) is wrong (IN-1).
