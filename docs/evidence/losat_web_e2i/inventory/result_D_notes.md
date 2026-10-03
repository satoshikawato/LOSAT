# Range D (lookup table, subject scan, ungapped extension): notes

Final code: LOSAT commit 90c5f0181 (the working tree differs only by a `help` attribute on `-task` in `args.rs`; the later commit af03129ee adds fixtures only). Binary: `/home/kawato/.cache/losat-web-gui-target/sd/bin/LOSAT-90c5f0181`, NCBI 2.17.0 `/home/kawato/micromamba/bin/blastn`. Scratch: `/home/kawato/.cache/losat-web-gui-target/sd-res-D` (inputs `in/`, outputs `out/`, scripts `gen*.py`, `cmp.sh`, `matrix*.sh`, `cmp_idx.py`, `cmp_scan.py`, `mkresult.py`).

## Result

63 inventory rows: 39 ported, 16 faithful, 4 n/a, 4 rejected, 0 exception, 0 GAP. 4 extra rows (X1-X4). No UNSURE rows.

## GAP rows

None. Nothing on the path of range D was found that is neither ported nor rejected.

## What I compared

Code, line by line against NCBI (CRLF sources, pinned 598d8ae6):

- Index functions: the 12 `DiscontigIndex_*` of `blast_nalookup.h` against `discontig_index_*` in `disc_lookup.rs` by script (`cmp_idx.py`): every (mask, shift) term, in order, and the lo/hi split, 12 of 12 equal. `ComputeDiscontiguousIndex` dispatch: each case calls the same-named function; the default returns 0. Enum values 0..12 and 0..2 equal.
- `s_GetDiscTemplateType`: the Rust match equals the C if-chain (including TwoTemplates -> Coding value, `next()` = template + 1).
- `s_FillDiscMBTable` against `build_disc_mb_lookup`: 1-based index = `query_offset + pos + 2 - tl` (C: `from = left - (tl - 2)`, first word added when `seq` reaches `seq0 + tl`), ambiguity reset (`val & 0xFC`) also restarts the tl-base count (C: `pos = seq + tl`), `accum` is a u64 shifted by 2 and never masked, reset per location, no skip for locations shorter than tl, chain order newest (largest offset) first, `next_pos`/`next_pos2` sized concatenated length + 1, PV bit set only when the head cell was empty, and the PV array is shared by both templates (`pv_set_shift` for `ecode2` on the same array, `for_each_hit2` tests the same PV), `longest_chain` = (c1+1)+(c2+1) with the same 2048-cell helper arrays. PV size/shift arithmetic = `compute_mb_pv_params` (the contiguous table's), hashsize 4^11 / 4^12.
- Scans: `s_MB_DiscWordScanSubject_1`, `_TwoTemplates_1`, `_11_18_1`, `_11_21_1` against the four Rust functions: fill loop, the `switch (index - (sr0 + tl))` on 1/2/3 extra bases (`accum >>= 8; s--` for 3), the phases `accum >> 0/6/4/2`, the break checks before each phase; for `_11_18_1` and `_11_21_1` the four index expressions (7,7,7,8 and 8,8,8,9 terms) and the lo/hi rotation statements compared by script (`cmp_scan.py`). `s_MBChooseScanSubject` disc branch = `choose_disc_scan_subject`. NCBI's `max_hits` batching with the offset-array `break` is replaced by streaming (`on_word`) in the same order; batching cannot change the pairs.
- `BlastNaWordFinder`: word_length = lut_word_length = template_length, `scan_range[2] = length - tl`, masked subject keeps the disc scanner with `scan_range[1] = left`, `scan_range[2] = right - tl`, `s_DetermineScanningOffsets` unchanged, `s_range = end + tl`.
- Two-hit extension (`s_BlastnDiagTableExtendInitialHit`, `s_BlastnDiagHashExtendInitialHit`, `s_BlastDiagHashInsert`, `BlastExtendWordNew`, `Blast_ExtendWordExit`) against run.rs: the window is `config.window_size` (40 for dc, 0 otherwise) in the diag array length (`>= qlen + window`), the initial offset, the `hit_len_array` allocation (only if window > 0), `DiagHashTable::new(window)`, the insert window `window + Delta + 1 = 41`, `advance_diag_table_offset(.., window)` at both call sites, `two_hits`, `hit_saved || s_end_pos > last_hit + window`, Delta = min(0, 40 - tl) = 0 (no off-diagonal search), `s_TypeOfWord` early return 1 (`word_length == lut_word_length`, now both = tl) so word_type is always 1, first hit `hit_ready = 0` with `hit_len = s_end_pos - s_off_pos` (u8), `last_hit`/`flag` update after a saved hit from the ungapped end, ungapped extension `s_NuclUngappedExtend` with `s_match_end = s_off + tl` (word_length >= 16 so never the exact variant), cutoff test `off_found || score >= cutoff`.

Runs (stdout, stderr and exit status compared; all equal unless noted below), about 590 NCBI/LOSAT pairs:

- dc-megablast, all 18 combinations of `-template_type` (coding, optimal, coding_and_optimal) x `-template_length` (16, 18, 21) x `-word_size` (11, 12), with `-lcase_masking` and lowercase in the subject (the masked-subject scan is reached: `-lcase_masking` changes NCBI's output on these subjects): 72 runs on small hand-made sets (`q_small`/`q_large`, `q2_*`, `q3_*` x plain/soft-masked subjects with 1 to 22 base lowercase gaps, lowercase at both ends), 20 runs query LC738874 vs LC738875 with 400 random lowercase intervals (1 to 1000 bases, plus both ends, 26,229 lowercase bases of 359,647), 36 runs of 60 mutated queries (up to 10,222 lines of output) vs the same subject, 72 runs on dust/lowercase queries.
- Both diag containers: single queries of 3998..4001 bases (concatenated length on both sides of 8000), small sets (array) and large sets (hash, concatenated > 8000).
- Stale state across subjects (window 40): duplicate and near-duplicate subjects (`dup_s.fa`), 1 and 2 templates, template lengths 16 and 21, outfmt 6 and 0 (11,251 lines), and LOSAT `-num_threads 2` and `4` equal NCBI.
- Aliasing of the diag array: hairpin subjects (rc of the first tl bases of the query, a filler, then the same bases forward) for query lengths 1010..1025 and tl 16/18/21, built so that a diag array of 2^11 entries instead of 2^12 would alias two words of one subject into a false two-hit (concatenated query length in (2^11 - 40, 2^11]); NCBI and LOSAT report no HSP; a control with two exact words on one diagonal reports the same single HSP in both. This test cannot fail with NCBI's size (the `+ window` exists to prevent exactly this), so it only guards against a window-less size.
- Degenerate inputs: queries of 10..40 bases, poly-A / AT repeats (dust), all-N, N/R/Y inside a query, lowercase query islands; subjects of 10..45 bases (plain and lowercase), all-lowercase subject, subjects with N and IUPAC codes: identical.
- megablast with a template (`-word_size 11` and `12`, window 0, one-hit): 42 and 27,196 lines, identical. `-task blastn -template_*` is an NCBI options error ("Invalid lookup table type for discontiguous Mega BLAST"), same text and exit 1 in both.
- Errors: dc-megablast `-word_size 10`, `7`, `28`, `13` ("word size must be either 11 or 12"), blastn-short with template args (word-size message, and "Invalid lookup table type for discontiguous Mega BLAST" with `-word_size 11` or `12`): identical. `-template_length 17` and `-template_length` without `-template_type` give clap text and exit 2 (NCBI USAGE, exit 1): the approved exception for argument-parser errors.
- blastn-short: 72 runs (300, 15, 1 short queries from LC738875 with mutations/N/both strands, plain and soft-masked 60 kb subject, default, outfmt 0, `-word_size` 4,5,6,8,9,10,11,12,13,16) and 20 runs with subjects of 4..20 bases: identical.
- Rejected options (`-window_size`, `-off_diagonal_range`, `-xdrop_ungap`, `-use_index`, `-index_name`, `-xdrop_gap`, `-soft_masking`, `-ungapped`) with both tasks: "the NCBI BLAST+ option -X is not supported by LOSAT's BLASTN", exit 2. `-strand plus/minus` is also rejected (not a range D row).

Not run: the `INT4_MAX / 4` (about 536 Mb of subjects) offset wrap in `Blast_ExtendWordExit`/`s_BlastDiagClear` (row 48), compared by reading only; `cargo test` was not run (read-only task), the unit test names in the TSV are from reading the test modules.

## Inventory corrections

- None changes a status. Small points: `BlastChooseNucleotideScanSubjectAny` is at blast_nascan.c:2993-3005 (inventory 2994-3007); the inventory's line numbers for the LOSAT code are of commit a92fa902f (as stated).
- Row 1 (dc) wording "LOSAT keys on the task string" was true before SD; the final code keys on `mb_template_length > 0` (`run.rs:6595`), so megablast with `-template_*` builds a disc table too (compared, row 1 evidence).
- Inventory row 3 (CreateTask) is a cross-range row; the final `task_defaults` covers it.

## Extra rows (X1-X4, in result_D.tsv)

- X1 `BlastInitialWordParametersNew` container choice (blast_parameters.c:166-233): faithful, `run.rs:7698`.
- X2 `LookupTableWrapInit_MT` / `EstimateNumTableEntries`: n/a for dc (entries only size the PV).
- X3 `CSetupFactory::InitializeMegablastDbIndex` (`-use_index` with a template): rejected with `-use_index`.
- X4 `JumperNaWordFinder` disc branch (mapper only): n/a.
