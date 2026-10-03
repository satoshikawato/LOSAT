# S08 inventory results: what the TBLASTX outfmt 0/7 port did with each NCBI row

Session S08 (E2b) ported TBLASTX `-outfmt 0` and `-outfmt 7`. Before the port, read-only
agents inventoried NCBI's path in 7 ranges (A to G): the TSVs
`/home/kawato/.cache/losat-web-gui-target/s08/inventory/<R>.tsv` and notes `<R>_notes.md`
(conventions: `COMMON.md` in the same directory). Each row's `status` was LOSAT's state
BEFORE the port. You now record, read-only, what the FINAL code does for each row of your
range, and you check that nothing on the path was left neither ported nor rejected.

Do not modify, build or commit anything in the repositories. Write only your result TSV and
notes (paths in your range line) and scratch files under your scratch directory.

## Sources
- NCBI C/C++ (the only authority): `/mnt/c/Users/genom/GitHub/ncbi-blast/c++` (pinned 598d8ae6; CRLF line endings; do not run git there).
- LOSAT final code (commit `0d533ba76`): `/mnt/c/Users/genom/GitHub/LOSAT-web-gui/LOSAT/src`. Main places of the port:
  - `algorithm/tblastx/blast_engine/run_impl.rs`: `run` (CLI order: threads, subject read, empty subject error, query read, `-out`, options, `Query is Empty!`), `run_local`, `search` (format checks, `check_report_titles`, `-window_size 0` rejection), `tblastx_query_batch_size`, `run_in_pool` (prologs, batch size 0 error, batches), `search_query_batch` (per-batch search, unsearched batch, query stats), `preliminary_subject_bases` (ncbi2na random resolution), `write_tblastx_outputs`, `write_pairwise_prologs`.
  - `algorithm/tblastx/report.rs` (new): display translation (`display_base`, `display_residue`, `displayed_rows`), identities/positives (`row_counts`, `set_displayed_identities`), SEG masks (`query_dna_masks`, `lowercase_query_row`), `pairwise_hits`, `fasta_defline`, outfmt 7 headers and `write_tabular`, `final_hit_order` (hit list).
  - `report/pairwise.rs`: TBLASTX section (`TblastxPairwiseQuery`, `TblastxPairwiseReport`, `write_tblastx_pairwise_prolog`, `write_tblastx_hsp_info`, `write_tblastx_alignment`, `write_tblastx_query_footer`, `write_tblastx_epilog`, `write_tblastx_pairwise_report`) and the shared `write_blastn_description_table` (new `show_sum_n`).
  - `algorithm/tblastx/sum_stats_linking/linking.rs` (`hsp_link_num`, `num`), `chaining.rs` (`UngappedHit` fields), `report/defline.rs` (titles, `ncbi_nucleotide_title_reads_past_end`), `cli.rs` (`ReportStream`, `NativeError`), `blastinput/value_parsers.rs` (`tblastx_outfmt`), `algorithm/tblastx/args.rs` (`max_target_seqs: Option`).
- LOSAT binary built from that code: `/home/kawato/.cache/losat-web-gui-target/s08/native/release/LOSAT` (`LOSAT tblastx -query Q -subject S -outfmt F ...`). NCBI 2.17.0 binaries (comparison only): `/home/kawato/micromamba/bin/tblastx`. Keep runs small; write inputs/outputs only in your scratch directory; do not set BATCH_SIZE, CTOOLKIT_COMPATIBLE etc. except in one command's environment when a row needs it.
- Frozen fixtures (what is already verified byte for byte against NCBI): `LOSAT/tests/outfmt0_manifest.tsv` (`tblastx.*` rows, column `covers`), `LOSAT/tests/fixtures/tblastx_regression/manifest.tsv` and the case comments in `LOSAT/tests/tblastx_regression_fixtures.py`.

## Known decisions (classify accordingly, do not re-decide)
- Approved exceptions (AGENTS.md): a non-default `-db_gencode` is applied to the subject (NCBI local `-subject` uses code 1); argument-parser errors use clap text and exit 2 (NCBI USAGE, exit 1); with `-subject` LOSAT honors `-num_threads` and prints no thread warning; an outfmt 6/7 write failure exits non-zero (NCBI aborts); an outfmt 0 write to a closed pipe exits 6.
- Explicit rejections in S08 (LOSAT stops with a message; not a gap): `-window_size 0`; outfmt 0/7 deflines that are empty, non-ASCII, have a control character or an empty id; outfmt 0 subject deflines with an HTML character reference or that NCBI's `x_CleanAndCompress` reads past; a `BATCH_SIZE` that is not an integer; outfmt specs other than 0, 6, 7 (also `6 qseqid`, `7 std`, `+6`, `06`).
- Deferred to S08+ (next session, recorded in the handoff; classify `deferred`): NCBI's argument parsing of `-threshold`/`-evalue` (`ncbi_double`), integers (`ncbi_integer`), `-outfmt` text; the `-subject`-missing error; NCBI's file error texts for `-query`/`-subject`; records without data ("Sequence contains no data"); `check_ncbi_application_settings` (DIAG_*, NCBI_CONFIG_*, `.ncbirc` keys) for TBLASTX; `Subject_N` ids; options LOSAT's TBLASTX does not accept (`-num_descriptions`, `-num_alignments`, `-sorthits`, `-sorthsps`, `-sum_stats`, `-line_length`, `-html`, ...).

## What to produce
For EVERY row of your range's TSV (in order), one result row. Columns (tab-separated, header first):
`range	row	ncbi_file	ncbi_lines	ncbi_function	branch	status_before	s08_result	losat_location	evidence	notes`
- `row`: the 1-based row number in your range's TSV (header excluded).
- `ncbi_file`, `ncbi_lines`, `ncbi_function`, `branch`, `status_before`: copied from the inventory row.
- `s08_result`, one of:
  - `ported`: the final code now does what NCBI does (it was missing/needs-param/divergent).
  - `reused`: a shared LOSAT function does it for TBLASTX unchanged.
  - `faithful`: TBLASTX already did it before S08 and still does.
  - `rejected`: LOSAT stops with an explicit message before this branch (cite the check and the message start).
  - `exception`: an approved exception above.
  - `deferred`: in the deferred list above.
  - `n/a`: no effect on outfmt 0/7 bytes, stderr or exit status (keep the inventory's reason).
  - `GAP`: the final code neither does it as NCBI nor rejects it, and it is not in the lists above. Give a concrete input that shows the difference (run NCBI and LOSAT on it if you can, and put the command and the differing bytes in `notes`).
- `losat_location`: `path:line(function)` relative to `LOSAT/src` in the FINAL code (`-` if none).
- `evidence`: a fixture id that covers the row (e.g. `tblastx.multi.0`, `env.batch700`), a unit test name, or an oracle command you ran (short), or `-`.
- No tabs or newlines inside a field.

Append each result row as soon as you have it (open, append, close), creating the file with the header if absent; continue after existing rows without duplicating. Write the notes file incrementally too: the GAP rows in detail (input, NCBI bytes, LOSAT bytes), rows you could not decide (start with `UNSURE:`), and anything the inventory got wrong.

Read the NCBI source for each row before deciding; read the LOSAT code; when a row is not covered by a fixture and you are unsure, run both binaries on a small input. Never guess.

When done, reply with: counts per `s08_result`, every GAP and UNSURE row (row number, function, one line), and anything in the inventory you found to be wrong.
