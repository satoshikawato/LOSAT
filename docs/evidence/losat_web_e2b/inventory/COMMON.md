# S08 inventory: NCBI TBLASTX outfmt 0 / 7 path (common part)

You inventory, read-only, the NCBI C/C++ code that produces the bytes of `tblastx -outfmt 0` (pairwise report) and `tblastx -outfmt 7` (tabular with comments), plus the stderr text and exit status that go with them, and map each piece to what LOSAT (a Rust port of NCBI BLAST) already has. LOSAT's TBLASTX currently accepts only `-outfmt 6`; outfmt 0 and 7 are being ported in this session, so much of the display path is expected to be `missing` for TBLASTX, but many pieces exist for BLASTN, TBLASTN, BLASTP or BLASTX and may be reusable.

Do not modify, build or commit anything in the repositories. The only files you write are your range's TSV and notes file (paths in your range part), and scratch inputs/outputs under your own scratch directory (given in your range part).

## Sources
- NCBI C/C++ (the only authority): `/mnt/c/Users/genom/GitHub/ncbi-blast/c++` (pinned commit 598d8ae6). Do not run `git status` or `git log` there (very slow). Files have CRLF line endings: count line numbers with `tr -d '\r' < FILE | cat -n` (or `grep -n` works too). Paths below are relative to `c++/`.
- NCBI BLAST+ 2.17.0 binaries, comparison oracle only: `/home/kawato/micromamba/bin/tblastx` (also `blastn`, `tblastn`). You may run them on small inputs to confirm a behaviour (write inputs and outputs only under your scratch directory). Keep each run short (use small FASTA, a few hundred to a few thousand nt). Do not set BATCH_SIZE, CHUNK_SIZE, CTOOLKIT_COMPATIBLE, BL2SEQ_LEGACY, OLD_FSC or create a `.ncbirc` unless your range part asks you to study that variable, and then only in that one command's environment.
- Example inputs (read-only): `/mnt/c/Users/genom/GitHub/LOSAT-web-gui/LOSAT/tests/fasta/` (e.g. `LC738874.fasta` vs `LC738875.fasta` gives 4319 tblastx HSPs in ~6 s; `blastn_parity_compact.fasta`, `outfmt0/*.fasta` are small nucleotide sets). You can cut subsequences into your scratch directory.
- LOSAT: `/mnt/c/Users/genom/GitHub/LOSAT-web-gui/LOSAT/src`. Relevant parts:
  - `report/pairwise.rs` (outfmt 0 writers: BLASTP `write_blastp_pairwise_report`, BLASTN `write_blastn_pairwise_report`, `write_blastn_description_table`, TBLASTN `write_tblastn_pairwise_report`/`write_tblastn_hsp_info`/`write_tblastn_alignment`, BLASTX `write_blastx_pairwise_report`/`write_blastx_alignment`/`write_blastx_query_footer`/`write_blastx_epilog`, shared `write_subject_header`, `write_subject_summary_table_with_sum_n`, `coordinate_width`, `write_sequence_row`, `write_translated_pairwise_intro`, `write_blastp_database_header(_spacing)`, `write_blastp_query_header`, `write_no_hits_found`, `write_tblastn_unsearched_query_footer`, `write_ncbi_ka_field`, `is_positive_match`, `ncbi_percent_match`, `PairwiseHit`).
  - `report/outfmt6.rs` (tabular rows, `write_outfmt7_header`, `write_outfmt7_grouped`, number formats `format_evalue_ncbi`, `format_bitscore_ncbi`, `format_evalue_ncbi_tabular`), `report/defline.rs` (`ncbi_nucleotide_title`, `has_html_character_reference`), `report/query_warnings.rs`.
  - `common.rs` (`Hit`, final ordering `write_output_ncbi_order_evalue_hsp_order_to_writer`, comparators).
  - TBLASTX engine: `algorithm/tblastx/` (`args.rs`, `blast_engine/run_impl.rs` ~3500 lines: grep it; the final hits are built near the `let out_hit = Hit {` line; `translation.rs`; `ncbi_cutoffs.rs`; `sum_stats_linking/`; `filtering/`), CLI parsing `blastinput/value_parsers.rs` (`tblastx_outfmt`), `cli.rs` (`NativeError`), `main.rs`.
  - BLASTN's outfmt 0/7 driver (a model for how a program hands data to the shared writers): `algorithm/blastn/pairwise.rs`, `algorithm/blastn/blast_engine/run.rs` (grep `write_blastn_pairwise_report`, `write_outfmt7`, `QueryWarnings`, `write_pairwise_prologs`).
- Earlier authority records (what was already established for BLASTN outfmt 0; reuse their findings, re-verify only what is TBLASTX-specific): `/mnt/c/Users/genom/GitHub/LOSAT-web-gui/docs/evidence/losat_web_e2a/AUTHORITY.md` (§A NCBI call path for blastn -outfmt 0, §B, §C reuse table, §G additions; §G.3 description-table rule, §G.5 pieces for TBLASTX), `/mnt/c/Users/genom/GitHub/LOSAT-web-gui/docs/evidence/losat_web_e2c/AUTHORITY.md` (§E, §M, §N: input reading and titles), `/mnt/c/Users/genom/GitHub/LOSAT-web-gui/docs/evidence/losat_web_e2g/README.md` (warnings order, prolog flush, write failure exit 6, CTOOLKIT_COMPATIBLE).
- Project rules: `/mnt/c/Users/genom/GitHub/LOSAT-web-gui/AGENTS.md`. Note the approved TBLASTX exception: for local `-subject` searches LOSAT applies a non-default `-db_gencode` to subject translation/search/reporting even where NCBI's local `-subject` does not. Differences caused solely by that are not defects; everything else must match NCBI.

## Path and option scope
NCBI `tblastx` 2.17.0 run as `tblastx -query Q -subject S` (never `-db`, never remote), with these options and any of their values: `-evalue -threshold -word_size -num_threads -out -query_gencode -db_gencode -max_target_seqs -seg -window_size -culling_limit -outfmt` (`-outfmt` 0, 6 or 7 only, without custom field lists). Inputs are nucleotide FASTA with any IUPAC letters and lowercase, one or many queries and subjects, including short sequences (fewer than 3 nt), all-N sequences and an empty query file. Environment variables and `.ncbirc` are unset. Out of scope: `-db`, `-html`, `-line_length`, `-num_descriptions`, `-num_alignments`, `-show_gis`, `-sorthits`, `-sorthsps`, `-parse_deflines`, `-lcase_masking`, `-query_loc`, `-subject_loc`, `-strand`, `-max_hsps`, XML/JSON/ASN/SAM/custom formats (LOSAT rejects them; mark a branch reached only through them `rejected` and cite the LOSAT check if you find it, otherwise `n/a` with the reason).

## What to produce
Enumerate every NCBI function (and branch) in your range that executes on that path and can change outfmt 0 or outfmt 7 bytes, stderr, or the exit status (include static helpers; skip pure allocation/freeing). Start from the entry points named in your range part and follow the calls. For a function whose behaviour branches on the program (tblastx vs others), on the format, on an option or on the input, write one row per branch that tblastx reaches (and a row for a non-tblastx branch only if LOSAT's existing shared code takes it for TBLASTX by mistake). For each row:
- Quote the decisive NCBI line(s) verbatim (short) in `ncbi_snippet`.
- Describe concretely what tblastx prints or does (`tblastx_behavior`): exact text, spacing, number format, which value is used (e.g. "the HSP's sum_n from the Seq-align score 'sum_n'"), conditions.
- Find LOSAT's counterpart and classify (`status`):
  - `reusable`: an existing LOSAT function does exactly this for another program and can be called for TBLASTX unchanged (cite it).
  - `needs-param`: an existing LOSAT function is close; say exactly what must change for TBLASTX (and whether the change would alter what it prints for the programs that use it now).
  - `missing`: no LOSAT code does this.
  - `faithful`: LOSAT's TBLASTX already does this (e.g. outfmt 6 rows).
  - `divergent`: LOSAT's TBLASTX (or the shared code TBLASTX already uses) does it differently; say what differs and an input that shows it.
  - `rejected`: LOSAT stops with an explicit error before the input reaches this branch (cite the check).
  - `n/a`: executes but cannot affect output/stderr/exit status for tblastx (give the reason).
- `impact` (how broad the inputs that depend on the row are): `high` = default options, common inputs; `medium` = common non-default option or input feature; `low` = rare inputs or extreme values; `none` for n/a.
- When an oracle run confirms a behaviour, give the command and the observed bytes in `notes` (short), and keep the input under your scratch directory.

## TSV format (tab-separated, one header line, then one line per row)
`range	ncbi_file	ncbi_lines	ncbi_function	branch	ncbi_snippet	tblastx_behavior	losat_location	status	impact	notes`
- `ncbi_file` relative to `c++/`; `ncbi_lines` as `start-end` (the decisive lines; CR-stripped numbering).
- `branch`: `-` if none, else a short condition.
- `losat_location`: `path:line(function)` relative to `LOSAT/src`, `;`-separated if several, `-` if none.
- No tabs or newlines inside a field; keep `ncbi_snippet` under ~200 characters (use ` / ` between quoted lines).

## Notes file
Also write a Markdown notes file (path in your range part) with: anything that does not fit a row (cross-function behaviour, the order of writes and flushes, open questions), the oracle commands you ran and what they showed, and a short list of "what the port must do" for your range, each item with the NCBI file:line it rests on.

## Working rules (important: usage limits)
- Append each row to the TSV as soon as you have classified it (open, append, close), so partial work survives if you stop. Create the file with the header line if it does not exist. If the file already has rows, continue after them without duplicating.
- Write the notes file incrementally too.
- Work function by function. If you cannot decide in reasonable time, write the row with notes starting `UNSURE:` explaining what to check.
- Never guess NCBI behaviour; read the NCBI source, and use the oracle to confirm.
- When done, reply with: the number of rows per status, the list of `missing`, `needs-param`, `divergent` and `UNSURE` rows (function, branch, impact, one line), and the 5 most important things the port must get right for your range.
