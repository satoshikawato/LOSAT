# E2i inventory: classification instructions (common part)

You classify, read-only, how faithfully LOSAT (a Rust port of NCBI BLAST) ports the NCBI functions that run when NCBI `blastn` uses `-task dc-megablast` (discontiguous megablast) or `-task blastn-short`. Do not modify, build or commit anything in the repository. The only file you write is your range's TSV (path given in your range part); scratch files go under the scratch directory given in your range part.

## Sources
- NCBI C/C++ (the only authority): /mnt/c/Users/genom/GitHub/ncbi-blast/c++ (pinned commit 598d8ae6). Do not run `git status` there (very slow).
- LOSAT: /mnt/c/Users/genom/GitHub/LOSAT-web-gui/LOSAT/src (BLASTN engine: algorithm/blastn/; `blast_engine/run.rs` has ~13,400 lines, so grep it instead of reading it whole; the task defaults are in `algorithm/blastn/coordination.rs` (`task_defaults`, `configure_task`, `choose_na_lookup_table`, `finalize_task_config`); CLI parsing is `cli.rs` and `blastinput/value_parsers.rs` (`blastn_task`); arguments `algorithm/blastn/args.rs`).
- The E2g inventory of the shared BLASTN path for the tasks `blastn` and `megablast`: /mnt/c/Users/genom/GitHub/LOSAT-web-gui/docs/evidence/losat_web_e2g/INVENTORY.tsv (columns id, range, ncbi_file, ncbi_lines, ncbi_function, branch, losat_location, status, impact, notes, e2g_action, e2g_result). Rows marked faithful there were checked for `blastn` and `megablast`; for such a row you only need to check whether dc-megablast or blastn-short takes a different branch or passes different values (template length instead of word length, window size 40, a different program enum, different defaults) and whether LOSAT handles that.
- Earlier authority notes: docs/evidence/losat_web_e2c/AUTHORITY.md (§G lists LOSAT's explicit rejections), docs/evidence/losat_web_e2f/AUTHORITY.md (query batches and splitting), docs/evidence/losat_web_e2g/stage2/REVIEW.md.
- NCBI BLAST+ 2.17.0 executables (comparison only, you may run them on small inputs): /home/kawato/micromamba/bin (`blastn`). Test inputs: /mnt/c/Users/genom/GitHub/LOSAT-web-gui/LOSAT/tests/fasta (e.g. LC738874.fasta and LC738875.fasta, ~300 kb each). Run them only from your scratch directory.

## Current LOSAT state (important)
LOSAT's CLI rejects both tasks today (`blastinput/value_parsers.rs` `blastn_task`), and rejects `-template_type`, `-template_length`, `-window_size`, `-xdrop_ungap`, `-xdrop_gap`, `-xdrop_gap_final`, `-off_diagonal_range`, `-no_greedy`, `-ungapped`, `-soft_masking`, `-use_index` and others as unported options (`cli.rs` `is_unported_blastn_arg`). Classify what LOSAT's code downstream of the CLI would do if the task string `dc-megablast` or `blastn-short` reached it (e.g. `coordination.rs` treats every task other than `megablast` like `blastn`; `run.rs` sets `discontig_template = args.task == "dc-megablast"`, which forces the megablast lookup table but there is no discontiguous table, scan or template). Options that LOSAT rejects today and that the session will keep rejecting (all of the list above except `-template_type` and `-template_length`) need one row each with status `rejected` only where the option is in your range.

## Path and option scope
NCBI `blastn` 2.17.0 run as bl2seq: `blastn -query Q -subject S` (never -db), with `-task dc-megablast` or `-task blastn-short`, and these options with any of their values: -template_type (coding, optimal, coding_and_optimal) -template_length (16, 18, 21) -word_size -num_threads -evalue -perc_identity -max_target_seqs -max_hsps -out -reward -penalty -gapopen -gapextend -dust -lcase_masking -subject_besthit -outfmt (0, 6, 7 only, without custom fields). Inputs are nucleotide FASTA with any IUPAC letters, lowercase, one or many queries and subjects, short queries (blastn-short is meant for queries under 50 bases). Environment variables BATCH_SIZE/CHUNK_SIZE/OVERLAP_CHUNK_SIZE and .ncbirc are unset. Everything else is out of scope except where the in-scope path still calls a function.

## What to produce
Enumerate every NCBI function in your range that executes on that path and whose behaviour for dc-megablast or blastn-short differs from what E2g classified for blastn/megablast, or that E2g did not list because blastn/megablast never reach it (follow the calls from the entry points named in your range part; include static helpers and macros with logic). A function that runs identically for all tasks (same branch, same values) needs no row unless E2g missed it. For a function whose behaviour branches on the task, the template type/length, the window size, the lookup table type, the word size or the input, write one row per branch. Prefix the branch with `dc:` (dc-megablast), `short:` (blastn-short) or `both:`. Then find LOSAT's counterpart, read both, and classify each row:

- `faithful`: LOSAT already handles this branch with the same logic and the same C details: integer widths and wraparound (Int2/Int4/Uint4/Uint8/size_t), shifts and masks, float-to-int conversions, rounding, pointer and offset arithmetic (off-by-one, 1-based lookup offsets, `get_data` positions), sort order and tie-breaks (glibc 2.39 `qsort` is a stable merge sort), call order and timing. A LOSAT implementation that differs only for speed or memory but provably gives the same output is `faithful` with note `equivalent: <how>`.
- `divergent`: LOSAT has code for it but some logic, value or C detail differs for this task. Say exactly what differs and which inputs show it.
- `unported`: on the path, no LOSAT counterpart for this branch (e.g. the discontiguous table), and it can change output, warnings, errors or exit status.
- `rejected`: LOSAT stops with an explicit "not supported by LOSAT" error before the inputs reach this branch, and the session keeps it rejected (the option list above). Cite the LOSAT check.
- `n/a`: executes on the path but cannot affect output, warnings, errors or exit status (allocation, freeing, pure bookkeeping, a speed filter whose result is re-checked). Give the reason.

Impact (how broad the inputs are that can show a difference; for faithful/n/a rows write `none`): `high` = the task's default options or common inputs; `medium` = a common non-default option or input feature; `low` = rare inputs or extreme values; `none`.

## TSV format (tab-separated, one header line, then one line per row)
`range	ncbi_file	ncbi_lines	ncbi_function	branch	losat_location	status	impact	e2g_row	notes`
- ncbi_file relative to c++/ (e.g. src/algo/blast/core/na_ungapped.c); ncbi_lines `start-end` of the function or branch at the pinned commit.
- branch: `dc: ...`, `short: ...` or `both: ...` (a short condition, e.g. `dc: two_templates`, `short: word_size 7, approx entries < 250 -> lut 6`).
- losat_location: `path:line(function)` relative to LOSAT/src, `;`-separated if several, `-` if none.
- e2g_row: the E2g INVENTORY.tsv id(s) of the same function, or `-`.
- notes: one line, concrete (what you compared; for divergent/unported the exact difference, the NCBI values for each task (e.g. window_size 40, template_length 18) and an input that shows it if you can tell). No tabs or newlines inside a field.

## Working rules (important: usage limits)
- Append each row to the TSV as soon as you have classified it (open, append, close), so partial work survives if you stop. Create the file with the header line if it does not exist. If the file already has rows, continue after them without duplicating.
- Work function by function. Do not spend long on one function: if you cannot decide in reasonable time, write the row with status `divergent` or `unported` and notes starting `UNSURE:` explaining what to check.
- Never guess NCBI behaviour; read the NCBI source. When you run the NCBI executable to confirm a value (for example the epilog), say so in the notes.
- When done, reply with: the number of rows per status, and the list of `divergent`, `unported` and `UNSURE` rows (function, branch, impact, one-line difference).
