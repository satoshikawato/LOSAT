# E2d (session S11) inventory: `-query_loc` / `-subject_loc` (common part)

You trace, read-only, what NCBI BLAST+ 2.17.0 does with `-query_loc` and `-subject_loc` in `blastn`, `blastp`, `tblastn` and `tblastx` (bl2seq only: `-query Q -subject S`), and where LOSAT (a Rust port of NCBI BLAST) must do the same. Today all four LOSAT programs reject both options as "not supported by LOSAT's <PROGRAM>". The session that reads your table will port every NCBI function on this path wholesale, so the table must be complete for your range and concrete enough to port from (values, coordinates, message texts, exit codes, call order). Do not modify, build or commit anything in the repository. You write only the files named in your range part, under `/home/kawato/.cache/losat-web-gui-target/s11/inventory/`; put scratch files (inputs, oracle outputs, helper scripts) under `/home/kawato/.cache/losat-web-gui-target/s11/inventory/scratch_<RANGE>/`.

## Sources
**Never read or write anything under `/mnt/c` (the Windows drive mount fails under load). Everything you need has a copy on the Linux filesystem.**
- NCBI C/C++ (the only authority): `/home/kawato/.cache/losat-web-gui-target/s08p/ncbi/c++` (`src/` and `include/` of the pinned commit 598d8ae6). Cite NCBI files relative to `c++/` with these line numbers.
- LOSAT source, frozen copy of the base commit `4fab73fdb` (branch `feature/losat-web-gui`): `/home/kawato/.cache/losat-web-gui-target/s11/base-src/LOSAT/src`. Cite LOSAT locations relative to `LOSAT/src` with line numbers of this copy. Adapter (ABI v2 `validate`/`run`): `.../base-src/web/adapter/src`, contract `.../base-src/docs/web/abi_v2.md`.
  - CLI: `cli.rs` (single-dash translation, `is_unported_{blastn,blastp,tblastn,tblastx}_arg` that reject the two options today, `NativeError`), `main.rs`, `blastinput/app.rs` (NCBI's app layer for BLASTP/TBLASTN/TBLASTX: `options_error`, `engine_error`, `check_options` order), `blastinput/value_parsers.rs`, `blastinput/query_batch.rs`.
  - BLASTN: `algorithm/blastn/args.rs`, `algorithm/blastn/blast_engine/run.rs` (`pub fn run`, `run_local`; NCBI's order of option checks, subject read, query read, batches), `algorithm/blastn/input.rs` (FASTA reading and checks), `algorithm/blastn/coordination.rs`, `algorithm/blastn/query_split.rs`, `algorithm/blastn/filtering/` (DUST, lowercase), `algorithm/blastn/hsp.rs` (`parse_blastn_output_format`), `report/`.
  - BLASTP: `algorithm/blastp/args.rs`, `algorithm/blastp/blast_engine.rs` (8,000+ lines; grep: `run`, `run_local`, tabular fields `blastp_tabular_field`), `algorithm/blastp/query_split.rs`, `algorithm/common/protein_query_split.rs`.
  - TBLASTN: `algorithm/tblastn/args.rs` (`TblastnArgs`, `validate`, `run`, `run_local`), `search_*.rs`, `stage_d_*.rs`, `stage_e_report.rs`.
  - TBLASTX: `algorithm/tblastx/args.rs`, `algorithm/tblastx/blast_engine/run_impl.rs` (`run`), `algorithm/tblastx/report.rs`, `translation.rs`.
  - BLASTX (`algorithm/blastx/`) is out of scope (session SX); do not classify it, but you may read it as a hint (it has query-split and remap comments).
- Earlier authority records (NCBI behaviours already traced for these programs): `docs/evidence/losat_web_e2e/AUTHORITY.md` (BLASTP/TBLASTN/TBLASTX app layer §A, argument parsing §B, option order §D, input §J, query split), `docs/evidence/losat_web_e2c/AUTHORITY.md` (BLASTN CLI and input), `docs/evidence/losat_web_e2f/` (BLASTN query batches and split), `docs/evidence/losat_web_e2a/AUTHORITY.md` and `docs/evidence/losat_web_e2b/AUTHORITY.md` (outfmt 0 for BLASTN/BLASTP/TBLASTN and TBLASTX). Paths relative to `.../s11/base-src`.
- Approved exceptions (cite, do not count as defects): `.../s11/base-src/AGENTS.md` (TBLASTX and TBLASTN local `-subject` non-default `-db_gencode`), `docs/product_decisions/PD-LOSAT-CLI-NONSEARCH-DIFFERENCES.md`, `docs/product_decisions/PD-LOSAT-NCBI-DEFECTS.md`.
- NCBI BLAST+ 2.17.0 executables (comparison oracle only; small inputs, only from your scratch directory): `/home/kawato/micromamba/bin/{blastn,blastp,tblastn,tblastx}`. LOSAT binary of the base commit (run it to observe current behaviour; do not build): `/home/kawato/.cache/losat-web-gui-target/s11/bin/LOSAT-before` (`LOSAT-before blastn -query Q -subject S ...`). Test inputs: `.../s11/base-src/LOSAT/tests/fasta` (protein `AvCLPV.faa`, `CoBV.faa`, ...; nucleotide `LC738874.fasta`, `LC738875.fasta` ~300 kb, ...). Cut small inputs into your scratch directory (random or real; a few hundred to a few thousand letters). Keep each oracle run under a few minutes. Unset `BLASTDB`, `BATCH_SIZE`, `CHUNK_SIZE`, `OVERLAP_CHUNK_SIZE` and run where no `.ncbirc` is present.

## Scope
- Programs: `blastn` (tasks megablast default, blastn, blastn-short, dc-megablast), `blastp` (blastp, blastp-fast, blastp-short), `tblastn` (tblastn, tblastn-fast), `tblastx`. bl2seq only (`-query`, `-subject`), never `-db`, never `-remote`.
- Options: `-query_loc`, `-subject_loc` (any string NCBI's parser accepts or rejects), combined with the options LOSAT already supports for each program (e.g. `-strand` is rejected by LOSAT's BLASTN today: record what NCBI does with `-strand plus|minus` together with `-query_loc`, but classify `-strand` itself as out of scope), `-lcase_masking`, `-dust`, `-seg`, `-soft_masking`, `-query_gencode`, `-db_gencode`, `-comp_based_stats`, `-max_target_seqs`, `-evalue`, `-outfmt` 0, 6, 7 (default fields and the custom field lists that LOSAT supports for the program), `-num_threads`.
- Inputs: one or many query records and one or many subject records (the range applies to every record that the input source reads); records shorter than the range (end past the end; start past the end); lowercase letters; ambiguous letters; empty records; standard input (`-`).

## What to produce
For your range, enumerate every NCBI function (and every branch of it that the range options select) that executes on the range path and can change output bytes, warnings, errors or exit status; follow the calls from the entry points named in your range part (include static helpers and macros with logic); find the LOSAT place that corresponds (where the same data is produced today, e.g. where the query sequence is encoded, where the subject length enters statistics, where output coordinates are written), and classify each row:

- `faithful`: LOSAT already does exactly this for ranged input (rare: e.g. a function that receives the already-cut sequence and needs no change). Say why no change is needed.
- `divergent`: LOSAT has code for the same step, but it would give a different result with a range (e.g. it uses the full record length where NCBI uses the range length). Say exactly what must change.
- `unported`: the range-specific NCBI logic has no LOSAT counterpart (e.g. parsing the range, cutting the sequence, shifting coordinates).
- `rejected`: LOSAT should keep an explicit not-supported error for this branch (only with a reason: an NCBI crash, a `-db`/`-remote`-only branch, out-of-scope option).
- `exception`: an approved exception covers the difference. Cite it.
- `n/a`: executes but cannot affect output, warnings, errors or exit status. Give the reason.

`impact`: `high` (any ranged search), `medium` (a common combination, e.g. range plus multiple records, range end past the record end), `low` (rare inputs or extreme values), `none` for faithful/n/a/exception.

`proposal`: `keep`, `port:S|M|L` (S under ~100 Rust lines, M a few hundred, L more), or `reject:<reason>`.

## TSV format (tab-separated, one header line, then one line per row)
`range	row	ncbi_file	ncbi_lines	ncbi_function	branch	losat_location	status_before	impact	proposal	evidence	notes`
- `row`: 1, 2, ... within your range.
- `ncbi_file` relative to `c++/`; `ncbi_lines` `start-end`.
- `branch`: prefix with the program(s): `N:` blastn, `P:` blastp, `T:` tblastn, `X:` tblastx (e.g. `NX: query range with strand both`, `T: subject range, minus frames`). One row per program when LOSAT's place differs per program.
- `losat_location`: `path:line(function)` relative to `LOSAT/src`, `;`-separated, `-` if none.
- `evidence`: scratch file names of oracle (and LOSAT-before) runs that show the behaviour, or `source` if read only.
- `notes`: one line, concrete: what NCBI does (values, formulas, message text, exit code), what LOSAT does today, and what must change. No tabs or newlines inside a field.

Also write `<RANGE>_notes.md`: the call path you followed (entry points, order), the coordinate conventions at each step (0-based/1-based, open/closed, which strand, nucleotide or amino-acid units), the oracle observations (command, what it shows), surprises, and anything the session must decide. Oracle observations matter as much as the source reading: confirm every non-obvious claim with a small NCBI run, and keep the input and output files in your scratch directory.

## Working rules (important: usage limits)
- Append each row to the TSV as soon as you have classified it (open, append, close), so partial work survives if you stop. Create the file with the header line if it does not exist. If it already has rows (you were restarted), continue after them without duplicating. Append to `<RANGE>_notes.md` as you go as well.
- Work function by function. Do not spend long on one function: if you cannot decide in reasonable time, write the row with notes starting `UNSURE:` saying what to check.
- Never guess NCBI behaviour; read the NCBI source. When you run an executable to confirm, say so in `evidence`.
- Do not kill processes you did not start; never use `pkill` or `killall`.
- When done, reply with: the number of rows per status, the list of `divergent`, `unported` and `UNSURE` rows (row, function, branch, impact, proposal, one-line difference), and the most important points of your notes (especially every place where NCBI uses the full record versus the range, and every error with its text, exit code and timing).
