# Range A (application and arguments): notes

Final code `90c5f0181`, binary `/home/kawato/.cache/losat-web-gui-target/sd/bin/LOSAT-90c5f0181`, NCBI 2.17.0 `/home/kawato/micromamba/bin/blastn`. Scratch: `/home/kawato/.cache/losat-web-gui-target/sd-res-A` (`cmp.sh` runs both binaries and compares stdout, stderr, exit status; `runs/<label>/` keeps the outputs; `sweep1..9.txt`, `rejected_run.txt` list the results).

## What was compared (both binaries, stdout + stderr + exit status)
About 2,100 runs on `LOSAT/tests/fasta/outfmt0/dc_*`, `short_*`, mask/lcase/edge inputs and the BLASTN regression inputs:
- tasks: all values, five task names with wrong case or spelling (parser errors);
- `-template_type`/`-template_length`: alone, both, invalid, with dc-megablast, megablast, blastn, blastn-short;
- `-word_size` 3..101 and non-integer forms x templates x scoring errors x tasks (504 + 240 runs);
- 35 reward/penalty/gap sets for dc-megablast and blastn-short; `-evalue`, `-perc_identity`, `-max_target_seqs`, `-max_hsps`, `-dust` (14 forms) x `-lcase_masking` x 6 input pairs x 4 tasks (672 runs), `-subject_besthit`, `-outfmt` 0/6/7, `-num_threads`, `-out`;
- 263 runs that combine two options errors (first-error order);
- dc-megablast word_size {default,11,12} x 9 template pairs x 6 option sets (162 runs, all with hits); blastn-short word sizes 4..28 x 7 scoring sets x 4 e-values (420 runs);
- a 2.33 Mb query (so that blastn-short splits into 1 Mb chunks and dc-megablast does not) vs a 30 kb subject, outfmt 0/6/7.

Everything is equal except: (a) parser errors (NCBI USAGE block, exit 1; LOSAT clap text, exit 2: approved exception), (b) `-num_threads 2` with `-subject` (approved exception), (c) the explicit rejections (`rmblastn`, the options of `cli.rs` `is_unported_blastn_arg`, `-reward 0`, `-outfmt 5/8`), (d) the GAP below, (e) `-evalue`/`-perc_identity` in hexadecimal (see OBSERVATION).

## GAP 1 (row 27 and X1): `-template_length` in hexadecimal is accepted
Input (any task; dc-megablast shown):
```
blastn -task dc-megablast -query dc_small_query.fasta -subject dc_subject.fasta -outfmt 6 -template_type coding -template_length 0x12
```
- NCBI: stdout empty, exit 1; stderr is the USAGE block, then
  `Error: Argument "template_length". Illegal value, expected Permissible values: '16' '18' '21' :  \`0x12'` and `Error:  (CArgException::eConstraint) Argument "template_length". Illegal value, ...`.
- LOSAT: exit 0, stderr empty, stdout 178 bytes (the same three HSP lines as `-template_length 18`).
- Same for `0x10`, `0X15`, `0x012` (16, 21 forms). `0x`, `0x0`, `0x1g` are errors in both (parser class). With blastn or blastn-short plus a hexadecimal length LOSAT reaches the lookup-table error (exit 1) where NCBI stops in the parser.
- Cause: NCBI's `CArgAllowIntegerSet::Verify` (blast_input_aux.hpp:239, `DEFINE_CARGALLOW_SET_CLASS(CArgAllowIntegerSet, int, NStr::StringToInt)`) converts the string with `NStr::StringToInt` in base 10 (no `0x` fallback); the `CArg_Integer` conversion before it accepts `0x`. LOSAT `value_parsers.rs:382 blastn_template_length` uses `ncbi_constrained_integer`, which accepts `0x` (right for the `>=`/`<=` constraints, whose `StringToDouble` takes hex: `-word_size 0x0b`, `-reward 0x2`, `-penalty 0x0` are equal in both). Fix: the template length parser must reject a `0x` prefix after the conversion succeeds (error 'Illegal value'). Impact low.
- Other spellings of `-template_length` are equal (error in both): `18.0`, ` 18`, `18 `, `1e1`, `-18`, `4294967314`, `99999999999999999999`, full-width digits, `0b10010`; accepted in both: `018`, `+18`, `000000000000000018`.

## OBSERVATION (not task specific, not counted as a range-A GAP)
`-evalue 0x10` and `-perc_identity 0x32` (hexadecimal doubles): NCBI accepts them (exit 0, 178 bytes), LOSAT exits 2 ('invalid value'). The doc comment of `ncbi_double` in `value_parsers.rs` says LOSAT rejects NCBI's hexadecimal forms on purpose. It applies to every task, so it was left to the caller to decide whether it is an accepted exception.

## UNSURE
None.

## Inventory corrections and remarks
- Row 63 and the help texts: `-help` is clap text for every option (approved exception), so the usage rows are `n/a`.
- Row 60 (`Window for multiple hits`): the inventory is right. One more fact on the same function (`CBlastFormat::PrintEpilog`, blast_format.cpp:2272): `m_Program` is `Blast_ProgramNameFromType(...)` = `blastn` for all four tasks (megablast too), so NCBI's zero-gap formula `(m_Program == "megablast" || m_Program == "blastn") && GetGapExtensionCost() == 0` applies to dc-megablast and blastn-short as well. LOSAT's `zero_gap_extension_formula` is `matches!(task, "megablast" | "blastn")`. Not observable: for dc-megablast and blastn-short (dynamic programming) gap extension 0 is always an option error ('BLASTN gap extension penalty cannot be 0' or 'Greedy extension must be used if gap existence and extension options are zero', equal text and exit 1 in `-gapopen 0 -gapextend 0` and `-gapopen 3 -gapextend 0` runs), and for blastn too, so the formula prints only for megablast. Nothing to do.
- Row 9: confirmed from blast_nucl_options.cpp:251-271 and disc_nucl_options.cpp (no override): dc-megablast `min_diag_separation` 6; LOSAT uses `MIN_DIAG_SEPARATION_MEGABLAST` for dc via `task_defaults`.
- Row 61: `GetQueryBatchSize(.., use_default false)` returns 0 (adaptive mixer) when `BATCH_SIZE` is unset; `ProgramNameToEnum("dc-megablast")` = `eDiscMegablast`, `("blastn-short")` = `eBlastn` (StartsWith "blastn"). LOSAT's chunk sizes: `task_uses_megablast_chunks`.
- rmblastn: NCBI runs it (exit 0) and LOSAT rejects it (exit 2, 'the task rmblastn is not supported by LOSAT's BLASTN ...'): an explicit rejection kept.
- `-reward 0` with any task is an explicit LOSAT limit (`scoring.rs:171 check_losat_limits`; NCBI runs and prints the 'Could not calculate ungapped Karlin-Altschul parameters' warning per query, exit 0): not task specific.

## Extra rows (X1..X7)
X1 is GAP 1. X2 (`CBlastOptionsLocal::Validate`), X3 (`s_DiscWordOptionsValidate` / `LookupTableOptionsValidate` template branch), X4 (`BlastInitialWordOptionsValidate`) are on the path of the option errors and equal. X5 (`CreateTask` case-insensitive compare), X6 (`CDiscNucleotideOptionsHandle::SetWordSize` throw, never reached from the CLI) and X7 (`-num_threads` with `-subject`) are listed for completeness.
