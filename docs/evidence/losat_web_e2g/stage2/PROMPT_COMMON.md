# E2g stage 2: classification instructions (common part)

You classify, read-only, how faithfully LOSAT (a Rust port of NCBI BLAST) ports the NCBI functions of one range of NCBI files on the BLASTN path. Do not modify, build or commit anything. The only file you write is your range's TSV (path given in your range part).

## Sources
- NCBI C/C++ (the only authority): /mnt/c/Users/genom/GitHub/ncbi-blast/c++ (pinned commit 598d8ae6). Do not run `git status` there (very slow).
- LOSAT: /mnt/c/Users/genom/GitHub/LOSAT-web-gui/LOSAT/src (BLASTN engine: algorithm/blastn/; `blast_engine/run.rs` has ~12,900 lines, so grep it instead of reading it whole).
- Stage 1 tables (LOSAT's "NCBI reference" annotations mapped to NCBI functions): /mnt/c/Users/genom/GitHub/LOSAT-web-gui/docs/evidence/losat_web_e2g/stage1_functions.tsv (ncbi_file, ncbi_function, ncbi_start, ncbi_end, refs, losat_locations) and stage1_refs.tsv (each annotation). A function that is not in stage 1 may still be ported without an annotation; grep LOSAT for its name, its constants and its logic before calling it unported.
- Earlier authority notes (what was already compared and fixed, and the explicit rejections): /mnt/c/Users/genom/GitHub/LOSAT-web-gui/docs/evidence/losat_web_e2c/AUTHORITY.md (§A-§R; §G lists LOSAT's explicit rejections) and docs/evidence/losat_web_e2f/AUTHORITY.md (§A-§E: query batches, query split, preliminary hit lists).

## Path and option scope
NCBI `blastn` 2.17.0 run as bl2seq: `blastn -query Q -subject S` (never -db), tasks `blastn` and `megablast` only, with these options and any of their values: -word_size -num_threads -evalue -perc_identity -max_target_seqs -max_hsps -out -reward -penalty -gapopen -gapextend -dust -lcase_masking -subject_besthit -outfmt (0, 6, 7 only, without custom fields). Inputs are nucleotide FASTA with any IUPAC letters, lowercase, one or many queries and subjects. Environment variables BATCH_SIZE/CHUNK_SIZE/ADAPTIVE_CBS and .ncbirc are unset. Everything else (database search, other tasks or programs, other options, composition-based statistics, PHI, RPS, culling, -ungapped) is out of scope except where the in-scope path still calls a function (then classify that call).

## What to produce
Enumerate every NCBI function in your range's files that executes on that path (follow the calls from the entry points named in your range part; include static helpers and macros with logic). For a function whose behaviour branches on an option value, a task, the lookup table type, the input (e.g. ambiguity, lowercase, many queries, long query), write one row per branch (for example `BlastChooseNaExtend` per extension method, `s_SmallNaChooseScanSubject` per scan routine, the lookup table choice per table type). Then find LOSAT's counterpart, read both, and classify each row:

- `faithful`: same logic and the same C details: integer widths and wraparound (Int2/Int4/Uint4/size_t), float-to-int conversions (`(Int4)` of out-of-range values is INT_MIN on the x86-64 oracle), rounding, pointer and offset arithmetic (off-by-one, `get_data` positions), sort order and tie-breaks (the oracle's glibc 2.39 `qsort` is a stable merge sort; a Rust `sort_unstable*` where the NCBI comparator can return 0 for distinct elements is a difference), call order and timing. A LOSAT implementation that differs only for speed or memory but provably gives the same output (packed scans, parallel work reduced back into NCBI order) is `faithful` with note `equivalent: <how>`.
- `divergent`: ported but some logic or C detail differs. Say exactly what differs and which inputs show it.
- `unported`: on the path, no LOSAT counterpart, and it can change output, warnings, errors or exit status.
- `rejected`: LOSAT stops with an explicit "not supported by LOSAT" error before the inputs reach this branch (cite the LOSAT check, e.g. `scoring.rs:check_losat_limits`).
- `exception`: an approved project exception (only AGENTS.md's approved exceptions; none of them is a BLASTN behaviour, so expect none).
- `n/a`: executes on the path but cannot affect output, warnings, errors or exit status (allocation, freeing, pure bookkeeping). Give the reason.

Impact (how broad the inputs are that can show a difference; for faithful/n/a rows write `none`): `high` = default options or common inputs; `medium` = a common non-default option or input feature; `low` = rare inputs or extreme values; `none`.

## TSV format (tab-separated, one header line, then one line per row)
`range	ncbi_file	ncbi_lines	ncbi_function	branch	losat_location	status	impact	notes`
- ncbi_file relative to c++/ (e.g. src/algo/blast/core/na_ungapped.c); ncbi_lines `start-end` of the function or branch.
- branch: `-` when the function has no in-scope branching, else a short condition (e.g. `eSmallNaLookupTable, word > lut_word, aligned`).
- losat_location: `path:line(function)` relative to LOSAT/src, `;`-separated if several, `-` if none.
- notes: one line, concrete (what you compared; for divergent/unported the exact difference and an input that shows it if you can tell). No tabs or newlines inside a field.

## Working rules (important: usage limits)
- Append each row to the TSV as soon as you have classified it (open, append, close), so partial work survives if you stop. Create the file with the header line if it does not exist. If the file already has rows, continue after them without duplicating.
- Work function by function. Do not spend long on one function: if you cannot decide in reasonable time, write the row with status `divergent` or `unported` and notes starting `UNSURE:` explaining what to check.
- Never guess NCBI behaviour; read the NCBI source.
- When done, reply with: the number of rows per status, and the list of `divergent`, `unported` and `UNSURE` rows (function, branch, impact, one-line difference).
