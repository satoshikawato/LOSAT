# Angle (a), second pass: BLASTP after the round 2 fixes (and checks of the b/c/d fixes)

Read `COMMON.md` first; everything there holds, with these changes:
- LOSAT binary: `FINAL2=/home/kawato/.cache/losat-web-gui-target/s08pb-audit/bin2/LOSAT` (the final gate's native build of the new HEAD; SHA-256 in `bin2/LOSAT.sha256`). Source of that commit: `/home/kawato/.cache/losat-web-gui-target/s08pb-audit/src2/`. The fixes since the first pass: `ref/round2_fixes.diff`.
- Work dir: `/home/kawato/.cache/losat-web-gui-target/s08pb-audit/a2/`. The first pass's work dir `.../s08pb-audit/a/` (harness, repro inputs, `out/`, `sweep*.tsv`) is read-only for you: copy what you need into `a2/` and point LOSAT at `$FINAL2`.
- The first pass's report: `src2/docs/evidence/losat_web_e2e/audit/round2/a_blastp.md` (findings R2A-1, R2A-2 and its harness totals).

What changed (S08+b, after the first pass):
1. R2A-1: the compressed lookup's scan reads the letters past a subject shorter than the word as NULLB (`LOSAT/src/algorithm/tblastx/lookup/compressed.rs`, `scan_subject`; NCBI `aa_ungapped.c:496-500`, `blast_aascan.c:264-286`). Fixture `e2e.blastp.compressed_short_subject`.
2. R2A-2: a negative Int4 ungapped length (one-hit extension) is passed to the gapped start as NCBI's Uint4: only the first window of 11 letters is scored (`blastp/blast_engine.rs`, `blastp_get_start_for_gapped_alignment_int4_length`; NCBI `aa_ungapped.c:1054,1083`, `blast_gapalign.c:3394-3437`). Fixture `e2e.blastp.one_hit_negative_width`.
3. D12 extended: an `-evalue` of DBL_MAX or more is rejected for BLASTP and TBLASTN (`blast_kappa.c:409,3687`, `blast_hits.c:3266`); `1.7976931348623156e308` runs.
4. A last `--` changes nothing for BLASTN, BLASTP, TBLASTN and TBLASTX (`cli.rs`; NCBI `ncbiargs.cpp:2866-2872`); a word after `--` stays a parser error (exception 1); `--` as an option's value is that value.
5. Web ABI v1 BLASTP keeps its former error order and tabular-field rejection (not runnable here; ignore).
6. The BLASTP stable sorts of S08+b are unchanged since the first pass.

Do:
1. Re-run the first pass's whole BLASTP harness (round 1 argv, sweeps, random comparisons, query-split grids, CHUNK_SIZE grids, tie-order sets, toolkit words) with `$FINAL2`, and report totals per class next to the first pass's.
2. R2A-1: the first pass's repros, and a wide search: subjects of 0-10 residues alone and mixed with normal subjects (first, middle, last record), `-word_size 5/6/7` where LOSAT accepts them, `-task blastp-fast`, `-window_size` 0/1/40/100, `-threshold` values, `-evalue` 10/1000/1e5, outfmt 0/6/7, X/B/Z/U/O/`*` letters in the short subjects.
3. R2A-2: the first pass's repros, and a wide search for negative one-hit widths: `-window_size 0` with `-word_size 3` and `-threshold` 1-12, `-word_size 5` and `-threshold` 1-15, `-task blastp-fast -window_size 0`, `-evalue` 10/1000/1e5/1e8, on the e2e inputs, short random pairs (5-60 residues), low-complexity and repetitive pairs, multi-query batches (several contexts), queries and subjects whose hits end within 11 letters of the sequence end, outfmt 0/6/7. Keep NCBI runs under `timeout 600`.
4. D12: `-evalue` `1.7976931348623157e308`, `1.7976931348623158e308`, `1.79769313486231570e308`, `+inf`, `1e999` (rejected) and `1.7976931348623156e308`, `1e308` (run, equal to NCBI) for blastp and tblastn on `e2e_many_subject.faa`/`.fna` and on `e2e_protein_subject.faa`/`e2e_tblastn_subject.fna`.
5. The last `--` for blastn, blastp, tblastn and tblastx: `… -outfmt 6 --`, `… --` (outfmt 0), `-- -outfmt 6`, `-- --`, `-out -- …` (NCBI writes a file named `--`), `-- ` after `-help`; compare with NCBI (stdout, stderr, exit, files created).
6. Report per COMMON.md: verdict per R2A finding and per round 1 BLASTP finding, totals, new findings, and the overall verdict for angle (a).
