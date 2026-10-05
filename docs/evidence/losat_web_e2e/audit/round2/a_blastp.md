# Round 2, angle (a) BLASTP: auditor's report

> S08+b. The final reply of the Sonnet auditor (read-only), as returned. Binary `d8d18ec0…cb116` (the gate native of `491292327`). Work dir `~/.cache/losat-web-gui-target/s08pb-audit/a/` (`REPORT.md`, `out/*.raw`, `sweep*.tsv`, `timing1.tsv`, `x.tsv`, inputs in `repro/`). Brief: [`brief/COMMON.md`](brief/COMMON.md), [`brief/ANGLE_A.md`](brief/ANGLE_A.md). How the findings were handled: [`../ROUND2.md`](../ROUND2.md).

## Overall verdict: UNSUPPORTED (narrowly)

Two reproducible LOSAT aborts (exit 134) on valid commands where NCBI exits 0 with a normal result. Everything else matched NCBI (query splitting, tie order, options, environment, toolkit words): about 3,600 sweep comparisons, 10,700 random comparisons and 1,443 re-run round 1 argv, apart from the two defects and accepted or approved items.

## New findings

**R2A-1 (high, defect): blastp aborts on a subject of 1 or 2 residues when the compressed lookup is used with two-hit.**
- Triggers: `-word_size 5` or `-task blastp-fast`, with `window_size` ≥ 1 (the default 40 counts). `-window_size 0` is SAME.
- NCBI: `aa_ungapped.c:492-500` clamps `scan_range` so one scan always runs; `blast_aascan.c:264-286` primes the index by reading past the short subject into the sentinel or padding, without aborting.
- LOSAT: `LOSAT/src/algorithm/tblastx/lookup/compressed.rs:781-786` `subject.get(s..end).expect("NCBI BLAST compressed scan prime word must be in range")`; caller `LOSAT/src/algorithm/blastp/blast_engine.rs:1813-1858`.
- Repro: `A=~/.cache/losat-web-gui-target/s08pb-audit/a/repro; blastp -query $A/r2a1_q.faa -subject $A/r2a1_s_A.faa -word_size 5 -outfmt 6` (also `-task blastp-fast`, `-outfmt 0`, `$A/r2a1_s_NA.faa`).
- Observed: NCBI exit 0 (outfmt 6 empty, outfmt 0 a 1,423-byte no-hit report); LOSAT exit 134 with the panic above. One 1-2 residue record anywhere in the subject file is enough; all 38 short-sequence pairs with a 1- or 2-residue subject failed; lengths 3 to 8 SAME.

**R2A-2 (medium, defect): blastp aborts in one-hit mode with the compressed lookup and a low word threshold, from a negative ungapped HSP length.**
- Triggers: `-word_size 5 -window_size 0` with `-threshold` ≤ 12 and a permissive e-value.
- NCBI: `aa_ungapped.c:1033,1054,1083`: `init_hit_width = q_right_off - q_left_off + 1` can be ≤ 0, and `hsp_len` is a signed Int4, so the search continues.
- LOSAT: `blastp/extension.rs:317` ports the arithmetic; `blastp/blast_engine.rs:4280-4281` `usize::try_from(length).expect(...)` panics.
- Minimal repro: `blastp -query $A/r2a2_q.faa -subject $A/r2a2_s.faa -word_size 5 -threshold 5 -window_size 0 -evalue 1e5 -outfmt 6` (query DWSYG, subject DNRAG): NCBI exit 0 with `q s 40.000 5 3 0 1 5 1 5 3.3 7.3`; LOSAT exit 134.
- Realistic size: `-query LOSAT/tests/fasta/outfmt0/e2e_protein_query.faa -subject …/e2e_protein_subject.faa -outfmt 6 -word_size 5 -window_size 0 -threshold 12 -evalue 1e8`: NCBI exit 0 with 2,626,578 bytes (about 16 s); LOSAT exit 134.

## Round 1 BLASTP findings, re-run

| ID | Verdict | Evidence |
|---|---|---|
| BP-1 | SAME | `-subject hdronly.faa` outfmt 0/6/7, 3/3 |
| BP-2 | LOSAT-REJECTS as decided | 14/14 values in 1073741799..2147483623: NCBI 139, LOSAT 1 explicit rejection; 1073741798 SAME |
| BP-3 | SAME | 8/8 (2147483624..2147483647) |
| BP-4 | LOSAT-REJECTS as decided (D12) | `+inf`, `1e999`: NCBI 139 or exit 0 (word size 3), LOSAT 1 with the D12 text; 1e10, 1e100, 1e308 SAME |
| BP-5 | SAME | 33/33 `frames`/`sframe` lists |
| BP-6 | SAME for 12 malformed values; LOSAT-REJECTS for the 10 values NCBI reads | explicit rejection, exit 1 |
| BP-7 | ACCEPTED | rejections still come before NCBI's later checks |
| BP-8 | LOSAT-REJECTS as decided; `--` ACCEPTED | toolkit words below |
| BP-9 | ACCEPTED | wording |
| BP-10 | ACCEPTED | `LOSAT_STARTUP_TRACE=1` prints 2 stderr lines |
| BP-11 | ACCEPTED (approved) | `-help`, `--help` |
| BP-12 | FIXED | residual `cli.rs:15`, `cli.rs:172` references now correct; S08+b sort comments cite the right NCBI lines |
| BP-13 | FIXED (code reading) | `validate` goes through `check_options` → `validate_threads` |

Harness re-run (about 1,000 distinct argv; 1,443 executions): SAME 993; SAME except stderr 12 (NCBI thread warning; 2 where both timed out); LOSAT-REJECTS where NCBI runs 273; LOSAT-REJECTS where NCBI also fails 96; parser exception 63; DIFF 6 = `-help` ×2, `--help` ×2, `-outfmt 6 --` ×1, `-task blastp-fast -threshold +inf` ×1.

## Query splitting (S08+a)

Thresholds confirmed from `split_query_aux_priv.cpp:73-138`: chunk 10000, overlap 100, split when L/9900 ≥ 2 (batch ≥ 19800). About 700 comparisons, all SAME including stderr: 19 single-query lengths 9800-39600 at the 19799/19800/19801 and 29699/29700/29701 edges with homologous fragments at chunk boundaries (`-evalue 1000`, outfmt 0/6/7, `-seg yes`, `-task blastp-fast`); 15 multi-query layouts × BATCH_SIZE (unset, 1, 5000, 20000, 30000, 1e6); 70 long-homology cases; 150-1,500 tiny or short queries in one batch; `-max_target_seqs` 1-100000, `-max_hsps`, `-num_threads 4`, `-comp_based_stats` D and t; a 1,500-subject set; mixed lower case, X, `U B Z J *`. CHUNK_SIZE/OVERLAP_CHUNK_SIZE grids (about 1,400 pairs): NCBI succeeds in 1,693 and LOSAT matches every one; in all 150 where NCBI fails with exit 3, LOSAT rejects explicitly. Timing: 250,000-residue query with 1,500 fragments: LOSAT 3.1-3.7 s (16-17 s at `-evalue 1000`), NCBI 1.8-3.2 s (14.8-17 s); memory LOSAT 42-95 MB, NCBI 71-76 MB.

## Tie order (S08+b stable sorts)

About 1,100 stable-sort comparisons and 780 random ones, all SAME (tandem and interspersed repeats, low complexity, duplicated queries, identical subjects, self-search; `-seg no`, `-evalue 1000/1e5`, blastp-fast, `-window_size 0`, `-max_hsps`, `-max_target_seqs 1/2`, `-num_threads 4`). No run was confirmed to have an output-visible tie.

## Toolkit words and X-drop options

58 argv with toolkit words in value positions: all explicit rejections, no file created. `-out -help`, `-out=-version`, `-out -version=1` behave as NCBI. `-xdrop_gap`, `-xdrop_gap_final`, `-xdrop_ungap` are rejected as unported for BLASTP. `-word_size 5 -threshold 1e9` or `-task blastp-fast -threshold +inf`: LOSAT rejects after 47-53 s and 6.7 GB; NCBI aborts after 84 s and 5.1 GB (justified, slow rejection).

## Not verified

The D11 boundary above 2^30 (allocations of several GB); the web adapter path; whether R2A-1 reaches TBLASTN or TBLASTX; `-seg` with SEG per chunk beyond the SAME cases.
