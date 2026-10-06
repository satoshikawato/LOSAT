# Round 2, angle (a) BLASTP, second pass: auditor's report

> S08+b. The final reply of the Sonnet auditor (read-only), as returned (the harness refused its report file; the written copy is `a2/FINDINGS.txt`). Binary `4d6c036f…8e3e` (the gate native of `8c905ccc5`, engine `75cc8e565`). Work dir `~/.cache/losat-web-gui-target/s08pb-audit/a2/` (logs `out/*.raw`, `sweep*.tsv`, `w1_*.tsv`, `w2_*.tsv`, `w3_*.tsv`, `wL3_41.tsv`, `g1.tsv`, `g2.tsv`, `d12*.tsv`, `dd.tsv`, `vg.log`, `an.log`, fuzz logs; repro inputs in `repro/`, `repro2/`). Brief: [`brief/COMMON.md`](brief/COMMON.md), [`brief/ANGLE_A2.md`](brief/ANGLE_A2.md). How the findings were handled: [`../ROUND2.md`](../ROUND2.md).

Every run used identical argv and a clean environment and compared stdout, stderr and exit status (and the created files for the `--` tests). No `vperf.lock` was present; never more than 3 processes at once.

## Overall verdict for angle (a): UNSUPPORTED (narrowly)

R2A-1 and R2A-2 are fixed, and D12 (extended) and the last-`--` rule behave as decided. The wide search for R2A-2 also exercised `-task blastp-fast`, whose NCBI preset turns on HSP chaining, and found one new reproducible silent defect:

- **R2A2-1:** with `-task blastp-fast` and an explicit low `-threshold` (13 or lower), LOSAT prints one extra short HSP that NCBI drops. It occurs in 68 of about 12,000 blastp-fast comparisons, never at the default threshold, and never in 800 comparisons with 100-500 residue sequences. It is not a regression: the first pass's binary has it too.
- **R2A2-2:** NCBI itself crashes on some valid commands, where LOSAT prints the result NCBI prints under valgrind.

Everything else matched NCBI or was an accepted or approved item.

## Verdicts on the first pass's findings

**R2A-1: FIXED.**
- The repros `repro/r2a1_q.faa` with `r2a1_s_A.faa` and `r2a1_s_NA.faa` (`-word_size 5` and `-task blastp-fast`, outfmt 6/7/0) are all SAME, exit 0, with outfmt 0 reports of 1443 and 1445 bytes.
- The fixture pair `e2e_protein_query.faa` with `e2e_short_mixed_subject.faa -word_size 5` is SAME (19,194 bytes).
- Wide search: `w1.py` 1,500 cases, the same with blastp-fast (1,500), and a 1,056-case grid of length x content x position. Subjects had 0-10 residues (random, X runs, mixes of X/B/Z/U/O/`*`, query substrings) and appeared alone or first, middle or last among normal subjects.
  - Options: `-word_size` 3/5, `-window_size` 0/1/40/100, thresholds, `-evalue` 10/1000/1e5, outfmt 0/6/7.
  - Result: no abort and no difference from a short subject.
  - `-word_size 6/7` is an explicit rejection ("compressed lookup table not supported"), which I accept as justified.

**R2A-2: FIXED.**
- The repros `repro/r2a2_q.faa` with `r2a2_s.faa` (`-word_size 5 -threshold 5 -window_size 0 -evalue 1e5`, outfmt 6 and 0) are SAME.
- At a realistic size, the e2e protein pair with `-word_size 5 -window_size 0 -threshold 12 -evalue 1e8 -outfmt 6` is SAME (2,626,578 bytes, about 15 s each). The `e2e_one_hit_*` fixture is SAME.
- Wide search: about 17,000 runs with no LOSAT abort.
  - Inputs: random pairs of 5-60 residues with planted hits at the sequence ends, low-complexity and repeat pairs, 1-5 queries, and e2e grids.
  - Options: word 3 with threshold 1-12, word 5 with threshold 1-15, blastp-fast with window 0, `-evalue` 10/1000/1e5/1e8, outfmt 0/6/7.

**D12: VERIFIED** (`d12.tsv`, `d12b.tsv`).
- blastp and tblastn were run on these file pairs:
  - blastp: `e2e_many_query.faa` with `e2e_many_subject.faa`; `e2e_protein_query.faa` with `e2e_protein_subject.faa`.
  - tblastn: `e2e_many_query.faa` with `e2e_many_subject.fna`; `e2e_protein_query.faa` with `e2e_tblastn_subject.fna`.
- Explicit LOSAT rejection (exit 1, "an -evalue of DBL_MAX or more"): `1.7976931348623157e308`, `1.7976931348623158e308`, `1.79769313486231570e308`, `+inf`, `1e999`. NCBI exits 0 on these four pairs.
- SAME, with byte-equal output: `1.7976931348623156e308`, `1.7e308`, `1e308`, `1e307`.
- Unsigned `inf` fails on both sides (NCBI 1, LOSAT 2; parser exception).
- The rejection is justified: on `e2e_protein_query.faa` with `e2e_many_subject.faa` (blastp) or `.fna` (tblastn), NCBI exits 139 for `1.7976931348623157e308` and `+inf`, and exits 0 for `1.7976931348623156e308` with output equal to LOSAT's.
- NCBI refs `blast_kappa.c:409,3687` and `blast_hits.c:3266`; LOSAT `blastp/blast_engine.rs:1917` and `tblastn/args.rs:508`.

**The `--` checks: VERIFIED** for blastn, blastp, tblastn and tblastx (35 argv each, `dd.tsv`). The class per argv is identical across the four programs.
- SAME:
  - `-outfmt 6 --`, `--` alone (outfmt 0), `-outfmt 7 --`, `-evalue 10 --`, `-outfmt=6 --`.
  - `-outfmt --` (both exit 1, identical stderr).
  - `-out --`, `-out -- -outfmt 6`, `-outfmt 6 -out --`, `-outfmt 0 -out -- --`, `-out=--`, `-out ---` (the created file `--` or `---` exists on both sides, with the same file list and bytes).
  - `-num_threads 2 --` (except for NCBI's thread warning).
- Parser exception, approved exception 1 (NCBI USAGE exit 1, LOSAT parser error exit 2): `-- -outfmt 6`, `-- --`, `-outfmt 6 -- --`, `-- extra`, `-outfmt 6 -- ""`, `-- -evalue 1`, `-- -help`, `-- -version`, `-- -dryrun`, `-outfmt 6 ---`.
- Both sides fail: `-evalue --` and `-max_target_seqs --` (NCBI 1, LOSAT 2).
- `-help --` prints NCBI's help vs LOSAT's help (approved). `-h --`, `-version --` and `-dryrun --` are explicit rejections (unported toolkit options).
- The first pass's only `--` difference (`-outfmt 6 --`, LOSAT exit 2) is gone. NCBI `ncbiargs.cpp:2866-2872`; LOSAT `cli.rs:120`.

## Round 1 BLASTP findings, re-run

| ID | Verdict |
|---|---|
| BP-1 | SAME (3/3) |
| BP-2 | LOSAT-REJECTS as decided (14/14; NCBI exit 139) |
| BP-3 | SAME (8/8) |
| BP-4 | LOSAT-REJECTS as decided (D12, now every value ≥ DBL_MAX); 1e10, 1e100, 1e308 SAME |
| BP-5 | SAME (33/33) |
| BP-6 | SAME for the 12 malformed values, LOSAT-REJECTS for the 10 values NCBI reads, as before |
| BP-7 | ACCEPTED |
| BP-8 | LOSAT-REJECTS as decided; `-outfmt 6 --` now SAME |
| BP-9 | ACCEPTED |
| BP-10 | ACCEPTED (2 stderr lines with `LOSAT_STARTUP_TRACE=1`) |
| BP-11 | ACCEPTED (approved): `-help`, `--help`; `-h` rejected |
| BP-12 | FIXED; the new comment citations checked against the NCBI source (`fasta.cpp:967`, `ncbiargs.cpp:2866-2872`, `blast_kappa.c:409,3687`, `blast_hits.c:3266`, `aa_ungapped.c:496-500`) |
| BP-13 | FIXED (code reading, unchanged) |

## Totals per class (first pass → this pass)

Argv harness, 1,443 executions:

| Class | First pass | This pass |
|---|---|---|
| SAME | 993 | 994 |
| SAME except stderr | 12 | 12 |
| LOSAT-REJECTS (NCBI runs) | 273 | 274 |
| LOSAT-REJECTS (NCBI also fails) | 96 | 95 |
| Parser exception | 63 | 63 |
| DIFF | 6 | 5 |

- The 5 DIFF are `-help` ×2 and `--help` ×2 (approved), and `-task blastp-fast -threshold +inf` (both sides time out at 60 s).
- The giant-`-threshold` rows (1e9, `+inf`, 2.14749e7, 2.1475e7, 3e7) change class between runs. Which side finishes inside 60 s on a loaded machine decides it: LOSAT's slow rejection or NCBI's abort. I do not treat them as findings.
- `-evalue 1.7976931348623157e308` is now a D12 rejection.
- Toolkit words (`tk.txt`, 51 argv): 44 explicit rejections, 3 parser exceptions, 4 SAME, no stray file.

Sweeps 1-18, 3,635 comparisons:

| Class | First pass | This pass |
|---|---|---|
| SAME | 3,246 | 3,287 |
| SAME except stderr | 38 | 38 |
| LOSAT-REJECTS | 310 | 310 |
| DIFF (the R2A-1 and R2A-2 aborts) | 41 | 0 |

Random comparisons (`cfuzz` 31/32, `ofuzz` 21/22, `pfuzz2` 51):
- First pass: about 10,700 compared, 15 DIFF.
- This pass: 8,271 compared, 129 rejected, 0 DIFF.
- `cfuzz.py 31` went from 13 DIFF to 0, and `ofuzz.py 21` from 2 DIFF to 0. The `pfuzz.py` panic search found 0 bad cases in 3,000.

New wide searches:

| Search | Cases | SAME | LOSAT-REJECTS | DIFF |
|---|---|---|---|---|
| R2A-1 short subjects, `-word_size` 3/5/6/7 | 1,500 | 1,056 | 444 (`-word_size 6/7`) | 0 |
| Same with `-task blastp-fast` | 1,500 | 1,496 | 0 | 4 (R2A2-1) |
| Short-subject grid | 1,056 | 1,056 | 0 | 0 |
| R2A-2 random short pairs | 5,000 | 4,505 | 469 (`-matrix BLOSUM62 -gapopen 9 -gapextend 2` pair unsupported) | 26 (21 NCBI SIGSEGV = R2A2-2, 5 = R2A2-1) |
| blastp-fast random short pairs | 5,000 | 4,943 | 0 | 57 (R2A2-1) |
| blastp-fast 100-500 residues | 800 | 800 | 0 | 0 |
| e2e grid | 256 | 254 | 0 | 2 (R2A2-1) |

The e2e grid used four inputs: the one_hit pair, the query vs `e2e_short_mixed_subject.faa`, the protein pair, and one query vs 300 subjects. Its options were word 5 with threshold 5-15, word 3 with threshold 8-12, and blastp-fast with threshold 5-13, over `-evalue` 10-1e8 and outfmt 6/0/7. The query-split, `CHUNK_SIZE` and tie-order sweeps came out the same as the first pass, minus the aborts.

## New findings

### R2A2-1 (medium, DEFECT): `-task blastp-fast` with an explicit low `-threshold`: LOSAT reports one extra short HSP that NCBI drops

NCBI:
- `c++/src/algo/blast/api/blast_options_handle.cpp:395-399`: blastp-fast sets word size 5, the compressed lookup and `SetChaining(true)`.
- `c++/src/algo/blast/core/blast_gapalign.c:3535-3558`: `s_ChainingAlignment` handles every context, including one holding a single ungapped HSP.
- `blast_gapalign.c:3628-3636`: the drop test `best_score - gap_score + word cutoff - 1 < hit cutoff` also applies to a node without partners.
- `blast_gapalign.c:3718-3724`: the call.

LOSAT: `LOSAT/src/algorithm/blastp/blast_engine.rs:3611-3613`, `chain_blastp_init_hsps`, has `if init_hsps.len() <= 1 { return init_hsps; }`. It returns before the drop test, so a lone ungapped HSP is never dropped.

Repros (inputs under `/home/kawato/.cache/losat-web-gui-target/s08pb-audit/`):
1. `blastp -query src2/LOSAT/tests/fasta/outfmt0/e2e_one_hit_query.faa -subject src2/LOSAT/tests/fasta/outfmt0/e2e_one_hit_subject.faa -task blastp-fast -window_size 0 -threshold 5 -evalue 1e5 -outfmt 6`
   - NCBI exits 0 with empty stdout; LOSAT exits 0 and prints `q s 40.000 5 3 0 1 5 1 5 3.3 7.3`.
   - The same pair with `-word_size 5` instead of `-task blastp-fast` is SAME, so the chaining causes the difference.
2. `blastp -query a2/repro2/r2a2n1_q.faa -subject a2/repro2/r2a2n1_s.faa -task blastp-fast -window_size 0 -threshold 6 -evalue 1000 -outfmt 6`
   - NCBI prints 10 rows; LOSAT prints the same 10 plus `q0 s1 100.000 1 0 0 5 5 11 11 24 8.1`.
3. With the default window: `blastp -query a2/repro2/r2a2n1b_q.faa -subject a2/repro2/r2a2n1b_s.faa -task blastp-fast -threshold 8 -evalue 1e5 -comp_based_stats D -outfmt 6`
   - LOSAT prints one extra row, `q0 s1 13.889 36 31 0 12 47 19 54 4.0 11.2`.

Observed in all 68 differences: LOSAT's rows are NCBI's rows plus exactly one extra row (two in one case; `an.log`).
- They occur only with an explicit `-threshold` of 1-13, never at the default 21 or at 14-25.
- The window can be any value, including the default; `-evalue` ranged from 10 to 1e8; the search spaces were tiny.
- A supporting experiment: when the subject of `r2a2n1b` is repeated twice (two ungapped HSPs), NCBI keeps both rows (2 rows with a 10-residue spacer, 3 with a 53-residue spacer) and LOSAT matches. With one copy, NCBI reports nothing and LOSAT reports one row.
- The root cause comes from code reading plus these experiments. I did not instrument NCBI or rebuild LOSAT without the early return.

Classification: DEFECT, a silent output difference with a narrow trigger.

### R2A2-2 (low; NCBI defect with a checkable valid result, pending the maintainer): NCBI SIGSEGV on valid one-hit searches with `-window_size 0 -threshold 1..3`

- **NCBI:** `c++/src/algo/blast/core/blast_gapalign.c:3393-3416` (`BlastGetStartForGappedAlignment`, called at `:3937-3943` with the Int4 length as Uint4). It reads the 11-letter window without a bound against the sequence ends.
  - Valgrind reports an invalid read of size 1, "0 bytes after a block of size 40" (query) and "of size 22" (subject), in this function.
  - The stray byte indexes `sbp->matrix->data[...]`, and gdb shows the SIGSEGV there.
- **LOSAT:** `LOSAT/src/algorithm/blastp/blast_engine.rs:4285`, `blastp_get_start_for_gapped_alignment_int4_length`, reads letters past a sequence as 0, the value a zeroed heap would give.
- **Repro:** `blastp -query a2/repro2/r2a2n2_q.faa -subject a2/repro2/r2a2n2_s.faa -word_size 5 -window_size 0 -threshold 3 -evalue 1e5 -outfmt 7 -comp_based_stats t`.
  - The query is SLFGLHDFLRHPLCWNGGWEAAHENYKEAHVGKNPES; the subject is NDCLPSNAEENYSEAPYTTR plus LKIQRGTGVMFHCQ.
  - `-comp_based_stats D`, `2` and `0` also crash; `-window_size 40` or `-threshold 11` do not.
- **Result:**
  - NCBI exits 139 with empty stdout: 3 of 3 runs, and with `MALLOC_PERTURB_` 0/1/85/170/255.
  - LOSAT exits 0 with 14 hits.
  - Under `valgrind -q`, NCBI exits 0 and its stdout is byte-identical to LOSAT's.
- **Rate:** 21 such cases in 5,000 random cases (all with threshold 1-3). All 21 are identical to LOSAT under valgrind (`vg.log`).
- **Classification:** an NCBI memory-safety defect. LOSAT prints the intended valid result, so this is a candidate approved exception under the NCBI defect policy; the maintainer needs to decide.

## Not verified

- D11 above 2^30 (allocations of several GB).
- The web adapter path.
- Whether R2A-1's shared `compressed.rs` fix is complete for TBLASTN and TBLASTX (other angles).
