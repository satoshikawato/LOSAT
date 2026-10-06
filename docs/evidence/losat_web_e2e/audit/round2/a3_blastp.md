# Round 2, angle (a) BLASTP, third pass: auditor's report

> S08+b. The final reply of the Sonnet auditor (read-only), as returned; the written copy is `a3/FINDINGS.txt`. Binary `6f070575…10bf91` (commit `3b84ce6c4`; the same bytes as the native of the final gate `run-20261004T163746Z`). Work dir `~/.cache/losat-web-gui-target/s08pb-audit/a3/` (`FINDINGS.txt`, `repros.tsv`, `w4_*.tsv`, `w4.py`, `vg.log`, `vg4.log`, `sweep*.tsv`, `out/*.raw`). Brief: [`brief/COMMON.md`](brief/COMMON.md), [`brief/ANGLE_A3.md`](brief/ANGLE_A3.md). How the findings were handled: [`../ROUND2.md`](../ROUND2.md).

Every comparison used identical argv and a clean environment and compared stdout, stderr and exit status. At most 3 search processes ran at once, and `/mnt/c` was never read. One old copy of the second pass's kept cases (`keep/`) was left in place because the safety check refused its `rm -rf`; it is harmless.

## Overall verdict for angle (a): SUPPORTED

R2A2-1 is fixed, and R2A2-2 is unchanged and verified. Nothing new was found.

## R2A2-1: FIXED

**Code reading.** `chain_blastp_init_hsps` (`LOSAT/src/algorithm/blastp/blast_engine.rs:3611-3618`) now returns early only for an empty list. That matches `BLAST_GetGappedScore` (`blast_gapalign.c:3707-3708`). The rest of the function still matches `s_ChainingAlignment`:
- the context walk, `blast_gapalign.c:3535-3558`;
- the chaining DP, `:3590-3626`;
- the drop test, `:3628-3636`;
- the sort, `:3658`.

When NCBI's drop test removes a lone HSP, the later loops in `BLAST_GetGappedScore` are bounded by `init_hitlist->total`, so a total of 0 is safe.

**The second pass's repros and the two new fixture pairs:** 160 comparisons, all SAME (`repros.tsv`). They covered:
- `a2/repro2/r2a2n1*`;
- the `e2e_one_hit` pair with `-task blastp-fast`, `-threshold` 1-30, `-window_size` 0/1/40/100 and `-evalue` 10/1e5;
- the new `e2e_fast_lone` and `e2e_fast_lone2` pairs in outfmt 0/6/7.

**New wide search of the chaining (`w4.py`):** 21,500 comparisons. 20,074 SAME, 1,424 LOSAT-REJECTS, 2 NCBI-CRASH, 0 DIFF.
- Every rejection is `-comp_based_stats 0`, which BLASTP rejects explicitly (`AUTHORITY.md` §K; `blast_engine.rs:1978`).
- The options were sampled from `-task blastp-fast`, `-threshold` 1-30 and default, `-window_size` 0/1/40/100/default, `-evalue` 1e-5/10/1000/1e5/1e8, `-comp_based_stats` 0/2/D/default, `-max_target_seqs` 1/5/default, and outfmt 0/6/7.

| Class | Comparisons | SAME | Rejected |
|---|---|---|---|
| (i) short random and planted pairs, 5-80 residues | 2,000 | 1,598 | 402 |
| (ii) each subject has exactly one planted ungapped alignment | 2,500 | 2,151 | 349 |
| (iii) multi-query batches, 1-copy and 2-4-copy queries mixed | 2,000 | 1,732 | 268 |
| (iii) dense short multi-query batches | 3,000 | 3,000 | 0 |
| (iv) 100-1,000 residue pairs | 1,000 | 860 | 140 |
| (iv) e2e inputs, `e2e_many_*` with 300 subjects included | 800 | 691 | 109 |
| (v) query split: one query of 9,801-19,900 residues, segments planted at the 9,800/9,900/10,000 boundaries | 100 | 84 | 16 |
| dense lone-HSP: queries of 5-40 residues, `-threshold` 1-13 | 9,000 | 9,000 | 0 |

**The fix is exercised.** Each case was also run with the second pass's binary. In 68 cases the old and new binaries differ, and in all 68 the new one is SAME with NCBI (the old one printed the extra row).

**Controls.** 1,100 cases ran without `blastp-fast`: `-task blastp`, and `-word_size 5` or 3. The new binary's output is byte-identical to the second pass's in all 1,100. Against NCBI: 958 SAME, 140 `-comp_based_stats 0` rejections, 2 NCBI-CRASH (R2A2-2).

**Before and after.** The second pass found 68 R2A2-1 differences in about 12,000 blastp-fast comparisons. This pass found 0 in about 28,000.

## R2A2-2: unchanged, NCBI-CRASH (pending the maintainer)

- **Crash rows:** NCBI exits 139 and LOSAT exits 0 in 23 cases: 21 from `w2` (14 in seed 21, 7 in seed 22) and 2 from `w4`.
- **Valgrind check:** under `valgrind -q`, LOSAT's stdout and exit status equal NCBI's in 23 of 23 (`vg.log`, `vg4.log`).
- **Code:** NCBI `blast_gapalign.c:3393-3416`, called at `:3937-3943`; LOSAT `blast_engine.rs:4285`.
- **Repro:** `blastp -query a3/repro2/r2a2n2_q.faa -subject a3/repro2/r2a2n2_s.faa -word_size 5 -window_size 0 -threshold 3 -evalue 1e5 -outfmt 7 -comp_based_stats t` (NCBI exit 139, LOSAT exit 0).

## Totals per class (second pass → this pass)

**Argv harness, 1,443 rows:**

| Class | Second pass | This pass |
|---|---|---|
| SAME | 994 | 983 |
| SAME except stderr | 12 | 10 |
| LOSAT-REJECTS | 369 | 382 |
| Parser exception | 63 | 63 |
| DIFF | 5 | 5 |

- The 13 rows that moved are giant `-threshold` values with `-word_size 5`. NCBI times out or aborts on them. LOSAT now gives its explicit overflow-bank rejection inside the time limit instead of also timing out, so only the class changed with machine load, not the behavior.
- The 5 DIFF are `-help` ×2 and `--help` ×2 (approved), and `-task blastp-fast -threshold +inf` (both sides time out).

**Other harnesses:**
- **Toolkit words (51):** 44 LOSAT-REJECTS, 3 parser exceptions, 4 SAME, as before.
- **Sweeps 1-18 (3,635):** 3,287 SAME, 38 SAME except stderr (`-num_threads` with `-subject`, approved), 310 LOSAT-REJECTS, 0 DIFF. Row for row the same as the second pass.
- **Random comparisons:** `cfuzz` 31/32, `ofuzz` 21/22 and `pfuzz2` 51 gave 8,271 comparisons, 129 rejected, 0 DIFF. The `pfuzz` panic search found 0 bad in 3,000.
- **`w1`:** 1,056 SAME, 444 rejected (word size 6/7), 0 DIFF.
- **`w1` fast:** 1,500 SAME (the second pass had 4 R2A2-1).
- **`g1`:** 1,056 SAME.
- **`w2` (5,000):** 4,510 SAME, 469 rejected (BLOSUM62 9/2), 21 NCBI-CRASH, 0 other DIFF (the second pass had 5 R2A2-1).
- **`w3` blastp-fast (5,000):** 5,000 SAME (the second pass had 57 R2A2-1).
- **`w2L` (800):** 800 SAME.
- **`g2` (256):** 256 SAME (the second pass had 2 R2A2-1).
- **`d12` (80) and `dd` (43):** identical to the second pass.
- **Query split, `CHUNK_SIZE` and tie-order sets:** in sweeps 1-18 and the `w4` split class, unchanged.

## Findings

- No new findings.
- Round 1 BLASTP findings:
  - SAME: BP-1, BP-3, BP-5.
  - LOSAT-REJECTS as decided: BP-2, BP-4 (D12), BP-8.
  - SAME or LOSAT-REJECTS as before: BP-6.
  - ACCEPTED: BP-7, BP-9, BP-10, BP-11.
  - FIXED: BP-12, BP-13.
- Pending the maintainer (not against support): D11, D12, D13, D14, R2A2-2.
- Not verified: D11 above 2^30, the web adapter path, TBLASTN and TBLASTX (other angles).
