# S08 linktie: TBLASTX equal-HSP exchange (rows 3612/3628, b5)

## Log (append-only)

### 1. Reproduction (done)
- `ncbi_full.tsv` / `losat_full.tsv` here are the full b5 outputs (4273 rows each); `diff` shows only rows 3612 and 3628.
- Pair: query nt 23655-23684 (frame +3, ctx 2) and 23653-23682 (frame +1, ctx 0), subject 1558-1587 (frame +1).
  Both have query protein offsets 7884..7894 and subject protein offsets 519..529, raw score 10 (22.1 bits). They are *fully comparator-equal*
  in s_RevCompareHSPsTbx (link_hsps.c:331-379): same context/3 group, same subject frame sign, same query.offset/end, same subject.offset/end.
  Only their `context` (0 vs 2) differs, which no comparator in link_hsps.c or ScoreCompareHSPs looks at.
- Result: link order of comparator-equal HSPs = their order in hsp_list->hsp_array (both qsorts are stable merge sorts in glibc 2.39).
- Smaller reproducer (query 1-25000 x subject 1-3000, `w/t1_25000_1_3000_*`) also shows the same pair reversed in the output (and several other
  equal-e-value groups ordered differently, same cause class). Full-size window shrink: 1-51000 x 1-51000 still shows the 4-line swap.

### 2. Root cause (proved)
**LOSAT never runs NCBI's `Blast_InitHitListSortByScore` (`score_compare_match`) on the per-chunk initial-HSP list, so HSPs that are
comparator-equal in every later sort (same context/3 group, same subject frame sign, same context-relative query/subject offsets and ends, same score) arrive
in link_hsps.c in *scan order* instead of NCBI's *absolute-query-offset (= context) order*.** Because link_hsps.c picks the best chain with `>=`, the last of two
equal candidates wins, so the arrival order decides which of the two HSPs gets the better chain (and hence which e-value).

NCBI call path (all lines in `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/`, CRLF, line numbers as `grep -n` prints them):
1. aa_ungapped.c:234 `Blast_InitHitListSortByScore(init_hitlist);` at the end of `BlastAaWordFinder` (the tblastx word finder, blast_engine.c:1019), once per subject chunk.
2. blast_extend.c:306-310 `Blast_InitHitListSortByScore` -> `qsort(..., score_compare_match)`; blast_extend.c:274-304 `score_compare_match`:
   score desc, `ungapped_data->s_start` asc, `length` desc, **`ungapped_data->q_start` asc**. `q_start` is the ABSOLUTE offset in the concatenated 6-frame query buffer
   (it is made context-relative later, in `s_AdjustInitialHSPOffsets`, blast_gapalign.c:4756-4758 via BLAST_GetUngappedHSPList 4719-4775), so two frames' copies of the
   same relative HSP are ordered by context (ctx 0 before ctx 2).
3. blast_gapalign.c:4719-4775 `BLAST_GetUngappedHSPList` walks init_hitlist in that order, `Blast_HSPListSaveHSP`, then `Blast_HSPListSortByScore` (blast_hits.c:1374-1381) which only calls
   qsort if `Blast_HSPListIsSortedByScore` (blast_hits.c:1355-1369) finds an inversion. `ScoreCompareHSPs` (blast_hits.c:1330-1353) only sees context-relative offsets, so the pair ties and keeps the init-list order.
4. link_hsps.c:485-486 `qsort(link_hsp_array, ..., s_RevCompareHSPsTbx)`; s_RevCompareHSPsTbx (link_hsps.c:330-378) compares `context/3`, sign of subject frame, query.offset desc, query.end desc,
   subject.offset desc, subject.end desc and never `context`: the pair is equal, so the (stable, glibc 2.39 merge) qsort keeps the arrival order (the NCBI binary here runs on Ubuntu glibc 2.39).
5. link_hsps.c:613-622 (`if(sum0>=max0)`, `if(sum1>=max1)`) and 759 / 887 (`if (new_sum >= maxscore)`): ties go to the LATER element of the sorted list. Chain removal then
   leaves the earlier element of the pair only the remaining (smaller) chain: first-selected chain length 12, other 11 (e-values 3.80e-99 and 1.62e-93).

LOSAT counterpart (commit 0d533ba76):
- run_impl.rs:2784 `init_hsps.push(init)` (scan order) -> run_impl.rs:2825 `get_ungapped_hsp_list(init_hsps, ...)` (blast_gapalign.rs:165-338): no `score_compare_match` sort exists anywhere in
  `algorithm/tblastx` (`grep -rn "score_compare_match\|InitHitListSortByScore"` finds nothing there, while blastp/blast_engine.rs:3017-3050 implements it and blastn/blast_engine/run.rs:10035 cites it).
- blast_gapalign.rs:337 `if !ungapped_hits_is_sorted_by_score_ncbi(...) { sort ... }` ports Blast_HSPListSortByScore faithfully, but with an unsorted-ness check on a comparator that ties
  the pair, so the scan order survives. The comment at blast_gapalign.rs:374-396 / linking.rs:257-276 ("NCBI delegates comparator-equal rows to the platform qsort") is right about the link
  sorts, but the incoming order those stable sorts preserve is wrong.
- linking.rs:195 `rev_compare_hsps_tbx`, linking.rs:620 `sort_hsps_by_ncbi_link_order` and the chain-selection code (DualMaximum, linking.rs:312-480) are correct: with the NCBI input order they reproduce NCBI exactly.

### 3. Proof (instrumented, private copies only)
NCBI side: an `LD_PRELOAD` qsort shim (`repro/shim/qshim.c`, identifies `s_RevCompareHSPsTbx` as libxblast.so base+0xe19d0 from `nm`) dumps the array handed to link_hsps.c's qsort.
LOSAT side: instrumented build of 0d533ba76 (`linktie_dump` / `linktie_full` in linking.rs; env `LOSAT_LINKTIE_DUMP`, `LOSAT_LINKTIE_FULL`; plus an env toggle `LOSAT_INITSORT` for the proposed sort).
Pair = (q 7884-7894 aa, s 519-529 aa, subject frame +1, raw score 42; ctx 0 = nt 23653-23682 and ctx 2 = nt 23655-23684):

| | pre-link array position | after s_RevCompareHSPsTbx | e-value after linking |
|---|---|---|---|
| NCBI (call 0) | 5294 ctx0, 5295 ctx2 | 2728 ctx0, 2729 ctx2 | ctx0 1.6225e-93, ctx2 3.8000e-99 (call 1: same) |
| LOSAT base | 5294 ctx2 (hsp_list_order 17), 5295 ctx0 (18) | 2728 ctx2, 2729 ctx0 | ctx2 1.6225e-93, ctx0 3.8000e-99 |
| LOSAT + sort | 5294 ctx0, 5295 ctx2 | 2728 ctx0, 2729 ctx2 | ctx0 1.6225e-93, ctx2 3.8000e-99 (= NCBI) |

LOSAT_TRACE_LINK_SELECTIONS: base selects chain_len=12 (3.800027e-99) with head q=23653-23682 (ctx0, idx 2073) in round 84 and chain_len=11 (1.622507e-93) with head q=23655-23684 (ctx2, idx 2072) in round 85;
with the sort the heads are exchanged (23655-23684 / 3.80e-99, 23653-23682 / 1.62e-93), as in NCBI. Files: `sel_base.txt`, `sel_fix.txt`.

Whole-array check (`full/`, `repro/*_prelink_order_call0.txt`; columns ctx q_start q_end s_start s_end s_frame score): NCBI vs LOSAT-base differ at exactly 2 positions of 6352 in the
prelink-link call (the b5 pair, plus a second comparator-equal pair `0 7882 7894 4189 4201 -2 45` that happens not to change the output) and 1 of 4275 in the second link call;
NCBI vs LOSAT+sort: 0 differing positions in both calls.

### 4. Proposed minimal change (not applied to any repo): `proposed_fix.diff` (89 changed lines incl. comments + a unit test), core:
diff -ru src_pristine/LOSAT/src/algorithm/tblastx/blast_engine/mod.rs src_fix/LOSAT/src/algorithm/tblastx/blast_engine/mod.rs
--- src_pristine/LOSAT/src/algorithm/tblastx/blast_engine/mod.rs	2026-10-02 18:33:45.000000000 +0900
+++ src_fix/LOSAT/src/algorithm/tblastx/blast_engine/mod.rs	2026-10-02 19:35:00.559694975 +0900
@@ -118,8 +118,8 @@
 
 // Import InitHSP and related functions from blast_gapalign module (NCBI blast_gapalign.c equivalent)
 pub(crate) use super::blast_gapalign::{
-    get_ungapped_hsp_list, sort_ungapped_hits_by_score_ncbi, trace_init_hsp_if_match,
-    ungapped_hits_is_sorted_by_score_ncbi, InitHSP,
+    get_ungapped_hsp_list, sort_init_hsps_by_score_ncbi, sort_ungapped_hits_by_score_ncbi,
+    trace_init_hsp_if_match, ungapped_hits_is_sorted_by_score_ncbi, InitHSP,
 };
 
 // Import subject scanning functions from blast_aascan module (NCBI blast_aascan.c equivalent)
diff -ru src_pristine/LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs src_fix/LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs
--- src_pristine/LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs	2026-10-02 18:33:45.000000000 +0900
+++ src_fix/LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs	2026-10-02 19:34:59.514694929 +0900
@@ -2809,6 +2809,15 @@
                 // ```
                 advance_tblastx_diag_offset(diag_offset, diag_array, window, chunk.length);
 
+                // NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:234-235
+                // ```c
+                // Blast_InitHitListSortByScore(init_hitlist);
+                // return status;
+                // ```
+                // BlastAaWordFinder sorts the chunk's init hit list by
+                // score_compare_match (blast_extend.c:274-310) before
+                // BLAST_GetUngappedHSPList sees it.
+                sort_init_hsps_by_score_ncbi(&mut init_hsps);
                 let mut hits = if init_hsps.is_empty() {
                     Vec::new()
                 } else {
diff -ru src_pristine/LOSAT/src/algorithm/tblastx/blast_gapalign.rs src_fix/LOSAT/src/algorithm/tblastx/blast_gapalign.rs
--- src_pristine/LOSAT/src/algorithm/tblastx/blast_gapalign.rs	2026-10-02 18:33:45.000000000 +0900
+++ src_fix/LOSAT/src/algorithm/tblastx/blast_gapalign.rs	2026-10-02 19:34:59.514694929 +0900
@@ -104,6 +104,58 @@
     }
 }
 
+// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:274-304
+// ```c
+// static int score_compare_match(const void *v1, const void *v2)
+// {
+//     ...
+//     if (0 == (result = BLAST_CMP(h2->ungapped_data->score,
+//                                  h1->ungapped_data->score)) &&
+//         0 == (result = BLAST_CMP(h1->ungapped_data->s_start,
+//                                  h2->ungapped_data->s_start)) &&
+//         0 == (result = BLAST_CMP(h2->ungapped_data->length,
+//                                  h1->ungapped_data->length)) &&
+//         0 == (result = BLAST_CMP(h1->ungapped_data->q_start,
+//                                  h2->ungapped_data->q_start))) {
+//         result = BLAST_CMP(h2->ungapped_data->length,
+//                            h1->ungapped_data->length);
+//     }
+//     return result;
+// }
+// ```
+// `q_start` is still the absolute offset in the concatenated query buffer here
+// (s_AdjustInitialHSPOffsets runs later, in BLAST_GetUngappedHSPList), so this
+// comparator orders HSPs that have identical context-relative coordinates in
+// different query frames by context. ScoreCompareHSPs (blast_hits.c:1330-1353)
+// cannot: it only sees context-relative offsets and ties.
+fn score_compare_init_hsps_ncbi(a: &InitHSP, b: &InitHSP) -> Ordering {
+    b.score
+        .cmp(&a.score)
+        .then_with(|| a.s_start.cmp(&b.s_start))
+        .then_with(|| (b.s_end - b.s_start).cmp(&(a.s_end - a.s_start)))
+        .then_with(|| a.q_start_absolute.cmp(&b.q_start_absolute))
+}
+
+// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:306-310
+// ```c
+// void Blast_InitHitListSortByScore(BlastInitHitList * init_hitlist)
+// {
+//     qsort(init_hitlist->init_hsp_array, init_hitlist->total,
+//           sizeof(BlastInitHSP), score_compare_match);
+// }
+// ```
+// NCBI reference: ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:234-235
+// ```c
+// Blast_InitHitListSortByScore(init_hitlist);
+// return status;
+// ```
+// (end of BlastAaWordFinder, the tblastx word finder, once per subject chunk.)
+// glibc's qsort is a stable merge sort and `slice::sort_by` is stable, so
+// comparator-equal records keep their scan order in both.
+pub(crate) fn sort_init_hsps_by_score_ncbi(init_hsps: &mut [InitHSP]) {
+    init_hsps.sort_by(score_compare_init_hsps_ncbi);
+}
+
 /// NCBI s_AdjustInitialHSPOffsets equivalent
 ///
 /// Reference: blast_gapalign.c:2384-2392
@@ -507,6 +559,24 @@
         }
     }
 
+    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_extend.c:274-310
+    // score_compare_match breaks score/s_start/length ties on the ABSOLUTE
+    // q_start, so equal context-relative HSPs of different query frames leave
+    // the word finder in context order whatever their scan order was.
+    #[test]
+    fn test_sort_init_hsps_orders_equal_relative_hits_by_absolute_query_start() {
+        let mut ctx2_first = make_init_hsp(20 + 7, 20 + 17, 519, 529, 42);
+        ctx2_first.ctx_idx = 2;
+        let mut ctx0_second = make_init_hsp(7, 17, 519, 529, 42);
+        ctx0_second.ctx_idx = 0;
+        let mut hsps = vec![ctx2_first, ctx0_second];
+        sort_init_hsps_by_score_ncbi(&mut hsps);
+        assert_eq!(
+            hsps.iter().map(|h| h.ctx_idx).collect::<Vec<_>>(),
+            vec![0, 2]
+        );
+    }
+
     #[test]
     fn test_get_ungapped_hsp_list_sorts_like_ncbi_score_compare_match() {
         let contexts = vec![make_context()];

### 5. Verification of the proposed change (private build `bin/LOSAT-fix`, built from 0d533ba76 + proposed_fix.diff, rustfmt clean)
- b5 full: byte-identical to NCBI (`repro/b5.losat-fix.tsv` == `repro/b5.ncbi.tsv`; base differs in rows 3612/3628).
- Windows of b5 where base != NCBI: `1-25000 x 1-3000` (38 diff lines), `1-51000 x 1-51000` (4), `1-80000 x 1-51000` (4), `1-51000 x 1-80000` (4), `1-25000 x 1500-1650` (96), `1-25000 x 1400-1700` (86),
  `1-25000 x 1001-1700` (58), `1-25000 x 1001-2200` (30), `1-24500 x 1001-2200` (2), `1-24000 x 1001-2200` (42), `1-25000 x 1-2000` (40): all 0 differences with the fix.
- No regressions (NCBI == base == fix, outfmt 6, default options): b1-b4, b6, b7, and LOSAT/tests/fasta pairs p02, p03, p04, p05, p07, p08, p09, p10, p12, p13 of tblastx_v010_parity_manifest.tsv
  (pre-link array identical to NCBI for b1-b4, b6, b7, p03, p07, p08, p12). p06 not completed (LOSAT base alone > 5 min on this loaded machine); regression sweep of the whole
  manifest (including non-default gencodes, multi-chunk subjects) should be done by the normal gates before merging.

- Unit tests (needs docs/evidence fixtures next to LOSAT/): `cargo test --lib tblastx::` 87 passed, including the new `test_sort_init_hsps_orders_equal_relative_hits_by_absolute_query_start`.

### 6. Notes / not checked
- Same bug class: any comparator-equal pair in the incoming order (a second one exists in b5: ctx 0/2 q 7882-7894 s 4189-4201 frame -2 score 45, no output effect there).
- blastp uses `sort_unstable_by` for its `score_compare_match` port (blast_engine.rs:3044-3050) while glibc qsort is a stable merge sort; irrelevant for this defect, but a latent tie-order risk. blastx has a doc-comment copy of score_compare_match (blastx/seed.rs:140-180) - not examined.
- NCBI side stability assumption: the oracle runs on glibc 2.39 (merge sort, stable). With the sort in place the full pre-link arrays match on 11 inputs, so the assumption holds empirically.

### 7. How to rerun
```
cd /home/kawato/.cache/losat-web-gui-target/s08/linktie/repro
N=/home/kawato/micromamba/bin/tblastx; L=/home/kawato/.cache/losat-web-gui-target/s08/bin/LOSAT-base; F=../bin/LOSAT-fix
diff <($N -query b5_q.fa -subject b5_s.fa -outfmt 6) <($L tblastx -query b5_q.fa -subject b5_s.fa -outfmt 6)          # rows 3612, 3628 (4 lines)
diff <($N -query small_q.fa -subject small_s.fa -outfmt 6) <($L tblastx -query small_q.fa -subject small_s.fa -outfmt 6)   # 96 lines (25 kb query x 151 nt subject)
diff <($N -query swap2_q.fa -subject swap2_s.fa -outfmt 6) <($L tblastx -query swap2_q.fa -subject swap2_s.fa -outfmt 6)   # two adjacent rows swapped
# $F (base + proposed_fix.diff) is identical to NCBI on all three
NCBI_LINKTIE_FULL=/tmp/x LD_PRELOAD=shim/qshim.so $N -query b5_q.fa -subject b5_s.fa -outfmt 6 >/dev/null   # NCBI pre-link arrays -> /tmp/x.0, /tmp/x.1
```
(small_/swap2_ are cuts of b5: query 1-25000 x subject 1500-1650, and query 1-24500 x subject 1001-2200; coordinates in the outputs are those of the cut, which starts at 1 for the query.)
