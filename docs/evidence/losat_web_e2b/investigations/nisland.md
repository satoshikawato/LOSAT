# S08 nisland: TBLASTX gap3 (9-N island) pair exchange

## Log (append-only)

### 1. Reproduction + first localisation (done)
- gap3: NCBI vs current binary differ on rows 4/6 (swap of the two 62.0-bit HSPs between the 2.84e-53 / 1.37e-52 chains); LOSAT-base = NCBI.
- NCBI side dumped with a new qsort shim `shim/qshim2.c` (dumps every array handed to qsort for ScoreCompareHSPs / s_RevCompareHSPsTbx / s_RevCompareHSPsTransl / s_FwdCompareHSPsTransl; symbol offsets from `nm libxblast.so`): files `qs/g3.NNN.*`.
- LOSAT side dumped with the built-in `LOSAT_DUMP_TBLASTX_STAGE=<dir>` (stage_cur/).
- Call-0 pre-link array (13 HSPs; includes the random-resolved 2na-translation extents, e.g. ctx2 q67-107 s78-118 score 135 and ctx2 q0-74 score 312): LOSAT after_initial_ungapped_extension == NCBI ncbi_pre.0 exactly (all 13 rows, same order). So the CRandom bases, prelim scores/extents and init-hit sort are NOT the defect.
- Post first-link hsp_list order (NCBI g3.008.ScoreCompareHSPs.in, 13 rows) == LOSAT after_prelim_link_hsps_before_reevaluate (same order, same e-values). First link is NOT the defect.
- NCBI then runs a qsort(ScoreCompareHSPs) on that list BEFORE re-evaluation (g3.008.ScoreCompareHSPs.in/out): this is `Blast_HSPListSortByScore(hsp_list)` at the end of BLAST_LinkHsps (link_hsps.c:1802-1803), with the PRE-reevaluation scores. LOSAT has no such sort after its first link (it goes straight to retain + reevaluate).

### 2. Root cause (proved by NCBI qsort dumps and an A/B build)
**LOSAT never runs NCBI's `Blast_HSPListSortByScore` at the end of the FIRST `BLAST_LinkHsps` call (link_hsps.c:1802-1803), i.e. between the preliminary link and the re-evaluation. LOSAT goes link (chain order) -> prelim-evalue reap -> re-evaluate -> ONE score sort. NCBI goes link -> sort by (pre-re-evaluation) score -> reap -> re-evaluate -> sort again -> second link.**
The two sorts differ whenever re-evaluation changes scores/extents so that HSPs with different pre-re-evaluation scores tie afterwards. `ScoreCompareHSPs` (blast_hits.c:1330-1353) ties on (score, subject.offset, subject.end, query.offset, query.end), all context-relative and blind to context/frame, so tied HSPs of different query frames keep their incoming order (qsort is the stable merge sort of glibc 2.39). That order is the input order of the second BLAST_LinkHsps (blast_engine.c:1515-1520); s_RevCompareHSPsTbx (link_hsps.c:330-378) cannot separate them either (it never looks at `context` beyond context/3 nor the exact subject frame), and the `>=` of the chain choices (link_hsps.c:613-622, 759, 887) gives the better chain to the later one.

The gap3 pair (columns ctx q_start q_end s_start s_end s_frame score; NCBI `qs/g3.00{7,8,9}.*`, LOSAT `stage_cur/`):
- call-0 pre-link array (NCBI `ncbi_pre.0`) == LOSAT `after_initial_ungapped_extension` in all 13 rows, i.e. CRandom bases, 2na translation, prelim scores/extents and the init-hit sort are all identical. Contains `2 67 107 78 118 -3 135` (extent reached with the random bases of the island) and `0 72 107 83 118 -1 129`.
- post first-link hsp_list order == LOSAT `after_prelim_link_hsps_before_reevaluate` in all 13 rows (chain order): `2 0 74 8 82 -3 312`, `0 72 107 83 118 -1 129` (same chain), ..., `2 67 107 78 118 -3 135`, `5 0 35 7 42 1 129`, ...  -> chain order has ctx0(129) BEFORE ctx2(135).
- NCBI then qsorts this list by ScoreCompareHSPs with the pre-re-evaluation scores (`g3.008.ScoreCompareHSPs.in/out`): ctx2(135) now precedes ctx5(129) and ctx0(129).
- re-evaluation (blast_hits.c:2609-2737, Blast_HSPReevaluateWithAmbiguitiesUngapped blast_hits.c:676-733) trims `2 67 107 78 118 -3 135` to `2 72 107 83 118 -3 129` (the N-island bases are ncbi4na N in the re-evaluation, X codons) -> it ties with `0 72 107 83 118 -1 129`.
- NCBI second ScoreCompareHSPs sort (`g3.009.*`): in = `.. 2(72-107,129), 5(0-35,129), 0(72-107,129) ..`, out = `5, 2, 0` (5 first by subject.offset; 2 before 0 by incoming order). Input of the second link (`ncbi_pre.1` == `g3.010.RevCompareHSPsTbx.in`): ctx2 before ctx0.
- LOSAT (current binary) `before_link_hsps` has ctx0 before ctx2 (the single score sort starts from chain order, where ctx0(129) preceded ctx2). Second-link input differs at exactly these 2 of 12 positions -> chains exchanged.
- LOSAT-base matched NCBI only by accident (scan-order of the init hits happened to put ctx2 first); the init-hit sort of this session exposed it (it is not the cause).

NCBI lines (CR-stripped numbering, /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/): link_hsps.c:1761-1810 (BLAST_LinkHsps; 1802-1803 the sort), blast_engine.c:871-875 (first BLAST_LinkHsps in s_BlastSearchEngineCore), 899 (s_Blast_HSPListReapByPrelimEvalue), 1493-1497 (Blast_HSPListReevaluateUngapped), 1515-1520 (second BLAST_LinkHsps), blast_hits.c:2733-2734 (sort at the end of the re-evaluation), 1355-1381 (Blast_HSPListIsSortedByScore / SortByScore: qsort only if an inversion exists).
LOSAT lines (working tree LOSAT-web-gui, uncommitted): blast_engine/run_impl.rs:3193-3212 (first link, no sort), 3213-3214 (`prelinked_ungapped_hits.retain`), 3234 (reevaluate_ungapped_hsp_list, which ends with the sort of mod.rs:451-453) and 3261-3263.
The sort at the end of the SECOND BLAST_LinkHsps is already modelled by `report::final_hit_order` (report.rs:765-800: per (query, subject) `score_compare_hsps` stable sort); verified: also sorting `linked` right after the second link changes no output on 4000 inputs.

### 3. Proposed change (not applied): `proposed_fix.diff` (18 added lines incl. comments; against the working-tree run_impl.rs as of 20:15 JST)
```diff
--- a/LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs
+++ b/LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs
@@ -3204,6 +3204,24 @@   (after the `apply_sum_stats_even_gap_linking_with_parallel(...)` of the preliminary link, before the stage dump / retain)
+        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/link_hsps.c:1802-1803
+        // ```c
+        //     /* Sort the HSP array by score */
+        //     Blast_HSPListSortByScore(hsp_list);
+        // ```
+        // ... (comment: the re-evaluation may change/trim scores; ties then keep the order of this first sort, ...)
+        if !ungapped_hits_is_sorted_by_score_ncbi(&prelinked_ungapped_hits) {
+            sort_ungapped_hits_by_score_ncbi(&mut prelinked_ungapped_hits);
+        }
```
Both helpers are already imported in run_impl.rs (blast_engine/mod.rs:121-122). `sort_ungapped_hits_by_score_ncbi` itself re-checks `ungapped_hits_is_sorted_by_score_ncbi` (blast_gapalign.rs:437-440), so the explicit check only mirrors Blast_HSPListSortByScore. (Alternative: put the sort at the end of `apply_sum_stats_even_gap_linking_with_parallel`, where NCBI has it; not tested, would also move the `after_link_hsps_before_output_conversion` dump order.)
Clean (env-free) build of exactly this diff: `nisland/target/release/LOSAT` (source `nisland/LOSAT`).

### 4. Verification of the clean fix (all vs NCBI tblastx 2.17.0; cur = native/release/LOSAT, fix = nisland/target/release/LOSAT)
Per-corpus table: runs, cur != NCBI, fix != NCBI, fix != cur.
- gap3 (`-seg no`, outfmt 6): cur 4 diff lines, fix 0; gap3_s_aaa/acg/ccc/rand: all equal. Stage arrays of the fix binary == NCBI in both link calls (`gap3fix/`).
- tblastx_ambig_query x tblastx_ambig_subject, outfmt 0/6/7 x {default, -seg no}: 6 runs, 0/0/0 (this fixture does not expose the defect).
- rscratch_E `rnd/q1..700 x s1..700` (random coding-like DNA, low-complexity, N runs, IUPAC): `-seg no` outfmt 6 700 runs: 1 (seed 683 = gap3) / 0 / 1; outfmt 7 default 700: 0/0/0; earlier outfmt 0 `-seg no`: 1/0/1; default outfmt 6: 0/0/0.
- new `gen2.py` ambiguity corpus (seeds 1-3000; low-complexity inserts, N islands 1-15 and IUPAC runs both inserted and overwriting, near inserts, rc subjects, 1-3 subjects, 30% IUPAC queries): `-seg no` outfmt 6: 3000 runs, 162 / 0 / 162; same with outfmt 0: 162 / 0 / 162; default SEG outfmt 6 and 7: 18 / 1 / 17 (the 1 = seed 1898, separate SEG issue, section 5).
- `gen2.py` seeds 3001-5000 with random option mixes (-seg no, -query_gencode 2/4/5/11, -max_target_seqs, -evalue, -threshold, -window_size; outfmt 0/6/7): 2000 runs, 50 / 0 / 50.
- `gen3.py` (subject ACGT only, query always with IUPAC/N): 2000 runs `-seg no`: 1 / 0 / 1 (seed 428).
- `gen4.py` (NO ambiguity letters anywhere): 4000 runs `-seg no`: 3 / 0 / 3 (seeds 975, 2233, 3474). => the defect does not need ambiguity letters; re-evaluation also trims two-hit extension tails (e.g. seed 975 `0 8 29 18 39 -2 29` -> `0 8 13 18 23 -2 35`); the N island only makes it frequent. So "no change without ambiguity letters" is true for the LC738874/LC738875 data, not in general: the fix changes the output on 3/4000 ambiguity-free random inputs, always towards NCBI.
- Real data: LC738874 x LC738875, 40 random windows 6-30 kb x {no ambiguity, subject N islands/IUPAC, subject+query} (`genlc.py`), default SEG outfmt 6: 120 runs 0/0/0; `-seg no`: 120 runs 2 / 0 / 2 (w16_1, w16_2); the 40 ambiguity-free windows are byte-identical cur vs fix in both modes; six 100-200 kb windows with ambiguity, default SEG: 12 runs 0/0/0; whole LC738874 x LC738875 default outfmt 6, 0, 7: byte-identical cur == fix == NCBI (the `-seg no` whole-genome NCBI run did not finish in 17 CPU-min, killed).
- multi-query/subject corpora rscratch_E `bq1-4, mq1-25, biq1-40` outfmt 0/6/7: 182 + 40 runs 0/0/0.
- Gates in the private copy: `tests/tblastx_regression_fixtures.py check` 46 cases 0 differ; `docs/evidence/losat_web_e2a/check_losat.py --programs tblastx` 26 fixtures, 24 same + 2 approved code4 exceptions, 0 differ; `cargo test --release --lib tblastx` 92 passed.
- Suggested regression fixtures (inputs exist): `inventory/rscratch_E/bs/gap3_q.fa` x `gap3_s.fa` (`-seg no`, outfmt 6) and the ambiguity-free `rnd4/q975.fa` x `rnd4/s975.fa` (`-seg no`).

### 5. Separate pre-existing defect seen on the way (NOT this bug)
Default SEG, seed 1898 of `gen2.py` (`rnd2/q1898.fa` x `rnd2/s1898.fa`, outfmt 6/7): differs from NCBI in the current binary, the fix and LOSAT-base; it reproduces with every ambiguity letter of the subject replaced by A (`v1898/s_noamb.fa`) and with the single query `v1898/q_only.fa` (36 diff lines), and `-seg no` is equal. First divergence is the pre-link array (call 0): NCBI has `0 54 113 57 116 2 148`, LOSAT has `0 54 84 57 87 2 99` instead; all other 35 rows are identical. The query region aa 84-113 is the `GAAGAAGAAG...` (E-run) low-complexity stretch, so this looks like a SEG-mask difference of the translated query (not analysed further).

### 6. Correction to the LOSAT-base sentence in section 2
LOSAT-base (checked with its stage dump, `gap3base/`) does not resolve the island with CRandom: its preliminary list already holds the final extents (`2 72 107 83 118 -3 129` and `0 72 107 83 118 -1 129` both at 129 from the start, no `2 67 107 78 118 -3 135`), so no re-evaluation-induced tie arises and the single post-re-evaluation sort starting from the init-hit order gives NCBI's order. It matches NCBI here because the missing first sort has nothing to reorder; the CRandom change (correct) makes the re-evaluation shorten HSPs and thereby exposes the missing sort. The init-hit sort (score_compare_match) is not involved: gap3 stage arrays are identical to NCBI up to and including the post-link order.
