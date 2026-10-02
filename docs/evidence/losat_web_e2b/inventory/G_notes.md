# Range G notes: how NCBI tblastx (`-subject`, ungapped, sum statistics) limits the subjects per query to the hit list size

All NCBI paths are relative to `c++/` (pinned 598d8ae6), line numbers CR-stripped. LOSAT paths are relative to
`LOSAT/src`. Oracle: `/home/kawato/micromamba/bin/tblastx` 2.17.0 (env unset). LOSAT binary:
`s08/bin/LOSAT-base`. Inputs and outputs of every oracle run are in `scratch_G/` (names below). TSV: `G.tsv` (32 rows).

## 0. The answer in six lines

1. The hit list size N is `-max_target_seqs` (default 500). It is written by `SetHitlistSize` at
   `blast_args.cpp:2978` (outfmt 6/7: `2960-2962`, outfmt 0: `2924-2928`). For tblastx the preliminary hit list size is
   exactly N (`GetPrelimHitlistSize`, `blast_hits.c:44-71`: no composition-based statistics, not gapped, so N is returned unchanged).
2. The only step that drops subjects is `Blast_HitListUpdate` (`blast_hits.c:3243-3299`), called by the HSP collector
   (`hspfilter_collector.c:149-161`) once per subject per query during the preliminary search. The traceback stage of an ungapped
   search (`blast_traceback.c:1505-1707`) re-inserts the kept lists into a second hit list of the same size N (never overflows) and
   `s_BlastPruneExtraHits` (`blast_traceback.c:877-893`) is a no-op.
3. The kept set is the N best subjects of a query by the total order `s_EvalueCompareHSPLists`
   (`blast_hits.c:3078-3107`): best e-value ascending (all values below 1e-180 are equal), then `hsp_array[0]->score` descending,
   then **oid descending (the later subject in the file wins a full tie)**. It does not depend on the order in which subjects arrive.
4. The comparator reads `hsp_array[0]->score` of the list **as it is**. In `Blast_HitListUpdate` every list is put in e-value order
   (`Blast_HSPListSortByEvalue`) as soon as the query's hit list overflows (more than N subjects with hits); before that, lists stay in score
   order (`BLAST_LinkHsps` ends with `Blast_HSPListSortByScore`, `link_hsps.c:1803`). The last subject sort
   (`Blast_HSPResultsSortByEvalue`, `blast_traceback.c:1770-1772`, active because `-subject` runs in database-scan mode) uses the same
   comparator on the same state. So: no overflow (no subject dropped) -> order by the highest-score HSP; overflow -> order and kept set by the
   lowest-e-value HSP's score. The two differ only when best e-values tie (equal doubles, or both below 1e-180) and the lists' highest score is
   not the first score after the e-value sort.
5. Consequently the kept set equals "the first N subjects of the default run's order" except in those ties (oracle: `tie2AB.fa`, `tie3*.fa`):
   there the default run (no overflow) prints SA before SB but `-max_target_seqs 1` keeps SB.
6. LOSAT parses `-max_target_seqs` and ignores it (all subjects kept; no stderr warning for N < 5), and its subject order uses a different
   state for `hsp_array[0]` than NCBI does (git HEAD/LOSAT-base: always the e-value state; the uncommitted edit in the working tree: always the
   score state).

## 1. The path, step by step, with the state of the HSP list

| # | NCBI step (file:line) | HSP list state after the step | `best_evalue` | `hsp_array[0]` |
|---|---|---|---|---|
| 1 | `CFormattingArgs::ExtractAlgorithmOptions`, `blast_args.cpp:2895-2978`: N from `-max_target_seqs`, else 500 (`blast_options.c:1453`); `opt.SetHitlistSize(N)`; warning if N < 5 (`2975-2977`) | - | - | - |
| 2 | `CreateHspWriter`, `setup_factory.cpp:343-358`: collector writer with `prelim_hitlist_size = GetPrelimHitlistSize(N, 0, FALSE) = N` (`hspfilter_collector.c:336-337`) | - | - | - |
| 3 | `BLAST_PreliminarySearchEngine`, `blast_engine.c:1412-1478`, subjects in increasing oid; per subject `s_BlastSearchEngineCore`: six frames, `BLAST_LinkHsps` (`870-873`), `s_Blast_HSPListReapByPrelimEvalue` (`899`, cutoff `-evalue`, `blast_parameters.c:868`) | score order (`link_hsps.c:1803`); reap keeps order | min over HSPs (`link_hsps.c:1806-1810`), not recomputed after the reap | highest raw score |
| 4 | same function, `1480-1542` (ungapped): `Blast_HSPListReevaluateUngapped` (re-sorts by score, `blast_hits.c:2733-2734`), `BLAST_LinkHsps` again (`1515-1520`), prelim e-value reap (`1535`), query coverage (`1538`), bit scores (`1541`) | score order | recomputed by the second `BLAST_LinkHsps` | highest raw score |
| 5 | `BlastHSPStreamWrite` (`1554`) -> `s_BlastHSPCollectorRun` (`hspfilter_collector.c:82-168`): HSPs split by query (order kept); empty list freed | score order | - | highest raw score |
| 6a | `Blast_HitListUpdate`, hit list not full (`blast_hits.c:3252-3266`): `best_evalue = s_BlastGetBestEvalue` (`3246`), list appended, no comparison | score order | min over HSPs | highest raw score |
| 6b | first overflow (`3269-3279`): all stored lists `Blast_HSPListSortByEvalue` + `best_evalue` recomputed, `s_CreateHeap` (root = worst by `s_EvalueCompareHSPLists`) | e-value order (e-value class, then `ScoreCompareHSPs`) | min over HSPs | HSP with the lowest e-value (highest score inside the lowest class) |
| 6c | full: new list `Blast_HSPListSortByEvalue`, compared with the root (`3284-3296`); root strictly better -> new list freed, else root replaced and sifted down | kept lists: e-value order | min | lowest-e-value HSP |
| 7 | `BlastHSPStreamClose` (`blast_hspstream.c:133-207`): lists of all queries sorted by decreasing oid, read from the end (`BlastHSPStreamBatchRead`, `568-614`) -> traceback sees increasing oid | unchanged | unchanged | unchanged |
| 8 | `BLAST_ComputeTraceback_MT` ungapped branch (`blast_traceback.c:1526, 1691-1707`): no traceback (`perform_traceback = gapped_calculation = FALSE`); `Blast_HSPListGetBitScores`; `Blast_HSPResultsInsertHSPList(..., hitlist_size = N)` -> second hit list, never full | unchanged (score order if no overflow in step 6, else e-value order) | recomputed (same value) | as in step 6 |
| 9 | `Blast_HSPResultsSortByEvalue` (`blast_traceback.c:1770-1772`; `blast_hits.c:3383-3402`), only because `BlastSeqSrcGetTotLen > 0` (database-scan mode of `-subject`, `blast_app_util.cpp:204-210`, `seqsrc_multiseq.cpp:175-181`): qsort of the subjects with `s_EvalueCompareHSPLists` | unchanged | - | as in step 6 |
| 10 | `s_BlastPruneExtraHits` (`1777`): no-op for tblastx | - | - | - |
| 11 | `BlastHitList2SeqAlign_OMF` (`blast_seqalign.cpp:1569-1577`): subjects in hit-list order; HSPs of each subject `Blast_HSPListSortByEvalue` **after** the subject order is fixed | e-value order | - | - |
| 12 | `x_PrintTabularReport` (`blast_format.cpp:813`) -> `PruneSeqalign(..., m_HitlistSize = N)`: first N runs of different subject ids (no-op after step 6) | - | - | - |

The hit list counts subjects (one `BlastHSPList` per subject per query), only subjects that have at least one HSP with e-value <= `-evalue`
after the second `BLAST_LinkHsps` reap, and each query has its own list (a 12000-nt query set that NCBI searches in two batches still gets N per query).

`hit_params->low_score` is NULL (`low_score_perc` is 0, `blast_parameters.c:819-822`), `worst_evalue` and `low_score` of the hit list are never read
for ungapped tblastx; query splitting is off for ungapped searches (`split_query_cxx.cpp:60-61`), so no `Blast_HitListMerge`;
`hsp_num_max` is 0 -> `INT4_MAX` (`blast_hits.c:213-228`): no per-subject HSP cap.

### Why the kept set does not depend on arrival order
The comparator is a strict total order over the lists of one query (oids are distinct) and, once the hit list overflows, every list that is
ever compared is in e-value order (stored lists are sorted at heap creation, new lists before the comparison). A keep-the-best-N heap over a total
order keeps the same N lists whatever the order. The state of the kept lists afterwards is also order independent: all e-value order if the
query had more than N subjects with hits, else all score order. Therefore LOSAT may run its (parallel) subject loop as now and apply the hit list
afterwards in oid order.

## 2. -max_hsps and hsp_num_max (question 4)
* NCBI tblastx **has** `-max_hsps <int >= 1>` (`blast_args.cpp:203-207`, `317-319`, listed by `tblastx -help`). It sets `max_hsps_per_subject` (default 0 = off), applied only in
  the traceback stage by `s_FilterBlastResults` -> `Blast_TrimHSPListByMaxHsps` (`blast_traceback.c:1763-1766`, `849-851`; `blast_hits.c:2049-2069`), which keeps the first
  `max_hsps` HSPs of every kept subject **in the list's current order** (score order, or e-value order if the query's hit list overflowed).
* `hsp_num_max` (`BlastHitSavingOptions::hsp_num_max`) is a different field, 0 by default, never set by any tblastx option, `BlastHspNumMax` returns `INT4_MAX`.
* LOSAT: `-max_hsps` is rejected (`unknown option`); out of scope per COMMON.md. Oracle: `tblastx -max_hsps 2` prints 2 rows per subject.

## 3. Oracle runs (all in `scratch_G/`; helper scripts `cmp.sh`, `summ.py`, `chk.py`, `mk*.py`)
Query `Q.fa` = `LC738874` positions 180001-183000 (id `Q1`). `multi.fa` = 18 windows of 2000 nt of `LC738875` (240001-242000, 242001-244000, 244001-246000,
246001-248000, 248001-250000, 241001-243000, 243001-245000, 245001-247000, 100001-102000, 50001-52000, 300001-302000, 346001-348000, 266001-268000, 273001-275000,
276001-278000, 10001-12000, 200001-202000, 330001-332000); 14 of them have hits, 4 none.

Commands: `tblastx -query Q.fa -subject S -outfmt 6 [-max_target_seqs N]` and `LOSAT-base tblastx -query Q.fa -subject S -outfmt 6 [-max_target_seqs N]`.

| Run | NCBI | LOSAT-base |
|---|---|---|
| `multi.fa`, default | order W05 W02 W00 W07 W06 W11 W01 W15 W03 W12 W13 W14 W16 W08 (14 subjects, 700 rows) | byte-identical (md5 d243019e...) |
| `multi.fa`, N=1 / 2 / 3 / 5 | 1 / 2 / 3 / 5 subjects = the first N of the default order, rows byte-identical to the default run's rows of those subjects; stderr `Warning: [tblastx] Examining 5 or more matches is recommended` for 1,2,3 (none for 5) | 700 rows for every N, no stderr |
| `multi_shuf.fa` (same subjects, shuffled file order) | same subject order as `multi.fa`, same prefixes | same as above |
| `Q2.fa` (two queries), N=1,2,3 | per query: Q1 W05; Q2 W10 (N=1); +W02 / +W06 (N=2); +W00 / +W01 (N=3) | all |
| `Q3.fa` (two 6000-nt queries = 2 batches), N=2,3 | 2 / 3 subjects per query | all |
| `zero.fa` (subjects Ca,D1,Cb,Cc,D2,Cd,W05,W02; Ca=Cb=Cc=Cd identical copies of the query window, D1/D2 1% mutated copies; E = 0.0 for the first six) | default `Cd Cc Cb Ca D1 D2 W05 W02` (first bit scores 1397x4, 1376, 1333); N=1 `Cd`; N=2 `Cd Cc`; N=3 `Cd Cc Cb`; N=5 `Cd Cc Cb Ca D1` | default byte-identical; N ignored |
| `copies600.fa` (600 identical copies `C000..C599` of a hitting 1000-nt window), default | 500 subjects `C599 .. C100` (first = last file record; the 100 lowest oids dropped) | 600 subjects |
| `copies300.fa`, default | 300 subjects `C299 .. C000` | same |
| `copies600.fa`, `-outfmt 0` default | 500 description rows, 250 `>` alignment blocks (`copies300.fa`: 250 blocks) | not ported |
| `copies600.fa`, `-outfmt 0 -max_target_seqs 7` | 7 description rows, 7 blocks | not ported |
| `tie2AB.fa` (SA oid 0, SB oid 1; equal best e-value 2.22e-167; SA has a 140-bit HSP with e-value 5.60e-80, first HSP after e-value sort 137 bits; SB 137/137), default | `SA SB` | `SB SA` (differs) |
| `tie2AB.fa`, N=1 | `SB` (higher oid wins the tie in the e-value state) | all |
| `tie2AB.fa`, N=2 | `SA SB` | all |
| `tie2BA.fa` (SB oid 0, SA oid 1), default / N=1 / N=2 | `SA SB` / `SA` / `SA SB` | default `SA SB`; N ignored |
| `tie3.fa` (SA0, SB1, SC2 weak), default / N=3 / N=2 / N=1 | `SA SB SC` / `SA SB SC` / `SB SA` / `SB` | `SB SA SC` for all (differs) |
| `tie3b.fa` (SC0, SA1, SB2) | `SA SB SC` / `SA SB SC` / `SB SA` / `SB` | `SB SA SC` |
| `tie3c.fa` (SB0, SC1, SA2) | `SA SB SC` / `SA SB SC` / `SA SB` / `SA` | `SA SB SC` |
| `tie3.fa`, `tie3c.fa` with `-culling_limit 100`, N=1/2/3 | tie3: `SB` / `SA SB` / `SA SB SC`; tie3c: `SA` / `SA SB` / `SA SB SC` | not compared |
| `-max_target_seqs 0`, `-1`, `abc`, `2147483648` | USAGE + error, exit 1 | error, exit 2 |
| `-max_target_seqs 2147483647` | 700 rows, no warning | 700 rows |
| `-outfmt 7 -max_target_seqs 2` | header `# 107 hits found` (HSPs of the two kept subjects), stderr warning | not ported |
| empty query, N=2 | stderr `Warning: [tblastx] Examining 5 or more matches is recommended` then `Warning: [tblastx] Query is Empty!`, exit 0 | - |

Reading of the tie runs: with `tie3.fa` and N=2 (three subjects have hits, so the hit list overflows) the e-value state decides and the later subject wins
(`SB SA`); with N=3 or the default nothing overflows, lists are in score order and the 140-bit HSP puts SA first (`SA SB SC`). The kept set for N=1 and N=2 is
independent of the file order (checked with three file orders).

## 4. Note on LOSAT's working tree (changing while this inventory ran)
The repository working tree is being edited (uncommitted; files modified 17:14-17:20 on 2026-10-02); `LOSAT-base` and git HEAD do not have these edits:
* `common.rs:899`: `best_score` changed from `hsp_hits.first()` after the e-value sort (HEAD line 913) to `hsp_hits.iter().map(raw_score).max()` with the comment
  "an HSP list is still in score order". That is NCBI's `hsp_array[0]` when the hit list never overflowed, so it fixes the `tie2AB`/`tie3` **default** orders,
  but it is wrong for the kept subjects of an overflowed query (step 6b-9 of the table above).
* `algorithm/tblastx/args.rs:75`: `max_target_seqs: Option<usize>`; `run_impl.rs:715-725` `process_options()` writes `few_matches_warning("tblastx")` for N < 5
  (called before `Query is Empty!`, lines 735/782, which is NCBI's order); `run_impl.rs:977` calls `report.rs:658 prune_to_hitlist_size(hits, N)` on the output of
  `ncbi_order_evalue_hsp_order`, i.e. **"the first N subjects of the final order"**, labelled as `s_BlastPruneExtraHits`.
  This simplification gives NCBI's rows for every oracle case without an exact tie (multi.fa, Q2, Q3, copies600 with all-equal scores keeps the highest 500 oids
  because the order puts higher oids first) but not for the ties of question 3: `tie3.fa` N=1 keeps SA (NCBI: SB) and prints N=2 as `SA SB` (NCBI: `SB SA`);
  `tie2AB.fa` N=1 keeps SA (NCBI: SB). The reason is that in NCBI the cut is made by `Blast_HitListUpdate` in the preliminary stage with lists in e-value
  order, not by `s_BlastPruneExtraHits` (a no-op for tblastx) on the final score-order sort. Suggested regression inputs: `scratch_G/tie3.fa`, `tie3b.fa`, `tie3c.fa`,
  `tie2AB.fa`, `tie2BA.fa`, `zero.fa`, `copies600.fa` with `Q.fa` (expected outputs: the `*.ncbi.*.tsv` files next to them).

## 5. What TBLASTX needs (question 5)
NCBI functions to port (all `HitList`-level pieces already exist generically in `algorithm/blastn/hsp.rs`):
* `GetPrelimHitlistSize` -> `get_prelim_hitlist_size(n, false, false)` (`hsp.rs:508`) = N.
* `s_BlastHSPCollectorRun` (`hspfilter_collector.c:82-168`): per subject, per query, the HSP list in score order; model `collect_prelim_hit_lists` (`blastn/blast_engine/run.rs:3489-3520`).
* `Blast_HitListNew` / `Blast_HitListUpdate` / `s_BlastHitListInsertHSPListInHeap` / `s_CreateHeap` / `s_Heapify` / `s_EvalueCompareHSPLists` / `s_EvalueComp`
  -> `HitList<L>::update`, `compare_hsp_lists`, `evalue_comp` (`hsp.rs:1023, 793, 549`), with a new `HitListEntry` impl for a TBLASTX subject list of `Hit`s
  (oid = `Hit.s_idx`, `update_best_evalue` = min `e_value`, `first_score` = `raw_score` of the first HSP in the current order, `sort_by_evalue` = `Blast_HSPListSortByEvalue`
  with `common.rs:380`; the list must start in `score_compare_hsps` order, `common.rs:283`).
* `BlastHSPStreamClose` + traceback loop (`blast_hspstream.c:133-207`, `blast_traceback.c:1505-1707`): second hit list of size N filled in increasing oid;
  state unchanged. Then `Blast_HSPResultsSortByEvalue` (`HitList::sort_by_evalue`, `hsp.rs:1096`) and `s_BlastPruneExtraHits` (`prune_by_size`, `hsp.rs:1115`, no-op).
* The output order must come from that final sort, or `ncbi_query_subject_groups` (`common.rs:795-905`) must take, per query, a flag "hit list overflowed" and use
  `best_score = max raw_score` if not, `best_score = first HSP's score after the e-value sort` if so. The HSP order inside each subject stays the e-value sort of
  `blast_seqalign.cpp:1577` (`common.rs:868-875`).

Where: `algorithm/tblastx/blast_engine/run_impl.rs`, function `search_query_batch` (line 871). The final HSPs of a subject are complete after the loop `for h in linked`
(built at `let out_hit = Hit {` line 3009, e-value filter at line 2847, pushed to `tx.send(final_hits)` / `st.hits.extend(final_hits)` at lines 3154-3156).
NCBI's `BlastHSPStreamWrite` of that same list is `blast_engine.c:1554`, i.e. after the second `BLAST_LinkHsps` + e-value reap, before the next subject.
Because the result is independent of arrival order, the faithful and simplest place is **after the parallel subject loop**, at `let hits = final_hits.unwrap_or_default();`
(line 3424), before `TblastxBatch { hits, .. }` is returned: `Hit.q_idx` is still batch-local (NCBI's `query_index`; `run_in_pool` adds the batch start at lines 848-852),
`Hit.s_idx` is the oid. There: group by (q_idx, s_idx) -> per query `HitList::new(N)`; `update()` in increasing s_idx; second pass / sort / prune; flatten in final order;
keep the per-query overflow flag for the writer.

Also for TBLASTX: `few_matches_warning("tblastx")` (`report/query_warnings.rs:119-123`) on stderr when N < 5 at option-extraction time (after the subject FASTA diagnostics, before
`Query is Empty!` and query diagnostics); outfmt 0 display caps (descriptions 500, alignments 250 unless `-max_target_seqs`; `algorithm/blastx/args.rs:813-921` resolves them for BLASTX).

## 6. What the port must get right (each item with its NCBI source)
1. Hit list size N per query, N = `-max_target_seqs` or 500, prelim size = N for tblastx, counting subjects that still have an HSP after the final `-evalue` reap
   (`blast_args.cpp:2960-2978`, `blast_hits.c:44-71`, `blast_engine.c:1535`).
2. Kept set = N best by (e-value class with 1e-180 equality, first-HSP score in the e-value-sorted list, later oid first) (`blast_hits.c:3078-3107`, `3243-3299`).
3. `hsp_array[0]` state: score order unless the query's hit list overflowed, then e-value order, for the subject order of the kept lists as well
   (`link_hsps.c:1803`, `blast_hits.c:3269-3285`, `blast_traceback.c:1770-1772`).
4. stderr `Warning: [tblastx] Examining 5 or more matches is recommended` for N < 5 in every format (`blast_args.cpp:2975-2977`).
5. outfmt 0: description rows up to N (500 default) but alignment blocks only up to 250 unless `-max_target_seqs` is given (`blast_args.cpp:2913-2928`, `blast_format.cpp:1551`).

## 7. Open points / UNSURE
* `-culling_limit` (hit list cut moves to the traceback stage; final subject sort sees the score state): two oracle runs only (row `CreateHspWriter / CreateHspPipe`).
* `PruneSeqalign` at the formatter counts runs of the same subject Seq-id; two adjacent subject records with the same id near the N-th place were not tried.
* Behaviour with `BL2SEQ_LEGACY` (hit lists would then be per (query, subject) sequence comparison) is out of scope (environment unset).
* Not tried: ties through `-num_threads` (ignored with `-subject`), `-max_target_seqs` together with `-seg`/`-window_size` (they change e-values/HSPs, not the hit list rule).
