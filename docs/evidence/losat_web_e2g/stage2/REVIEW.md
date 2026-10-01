# E2g stage 2: orchestrator review of the agents' findings

Each agent's divergent / unported / rejected / UNSURE rows, checked against the NCBI source and LOSAT by the orchestrator (Opus). Decision: `transpile` (port in the E2g transpile step), `keep-rejected` (with reason), `no-change` (agent was wrong, with reason).

## E1 (engine and gapped alignment): 68 rows; faithful 52, n/a 11, divergent 3, unported 1, rejected 1

| # | Row | Check | Decision |
|---|---|---|---|
| E1-1 | init hit order (`score_compare_match`, blast_extend.c:274-310; sorted at na_ungapped.c:1690-1691) | Confirmed in part: LOSAT `run.rs` `score_compare_ungapped_hits` compares `qs`, which is context-local in at least one producer (`run.rs` ~8439 `ungapped.q_start - q_context_start`), whereas NCBI's `ungapped_data->q_start` is in the concatenated query; and `sort_unstable_by` where glibc 2.39 `qsort` is a stable merge sort. Two hits on different strands with equal score, subject start and length order differently (NCBI: plus strand first). | transpile: compare the concatenated query offset, stable sort |
| E1-2 | `BlastGetStartForGappedAlignmentNucl` (blast_gapalign.c:3323-3389) | Confirmed: `gapped.rs` `blast_get_start_for_gapped_alignment_nucl` computes `(s_gapped_start - s_offset).min(q_gapped_start - q_offset)` in `usize`; NCBI's `Int4 offset = MIN(...)` can be negative (start left of the HSP start). Release wraps and picks the other operand; debug panics. | transpile: Int4 arithmetic for the whole function |
| E1-3 | `s_Blast_HSPListReapByPrelimEvalue` `evalue > cutoff` vs LOSAT `e <= threshold` | Differ only for NaN; LOSAT rejects NaN and infinite `-evalue` (E2c §G), so unreachable. | transpile anyway (same comparison as C, no output change) |
| E1-4 | `BL2SEQ_LEGACY` (`db_length == 0` path, `BLAST_OneSubjectUpdateParameters`) | The variable also changes the C++ formatting (bl2seq legacy report; `run_oracle.py` refuses to run with it). Porting means the legacy report layer. | keep-rejected: add an explicit rejection of the environment variable, as for `BATCH_SIZE`/`CHUNK_SIZE` (check with ranges A and F which other variables change blastn output) |
| E1-5 | greedy gap costs above 32767 | Existing explicit rejection (E2c §G). | keep-rejected |

## E2 (traceback and hit saving): 119 rows; faithful 85, n/a 33, divergent 1

| # | Row | Check | Decision |
|---|---|---|---|
| E2-1 | `Blast_TracebackFromHSPList` HSP order (blast_traceback.c:358-365) | Confirmed: NCBI sorts by score only under `#ifdef _DEBUG` (the 2.17.0 oracle is a release build; `ASSERT` is empty), so a list that `Blast_HitListUpdate` heapified (more than `prelim_hitlist_size` subjects) is traced in e-value order (`Blast_HSPListSortByEvalue`), others in the order of the preliminary stage. LOSAT `run.rs` ~10530 `sort_prelim_hits_by_score` re-sorts every list by score. Differs when e-value order is not score order inside one list (contexts with different gapped blocks: gap costs above the table). | transpile: trace in the stored list order. NCBI's stored order is score order for a list that was not heapified (`Blast_HSPListSortByScore` per subject chunk at blast_engine.c:555, then `s_BlastHSPListsCombineByScore` merges and order-preserving reaps) and e-value order for a heapified one; LOSAT's re-sort by score is a no-op for the first, so drop it and keep the heap's e-value order |
| E2-note | one interval tree for all queries of a subject (NCBI: one per query list) | Agent argues equivalence (contexts never overlap); a structural effect needs two HSPs with the same endpoint. | check during the transpile of E2-1 (per-query trees if cheap) |

## D (lookup, scan, ungapped): 88 rows; faithful 64, n/a 17, divergent 7

| # | Row | Check | Decision |
|---|---|---|---|
| D-1..4 | small and standard lookup table cell order (`BlastLookupAddWordHit` appends, blast_lookup.c:60-75; `BlastLookupIndexQueryExactMatches` walks the query forward; `s_BlastSmallNaLookupFinalize` copies in order) | Confirmed: NCBI's eSmallNaLookupTable/eNaLookupTable cells list query offsets ascending; NCBI's eMBLookupTable chain is newest first (`next_pos[index] = hashtable[ecode]; hashtable[ecode] = index`, blast_nalookup.c:1090-1091), and LOSAT stores the small table (lut < word) in the megablast chain, so its small-table seeds at one subject offset come in reverse order. Can matter with the diagonal hash (query block > 8000: default blastn with a query above ~4 kb), where `s_BlastDiagHashInsert` reuses a stale cell of another diagonal. | transpile: ascending cell order for the small and standard tables |
| D-5 | diagonal container (blast_parameters.c:225-231) | Confirmed: NCBI uses the array when the query block (last context offset + length) is at most 8000 (`kQueryLenForHashTable`), whatever the number of queries; LOSAT `run.rs` ~7239 also requires `queries.len() == 1` (E2c §J left it as "almost the same"). | transpile: array whenever NCBI uses it |
| D-6, D-7 | init hit sort (= E1-1) | see E1-1 | transpile |
| D-scan | `s_MBChooseScanSubject` never selects 9_1 for lut 9 (missing `else`); LOSAT selects 9_1 | Agent: both return the same words. | no-change (equivalent; note in INVENTORY) |

## B (C++ API layer): 175 rows; faithful 127, n/a 41, rejected 5, divergent 1, unported 1

| # | Row | Check | Decision |
|---|---|---|---|
| B-1 | `SplitQuery_GetOverlapChunkSize` with `OVERLAP_CHUNK_SIZE` (split_query_aux_priv.cpp:51-61) | Confirmed: NCBI reads `OVERLAP_CHUNK_SIZE` with `NStr::StringToInt` (when not blank); LOSAT keeps 100 and does not reject the variable (`run.rs` ~5462 checks only `BATCH_SIZE` and `CHUNK_SIZE`). | transpile, together with `BATCH_SIZE` (blast_input_aux.cpp:85-91) and `CHUNK_SIZE` (local_blast.cpp:54-62): read the three variables as NCBI does now that batches and splits are ported (removes two E2c §G rejections); values `NStr::StringToInt` does not accept, and values whose NCBI behaviour is not reproduced, stay explicit rejections |
| B-2 | `CMultiSeqInfo` / `CLocalDbAdapter(..., false)` with `BL2SEQ_LEGACY` (blast_app_util.cpp:206-211) | Confirmed: the variable turns off dbscan mode, which also changes `CLocalBlast` (local_blast.cpp:189, 289) and per-subject statistics (= E1-4). A separate legacy search mode. | keep-rejected: new explicit rejection of `BL2SEQ_LEGACY` |
| B-rej | CHUNK_SIZE, K-A failure with invalid queries (2 branches), subject without residues, subject total ≥ 2^31, other tasks | Existing E2c/E2f rejections. CHUNK_SIZE: see B-1. The K-A branch "first batch has only invalid queries" is portable (NCBI writes that batch, then errors on the next batch with valid queries); the mixed branch (NCBI crashes) stays. | B-1; port the first K-A branch; others keep-rejected |

## A (application, arguments, input): 236 rows; faithful 118, n/a 62, rejected 42, divergent 6, unported 8

| # | Row | Check | Decision |
|---|---|---|---|
| A-1 | argument errors (NCBI USAGE + exit 1, clap exit 2), `-help` text, output write failure (exit 6 / abort vs 1), non-UTF-8 `-subject` path in the `Database:` tag | Known and common to every program's CLI; the plan already moves them to S08+ (E2c §I, §N: approved exception or NCBI text/exit, a maintainer decision). | deferred to S08+ (maintainer decision pending; not a BLASTN engine path) |
| A-2 | no `-subject` (`CBlastDatabaseArgs::ExtractAlgorithmOptions`: "Either a BLAST database or subject sequence(s) must be specified", exit 1) | LOSAT's clap makes `-subject` required (exit 2). Specific to the BLASTN option set and small. | transpile (check the exact NCBI text and when it is raised) |
| A-3 | warning timing per query batch (`CFastaReader` title warnings when the batch is read; invalid-query warnings with the batch's results) | With two or more batches the order of the stderr lines differs (title warnings of all queries first in LOSAT). | transpile: write the warnings per batch in NCBI's order (with per-batch output, see the K-A item) |
| A-4 | `-num_threads` above the CPU count and with `-subject` (warnings) | Known E2c §I difference, common to all programs, recorded as an open plan item. | deferred to S08+ (maintainer decision pending) |
| A-5 | `CTOOLKIT_COMPATIBLE` (showdefline.cpp `kBits`: "(bits)" instead of "(Bits)" in the outfmt 0 description header; the other uses are under `#ifdef CTOOLKIT_COMPATIBLE`, not defined in the oracle build) | Oracle-confirmed by the agent. | transpile |
| A-6 | `PRE_FETCH_SEQS_LIMIT` (an integer only toggles prefetching; a non-integer is a CStringException, exit 255, before results) | Oracle-confirmed by the agent. | transpile the error (check the text) |
| A-7 | `BL2SEQ_LEGACY` | = E1-4, B-2 | keep-rejected (new explicit rejection) |
| A-8 | `OVERLAP_CHUNK_SIZE` (a non-integer exits 255 even without a split) | = B-1 | transpile with B-1 (also the error of a non-integer value) |
| A-9 | out of memory ("BLAST ran out of memory", exit 4) | A Rust allocation failure aborts; reproducing it would need fallible allocation everywhere. | deferred (maintainer decision: accepted limitation) |
| A-10 | NCBI toolkit diagnostics variables (`DIAG_POST_LEVEL` hides the `Warning: [blastn]` lines; `Trace` adds Info lines) | The CNcbiDiag layer. | keep-rejected: new explicit rejection of the diagnostics variables that change blastn's stderr (list them from ncbidiag.cpp) |
