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
| E2-1 | `Blast_TracebackFromHSPList` HSP order (blast_traceback.c:358-365) | Confirmed: NCBI sorts by score only under `#ifdef _DEBUG` (the 2.17.0 oracle is a release build; `ASSERT` is empty), so a list that `Blast_HitListUpdate` heapified (more than `prelim_hitlist_size` subjects) is traced in e-value order (`Blast_HSPListSortByEvalue`), others in the order of the preliminary stage. LOSAT `run.rs` ~10530 `sort_prelim_hits_by_score` re-sorts every list by score. Differs when e-value order is not score order inside one list (contexts with different gapped blocks: gap costs above the table). | transpile: trace in the stored list order (check that LOSAT's stored order of non-heapified lists is NCBI's) |
| E2-note | one interval tree for all queries of a subject (NCBI: one per query list) | Agent argues equivalence (contexts never overlap); a structural effect needs two HSPs with the same endpoint. | check during the transpile of E2-1 (per-query trees if cheap) |

## D (lookup, scan, ungapped): 88 rows; faithful 64, n/a 17, divergent 7

| # | Row | Check | Decision |
|---|---|---|---|
| D-1..4 | small and standard lookup table cell order (`BlastLookupAddWordHit` appends, blast_lookup.c:60-75; `BlastLookupIndexQueryExactMatches` walks the query forward; `s_BlastSmallNaLookupFinalize` copies in order) | Confirmed: NCBI's eSmallNaLookupTable/eNaLookupTable cells list query offsets ascending; NCBI's eMBLookupTable chain is newest first (`next_pos[index] = hashtable[ecode]; hashtable[ecode] = index`, blast_nalookup.c:1090-1091), and LOSAT stores the small table (lut < word) in the megablast chain, so its small-table seeds at one subject offset come in reverse order. Can matter with the diagonal hash (query block > 8000: default blastn with a query above ~4 kb), where `s_BlastDiagHashInsert` reuses a stale cell of another diagonal. | transpile: ascending cell order for the small and standard tables |
| D-5 | diagonal container (blast_parameters.c:225-231) | Confirmed: NCBI uses the array when the query block (last context offset + length) is at most 8000 (`kQueryLenForHashTable`), whatever the number of queries; LOSAT `run.rs` ~7239 also requires `queries.len() == 1` (E2c §J left it as "almost the same"). | transpile: array whenever NCBI uses it |
| D-6, D-7 | init hit sort (= E1-1) | see E1-1 | transpile |
| D-scan | `s_MBChooseScanSubject` never selects 9_1 for lut 9 (missing `else`); LOSAT selects 9_1 | Agent: both return the same words. | no-change (equivalent; note in INVENTORY) |
