# LOSATN (blastn) Parallelization Opportunities

Scope: parallelism opportunities in blastn while keeping bit-perfect NCBI output.

## Constraints
- Output must be bit-for-bit identical to NCBI BLAST+.
- Do not change algorithm order/timing in a way that changes results.
- Determinism matters: keep stable ordering and explicit tie-breakers.

## Current Parallelization
- Subject-level parallelism via Rayon: `subjects.par_iter().for_each_init(...)`. `LOSAT/src/algorithm/blastn/blast_engine/run.rs`
- Writer thread aggregates hits and writes output. `LOSAT/src/algorithm/blastn/blast_engine/run.rs`

## Opportunities (Low/Medium Risk)
- Parallelize gapped extension after ungapped hit collection: compute gapped alignments in parallel with per-thread `GapAlignScratch`, then sort by score/tie-breakers and run containment check sequentially in the same order as today. `LOSAT/src/algorithm/blastn/blast_engine/run.rs`
- Parallelize post-processing in the writer thread: apply `max_hsps_per_subject` and `hitlist_size` per `(q_idx, s_idx)` group in parallel, then merge deterministically before output ordering. `LOSAT/src/algorithm/blastn/blast_engine/run.rs`
- Parallelize reevaluation of trimmed HSPs (Phase 2 step 3) inside a subject: each HSP re-eval is independent; collect results, then re-sort and continue with the same filtering order. `LOSAT/src/algorithm/blastn/blast_engine/run.rs`
- Parallelize formatting of tabular output lines (outfmt 6/7) after final ordering; write sequentially. `LOSAT/src/common.rs` `LOSAT/src/report/outfmt6.rs`

## Implementation Sketches (NCBI Refs)
- Parallel gapped extension, sequential containment:
  - Step 1: Collect and sort `UngappedHit` by score DESC (NCBI: `blast_gapalign.c:3824` init_hsp_array sorted by score).
  - Step 2: `par_iter()` over sorted hits, compute gapped alignments using per-thread `GapAlignScratch` (NCBI: `blast_gapalign.c:313-319`).
  - Step 3: Sequential pass over results in the *original score order* to apply containment via interval tree (NCBI: `blast_gapalign.c:3811-3831`, `blast_traceback.c:679-692`).
- Parallel Phase 2 reevaluation:
  - Step 1: After purge pass 1 (NCBI: `blast_traceback.c:637-638`), split `local_hits[extra_start..]`.
  - Step 2: `par_iter()` each hit through `reevaluate_hsp_with_ambiguities_gapped_ex` (NCBI: `blast_traceback.c:647-665`).
  - Step 3: Collect survivors and re-sort by ScoreCompareHSPs before containment/filtering (NCBI: `blast_hits.c` ScoreCompareHSPs comparator).
- Parallel hitlist pruning in writer thread:
  - Step 1: Group hits by `(q_idx, s_idx)` and per-subject truncate to `max_hsps_per_subject` in parallel (NCBI: `blast_hits.c:2049-2067`).
  - Step 2: For each query, compute best_evalue/best_score and sort subjects using s_EvalueCompareHSPLists (NCBI: `blast_hits.c:3078-3095`).
  - Step 3: Merge groups in deterministic query order, then output.
- Parallel outfmt formatting after final ordering:
  - Step 1: Use `write_output_ncbi_order` to produce ordered hit list (NCBI ordering: `blast_hits.c` ScoreCompareHSPs + s_EvalueCompareHSPLists).
  - Step 2: `par_iter()` format lines (no reordering), then write sequentially to preserve order.

## High-Risk / Likely Not Safe
- Split the subject scan across chunks: the two-hit diagonal table is order-dependent and shared across the whole subject; chunking can alter hit detection. `LOSAT/src/algorithm/blastn/blast_engine/run.rs`
- Parallelize the interval-tree containment pass: containment depends on processing hits in strict score order. `LOSAT/src/algorithm/blastn/blast_engine/run.rs`

## Notes for Parity
- Any parallel step must preserve the same sort keys and tie-breakers as NCBI.
- If parallel work changes the order of equal-score HSPs, output can diverge; re-sort with explicit keys after parallel steps.
