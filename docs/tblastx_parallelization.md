# TBLASTX Parallelization Opportunities

Scope: parallelism opportunities in tblastx while keeping bit-perfect NCBI output.

## Constraints
- Output must be bit-for-bit identical to NCBI BLAST+.
- Preserve NCBI ordering and tie-breaking; avoid nondeterministic merges.
- Avoid changing two-hit/diag-table semantics.

## Current Parallelization
- Subject-level parallelism in backbone mode: `subjects_raw.par_iter().for_each_init(...)`. `LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs`
- Sum-stats linking parallelizes frame groups with order restoration. `LOSAT/src/algorithm/tblastx/sum_stats_linking/linking.rs`
- Neighbor-map mode parallelizes subject frames and per-subject linking. `LOSAT/src/algorithm/tblastx/blast_engine/run_neighbor_map.rs`

## Opportunities (Low/Medium Risk)
- Neighbor-map mode: parallelize subjects (currently sequential) and collect per-subject hit lists; remove the global `Mutex<Vec<_>>` and merge deterministically. `LOSAT/src/algorithm/tblastx/blast_engine/run_neighbor_map.rs`
- Parallelize per-subject identity/bit-score computation across final hits (pure computations), then keep output ordering via `write_output_ncbi_order`. `LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs`
- Parallelize HSP reevaluation within a subject (each HSP is independent), then re-sort as done today. `LOSAT/src/algorithm/tblastx/reevaluate.rs` `LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs`
- Parallelize output formatting (outfmt 6/7) after final ordering; write sequentially. `LOSAT/src/common.rs` `LOSAT/src/report/outfmt6.rs`

## Implementation Sketches (NCBI Refs)
- Neighbor-map mode subject parallelism:
  - Step 1: Convert subject loop to `par_iter()` and return per-subject `Vec<UngappedHit>` (no global `Mutex`). 
  - Step 2: After join, group by `s_idx` and run per-subject linking (NCBI: `link_hsps.c:553-983` per-subject cutoffs and linking).
  - Step 3: Ensure deterministic ordering by sorting final hits with NCBI comparators (NCBI: `link_hsps.c:990-1000`, `blast_hits.c` ScoreCompareHSPs).
- Parallel reevaluation within subject:
  - Step 1: After merge of per-frame init HSPs, `par_iter()` reevaluate each HSP (`Blast_HSPReevaluateWithAmbiguitiesUngapped` equivalent). NCBI: `blast_hits.c:675-733`.
  - Step 2: Collect survivors, then keep the same ordering/sorting as current pipeline before linking.
- Parallel identity/bit-score computation:
  - Step 1: After linking, build an index list of hits in final order.
  - Step 2: `par_iter()` compute identity (NCBI: `blast_hits.c:746-792`) and bit-score (NCBI: `blast_hits.c:1833-1918`).
  - Step 3: Write using `write_output_ncbi_order` to preserve ordering.
- Parallel outfmt formatting:
  - Step 1: Use `write_output_ncbi_order` to produce ordered hits (NCBI ordering: `link_hsps.c` + `blast_hits.c` comparators).
  - Step 2: Format lines in parallel, then write sequentially.

## High-Risk / Likely Not Safe
- Parallelize subject frames in backbone mode: the diagonal table is shared across frames and order-dependent. `LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs`
- Parallelize the two-hit scan within a frame by chunking subject ranges; diagonal state and hit gating are order-dependent. `LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs` `LOSAT/src/algorithm/tblastx/blast_aascan.rs`

## Notes for Parity
- Use explicit, deterministic sort keys after any parallel step.
- Do not parallelize steps that rely on monotonic updates to the diagonal table or interval/tree structures.
