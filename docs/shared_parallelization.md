# Shared Parallelization Opportunities (blastn + tblastx)

Scope: shared components/architectural steps where parallelization can improve throughput without changing output.

## Constraints
- Output must remain bit-for-bit identical to NCBI BLAST+.
- Preserve deterministic ordering and tie-breaking.
- Keep output formatting identical (same field ordering and precision).

## Shared Targets
- Output ordering/grouping: build per-query/per-subject groups in parallel, then merge in deterministic `q_idx` order before `write_output_ncbi_order` emits results. `LOSAT/src/common.rs`
- Output formatting: compute formatted lines in parallel after ordering is final, then write sequentially to preserve order. `LOSAT/src/common.rs` `LOSAT/src/report/outfmt6.rs`
- Post-processing filters (evalue/identity/mismatch) can be parallelized once HSP ordering is no longer needed; re-sort or preserve original indices. `LOSAT/src/common.rs` `LOSAT/src/algorithm/tblastx/reevaluate.rs`

## Implementation Sketches (NCBI Refs)
- Parallel group build, deterministic merge:
  - Step 1: `par_iter()` over hits to bucket by `(q_idx, s_idx)` using thread-local maps/vectors.
  - Step 2: Merge buckets in deterministic `q_idx` then `s_idx` order.
  - Step 3: Apply NCBI comparators for subject/HSP ordering (NCBI: `blast_hits.c` s_EvalueCompareHSPLists, ScoreCompareHSPs).
- Parallel formatting after ordering:
  - Step 1: Build a fully ordered hit list using `write_output_ncbi_order` semantics.
  - Step 2: `par_iter()` format each line using the same precision rules as outfmt6/outfmt7, then write sequentially.
  - NCBI refs: output ordering matches `BLAST_LinkHsps()` + `s_EvalueCompareHSPLists()` + `ScoreCompareHSPs()` in `blast_hits.c`.
- Parallel post-filter computations:
  - Step 1: After ordering, `par_iter()` compute identity/length/mismatch or e-value gate checks.
  - Step 2: Preserve stable indices and reassemble in original order before output.
  - NCBI refs: identity calculation in `blast_hits.c:746-792`, E-value comparison in `blast_hits.c:1390-1403`.

## Shared Cautions
- Do not parallelize stages that depend on ordered, stateful scans (diagonal tables, interval trees).
- Avoid parallel map iteration that changes tie ordering; re-sort with explicit keys.
