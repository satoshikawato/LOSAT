# Shared Single-Thread Speed Opportunities (blastn + tblastx)

Scope: architectural/components/algorithms shared by blastn and tblastx where single-thread speed can improve without changing output.

## Constraints
- Output must remain bit-for-bit identical to NCBI BLAST+.
- Do not change algorithm timing/order; focus on allocation, I/O, and dispatch overhead.

## High-Impact Shared Targets
- Output ordering/grouping uses general-purpose hash containers and extra passes; replace with indexed vectors keyed by `q_idx`/`s_idx` or switch to `FxHashMap`/`FxHashSet` to reduce hashing and allocation. `LOSAT/src/common.rs`
- Per-hit `String` cloning of IDs is common in both pipelines; store IDs by index (or `Arc<str>`) in `Hit` and resolve only during output. `LOSAT/src/common.rs` `LOSAT/src/algorithm/blastn/blast_engine/run.rs` `LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs` `LOSAT/src/algorithm/tblastx/blast_engine/run_neighbor_map.rs`
- Tabular output builds a new `String` per hit; add a writer-based formatter or reuse a scratch buffer to cut per-hit allocation. `LOSAT/src/common.rs` `LOSAT/src/report/outfmt6.rs`
- Single-thread fast path: when `num_threads == 1`, bypass Rayon, channels, and progress bars, then write output directly. `LOSAT/src/algorithm/blastn/blast_engine/run.rs` `LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs` `LOSAT/src/algorithm/tblastx/blast_engine/run_neighbor_map.rs`
- Progress bar and debug logging overhead: ensure all hot-loop logging is strictly gated by `verbose`/debug env vars. `LOSAT/src/algorithm/blastn/blast_engine/run.rs` `LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs`

## Medium-Effort Shared Targets
- Pre-reserve or reuse per-subject/per-query vectors (hit collection, grouping, subject lists) to reduce realloc churn. `LOSAT/src/common.rs` `LOSAT/src/algorithm/blastn/blast_engine/run.rs` `LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs`
- Cache CPU feature detection (`is_x86_feature_detected!`) in a static once-cell to avoid per-call overhead in identity/reevaluation helpers used by both pipelines. `LOSAT/src/algorithm/tblastx/reevaluate.rs` `LOSAT/src/algorithm/blastn/extension.rs`

## Notes for Parity
- Do not alter sorting keys or tie-breaking in output ordering.
- Any optimization that changes data layout or buffering must preserve exact output ordering and formatting.
