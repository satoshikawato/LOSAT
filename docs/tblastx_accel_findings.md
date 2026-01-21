# TBLASTX Acceleration Findings (Single-Thread, Parity-Safe)

Scope: targeted micro-optimizations that preserve exact NCBI control flow,
tie-breaking, and floating-point behavior. Findings are based on the current
LOSAT codebase and do not include profiling results.

## SIMD Opportunities

- Two-hit extension is scalar in `LOSAT/src/algorithm/tblastx/extension/two_hit.rs`.
  Reuse the SIMD-aware left/right helpers from
  `LOSAT/src/algorithm/tblastx/extension/ungapped.rs` so the per-residue
  accumulation and X-drop checks remain strictly sequential.
- Optional SSE2 fallback for ungapped extension (if AVX2 is not available).
  This mirrors the identity count approach in
  `LOSAT/src/algorithm/tblastx/reevaluate.rs`, but must preserve the same
  loop order and stop conditions.
- Keep gather-based SIMD limited to score lookup only. Any SIMD that changes
  branch timing, early-termination, or accumulation order risks parity.

## Reallocation and Memory Reuse

- Pre-allocate output vectors:
  - `kept_hits` in `LOSAT/src/algorithm/tblastx/blast_engine/mod.rs` can
    reserve `ungapped_hits.len()`.
  - `ungapped_hits` in `LOSAT/src/algorithm/tblastx/blast_gapalign.rs` can
    reserve `init_hsps.len()`.
- Hoist per-run data that does not depend on subject:
  - `context_params` in `LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs`.
  - `x_dropoff_per_context` in
    `LOSAT/src/algorithm/tblastx/blast_engine/run_neighbor_map.rs`.
- Reuse per-frame scratch in neighbor-map mode:
  - Build `ctx_base` once and reuse.
  - Build `temp_contexts` once and reuse.
  - Allocate `diag_array` in `map_init` so each thread reuses a single buffer.
- Consider changing `get_ungapped_hsp_list` to accept a mutable buffer and
  clear it in place (affects call sites in
  `LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs` and
  `LOSAT/src/algorithm/tblastx/blast_engine/run_neighbor_map.rs`).

## Parity Constraints (Non-Negotiable)

- Preserve per-residue iteration order, branch timing, and X-drop logic.
- Avoid fast-math or FMA changes; do not alter floating-point precision.
- Keep sort keys and stability identical to NCBI equivalents.
- Ensure any SIMD only replaces scalar score lookup, not control flow.
