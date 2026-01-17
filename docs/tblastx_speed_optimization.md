# TBLASTX Single-Thread Speed Optimization Findings

Scope: tblastx runtime hot paths with emphasis on single-thread speed. Output must remain bit-for-bit identical to NCBI BLAST+.

## Constraints
- No output changes: bit-perfect parity required.
- Prefer single-thread speedups over parallelization.
- Preserve NCBI timing/order for algorithmic steps (scan -> extend -> reevaluate -> link -> output).

## Hot Path Summary (where CPU time concentrates)
- Subject scan + two-hit gating: `LOSAT/src/algorithm/tblastx/blast_aascan.rs` and inner loop in `LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs:559`
- Two-hit ungapped extension: `LOSAT/src/algorithm/tblastx/extension/two_hit.rs`
- Ungapped reevaluation and identity counting (SIMD already used): `LOSAT/src/algorithm/tblastx/reevaluate.rs`
- Sum-stats linking with multiple sorts: `LOSAT/src/algorithm/tblastx/sum_stats_linking/linking.rs`
- Translation and frame generation per subject: `LOSAT/src/algorithm/tblastx/translation.rs`

## Quick Wins (Low Risk, High Impact)
- Gate or remove unconditional debug eprintlns in the hot path; they execute for every subject and dominate runtime on large inputs. `LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs:1256`
- Avoid per-hit `String` clones by storing `query_id`/`subject_id` as indices (or `Arc<str>`), and resolve only at final output. `LOSAT/src/common.rs:79` `LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs:1225` `LOSAT/src/algorithm/tblastx/blast_engine/run_neighbor_map.rs:840`
- Skip redundant E-value filtering in writer thread (already filtered before send); keeps output unchanged while saving a pass. `LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs:288`
- Add a single-thread fast path to bypass Rayon/channel/progress-bar overhead (both main run and linking). `LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs:350` `LOSAT/src/algorithm/tblastx/sum_stats_linking/linking.rs:361`
- Reuse per-subject vectors instead of re-allocating each subject: `cutoff_scores`, `length_adj_per_context`, `eff_searchsp_per_context`, `combined_ungapped_hits`, `init_hsps`. `LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs:404`

## Medium/High Effort (Likely Biggest Single-Thread Gains)
- Diagonal array reset is O(diag_array_size) per subject; consider generation-stamp or lazy-init to avoid full clears while preserving NCBI semantics. `LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs:360`
- Scan loop micro-optimizations in `s_blast_aa_scan_subject`: special-case the common unmasked range to avoid `determine_scanning_offsets`, and use wordsize=3 constants to reduce shifts/masks. `LOSAT/src/algorithm/tblastx/blast_aascan.rs:56`
- Two-hit extension currently uses scalar left/right loops; consider reusing SIMD-aware helpers from `extension/ungapped.rs` to reduce per-hit cost while keeping NCBI control flow. `LOSAT/src/algorithm/tblastx/extension/two_hit.rs:10`
- Cache CPU feature detection (`is_x86_feature_detected!`) once and reuse in hot helpers to cut per-call overhead for identity and reevaluation. `LOSAT/src/algorithm/tblastx/reevaluate.rs:481`
- Reduce per-subject recomputation of `gap_trigger` (context-only) by caching per context outside the subject loop. `LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs:422`
- Translation allocates reverse complements per subject; consider streaming reverse translation or reuse buffers to reduce allocations. `LOSAT/src/algorithm/tblastx/translation.rs:98`

## Neighbor-Map Mode Specific Notes
- `all_ungapped` uses a `Mutex` even though subject loop is sequential; use a plain `Vec` to avoid lock overhead. `LOSAT/src/algorithm/tblastx/blast_engine/run_neighbor_map.rs:192`
- Per-frame diag_array is `Vec<Vec<DiagStruct>>` built every frame; flatten and reuse to reduce allocation and improve cache locality. `LOSAT/src/algorithm/tblastx/blast_engine/run_neighbor_map.rs:348`
- `ctx_base` mapping is rebuilt inside each frame loop; precompute once outside the parallel frame loop. `LOSAT/src/algorithm/tblastx/blast_engine/run_neighbor_map.rs:338`
- Add a single-thread path for frame processing and linking to avoid Rayon overhead when `num_threads==1`. `LOSAT/src/algorithm/tblastx/blast_engine/run_neighbor_map.rs:306` `LOSAT/src/algorithm/tblastx/blast_engine/run_neighbor_map.rs:673`

## Notes for Parity
- Any change must preserve NCBI control flow and exact scoring/threshold logic.
- Do not change sort orders or tie-breaking in linking/output.
- Gate debug logging via existing env flags (`LOSAT_*`) rather than removing if needed.
