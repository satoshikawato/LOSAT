# LOSATN Single-Thread Speed Optimization Findings

Scope: blastn (losatn) runtime hot paths with emphasis on single-thread speed. Output must remain bit-for-bit identical to NCBI BLAST+.

## Constraints
- No output changes: bit-perfect parity required.
- Prefer single-thread speedups over parallelization.
- Keep NCBI behavior and timing/order intact.

## Quick Wins (Low Risk, High Impact)
- Remove or gate hot-loop debug timing and logging in two-stage scan. These run per-kmer and cost real time even when not debugging. `LOSAT/src/algorithm/blastn/blast_engine/run.rs:1317` `LOSAT/src/algorithm/blastn/blast_engine/run.rs:1322` `LOSAT/src/algorithm/blastn/blast_engine/run.rs:1989`
- Avoid per-hit `String` clones by storing `query_id`/`subject_id` as indices or `Arc<str>` and resolving at output. Reduces allocation churn for large hit sets. `LOSAT/src/common.rs:79` `LOSAT/src/algorithm/blastn/blast_engine/run.rs:2906`
- Group endpoint purge by subject index (`s_idx`) instead of `subject_id` string to avoid hashing/cloning. `LOSAT/src/algorithm/blastn/filtering/purge_endpoints.rs:360`
- Reuse/pre-allocate per-subject vectors (hits, ungapped hits, cutoff arrays, diag structures) to reduce allocator pressure. `LOSAT/src/algorithm/blastn/blast_engine/run.rs:1059` `LOSAT/src/algorithm/blastn/blast_engine/run.rs:1142` `LOSAT/src/algorithm/blastn/blast_engine/run.rs:1268`
- For single-thread runs (`num_threads == 1`), bypass Rayon channel/progress bar overhead and run a sequential loop. `LOSAT/src/algorithm/blastn/blast_engine/run.rs:1039` `LOSAT/src/algorithm/blastn/blast_engine/run.rs:3230`

## Medium/High Effort (Potential Largest CPU Savings)
- Optimize k-mer scan to reduce per-base div/mod and branching. Consider rolling byte-level access or SIMD for packed data; keep exact NCBI scan semantics. `LOSAT/src/algorithm/blastn/blast_engine/run.rs:163` `LOSAT/src/algorithm/blastn/blast_engine/run.rs:189`
- Ungapped extension loops are SIMD-friendly on 4-base blocks; reduce bounds checks with safe preconditions or use vectorized loads plus LUT. `LOSAT/src/algorithm/blastn/extension.rs:88` `LOSAT/src/algorithm/blastn/extension.rs:227`
- Pre-size and reuse gapped alignment scratch buffers (trace rows, DP arrays) more aggressively to reduce realloc during long alignments. `LOSAT/src/algorithm/blastn/alignment/gapped.rs:112`
- Avoid rebuilding the blastna scoring matrix during endpoint purge; pass or cache the prebuilt matrix from the run context. `LOSAT/src/algorithm/blastn/filtering/purge_endpoints.rs:1030` `LOSAT/src/algorithm/blastn/blast_engine/run.rs:1006`

## Suggested Order of Work
1) Gate/remove hot-loop debug/perf logging.
2) Single-thread fast path without Rayon/channels.
3) Reduce allocations for hits and subject grouping (IDs by index, group by `s_idx`).
4) Pre-allocation of per-subject buffers and scratch.
5) K-mer scan + ungapped extension micro-optimizations (SIMD/bit tricks).

## Notes for Parity
- All changes must preserve ordering and exact scoring/evalue logic.
- Any algorithmic change must be validated against NCBI source; avoid speculative optimizations.
