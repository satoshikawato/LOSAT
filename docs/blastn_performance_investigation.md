# BLASTN Performance Investigation (LOSAT vs NCBI BLAST+)

## Context
- Benchmark summary provided:
  - Megablast EDL933_vs_Sakai: BLAST+ 0.998s vs LOSAT 3.359s (3.37x)
  - Megablast NZ_CP006932_Self: BLAST+ 0.229s vs LOSAT 0.685s (2.99x)
  - Megablast Sakai_vs_MG1655: BLAST+ 0.888s vs LOSAT 3.704s (4.17x)
- Scope of review:
  - LOSAT `blastn` pipeline (`LOSAT/src/algorithm/blastn/blast_engine/run.rs`, `LOSAT/src/algorithm/blastn/coordination.rs`, `LOSAT/src/report/outfmt6.rs`, `LOSAT/src/common.rs`).
  - NCBI BLAST+ core references under `/mnt/c/Users/genom/GitHub/ncbi-blast/` (notably `blast_engine.c`, `blast_nascan.c`, `blast_nalookup.h`, `lookup_wrap.c`, `blast_hits.h`, `blast_traceback.c`, `blast_gapalign.c`).

## Findings (Likely Performance Drivers)

### 1) Hit storage and pruning are heavier than NCBI
- LOSAT stores `query_id` and `subject_id` as `Arc<str>` in every `Hit` (per-hit atomic refcount work and larger hit objects).
  - LOSAT: `LOSAT/src/common.rs:92-128`.
  - NCBI: `BlastHSPList` stores only `oid` and `query_index` (no per-hit string ownership) in `ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-180`.
- LOSAT defers `max_hsps_per_subject` and `hitlist_size` pruning to a global post-pass using `FxHashMap` grouping.
  - LOSAT: `LOSAT/src/algorithm/blastn/blast_engine/run.rs:3047-3174`.
  - NCBI: per-query/per-subject pruning in traceback phase (`blast_traceback.c:850-892`) and `Blast_TrimHSPListByMaxHsps` (`blast_hits.c:2049-2067`).
- Impact: higher memory pressure, more hashing/sorting, and extra per-hit work compared to NCBI's in-structure pruning model.

### 2) Offset-pair buffers are reallocated per chunk
- LOSAT allocates a new `Vec<OffsetPair>` per subject chunk, then drains it.
  - LOSAT: `LOSAT/src/algorithm/blastn/blast_engine/run.rs:4219-4230`.
- NCBI sizes `offset_pairs` once using `GetOffsetArraySize` and reuses it across scans.
  - NCBI: `lookup_wrap.c:255-288`.
- Impact: repeated allocations and deallocations in the hottest scan/extend loop.

### 3) Scan dispatch is more dynamic than NCBI
- LOSAT chooses scanning behavior via `match` in the scan function for each invocation, using `FnMut` callbacks.
  - LOSAT: `LOSAT/src/algorithm/blastn/blast_engine/run.rs:1900-2230`.
- NCBI selects the scan routine once via `s_MBChooseScanSubject` and uses a function pointer for repeated calls.
  - NCBI: `blast_nascan.c:2602-2677`.
- Impact: extra branching and closure overhead vs. NCBI's static dispatch model.

### 4) Diagonal table clearing is potentially more expensive
- LOSAT clears `hit_level_array`/`hit_len_array` with `resize` + `fill` per chunk.
  - LOSAT: `LOSAT/src/algorithm/blastn/blast_engine/run.rs` around `hit_level_array.resize(..); hit_level_array.fill(..)` in the subject loop.
- NCBI relies on `diag_table->offset` updates and explicit clear logic (`blast_extend.c:84-105`, `blast_extend.c:159-176`) to avoid full clears where possible.
- Impact: large, repeated memory writes on long queries/subjects.

### 5) Output formatting allocates per hit
- LOSAT builds an intermediate `String` for each hit line, then `writeln!`.
  - LOSAT: `LOSAT/src/report/outfmt6.rs:229-289`.
- NCBI uses direct `snprintf` formatting into buffers (e.g., `align_format_util.cpp:986-1011`) and writes without intermediate heap allocations.
- Impact: significant overhead when output is large (megablast tends to be output-heavy).

### 6) Subject encoding and chunk handling are more eager
- LOSAT now streams subject records in single-thread runs and encodes per subject in `run.rs`; limit_lookup/parallel still preload subjects.
  - LOSAT: `LOSAT/src/algorithm/blastn/blast_engine/run.rs`, `LOSAT/src/algorithm/blastn/coordination.rs`.
- NCBI streams subject chunks via `s_GetNextSubjectChunk` and processes incrementally.
  - NCBI: `blast_engine.c:478-500`.
- Impact: reduced upfront memory bandwidth and allocation cost; remaining eager behavior is limited to limit_lookup/parallel paths.

### 7) Parallel path has aggregation overhead
- LOSAT's parallel path sends `Vec<Hit>` through a channel to a writer thread which concatenates all hits.
  - LOSAT: `LOSAT/src/algorithm/blastn/blast_engine/run.rs:6718-6790`.
- NCBI maintains `BlastHitList`/`BlastHSPList` structures per query/subject, avoiding channel fan-in.
  - NCBI: `blast_hits.h:153-180`.
- Impact: extra allocations, channel synchronization, and a large final merge step.

## Suggestions (NCBI-anchored, parity-safe)

### A) Replace per-hit ID storage with indices, resolve strings at output time
- Align with NCBI `BlastHSPList`/`BlastHitList` (no per-hit string ownership).
  - NCBI references: `blast_hits.h:153-180`.
- Expected benefit: smaller `Hit` objects, less `Arc` cloning, lower memory traffic.
- LOSAT targets:
  - `LOSAT/src/common.rs` (replace `query_id`/`subject_id` with indices only).
  - `LOSAT/src/report/outfmt6.rs` (look up IDs on write).

### B) Introduce NCBI-like `BlastHitList`/`BlastHSPList` containers for pruning
- Apply `max_hsps_per_subject` and `hitlist_size` in the same phase/order as NCBI.
  - NCBI references: `blast_hits.c:2049-2067`, `blast_traceback.c:850-892`.
- Expected benefit: reduce `FxHashMap` grouping, shrink the number of stored hits earlier.
- LOSAT targets:
  - `LOSAT/src/algorithm/blastn/blast_engine/run.rs` post-processing section.

### C) Reuse `offset_pairs` buffer across chunks
- Mirror `GetOffsetArraySize` usage in NCBI and keep the buffer in per-thread scratch.
  - NCBI: `lookup_wrap.c:255-288`.
- Expected benefit: reduce per-chunk allocation churn in the seed scan loop.
- LOSAT targets:
  - `LOSAT/src/algorithm/blastn/blast_engine/run.rs` (move buffer into `SubjectScratch`).

### D) Pre-select scan routine once, avoid per-call dispatch
- Match NCBI's `s_MBChooseScanSubject` strategy.
  - NCBI: `blast_nascan.c:2602-2677`.
- Expected benefit: reduce branching and closure overhead in hot scanning loops.
- LOSAT targets:
  - `LOSAT/src/algorithm/blastn/blast_engine/run.rs` scan dispatcher.

### E) Optimize diagonal table clearing
- Follow NCBI's offset-based clearing to avoid full `fill` on large arrays.
  - NCBI: `blast_extend.c:84-105`, `blast_extend.c:159-176`.
- Expected benefit: less memory bandwidth in two-hit tracking.
- LOSAT targets:
  - `LOSAT/src/algorithm/blastn/blast_engine/run.rs` (`hit_level_array`/`hit_len_array`).

### F) Stream subject encoding per chunk (or keep packed-only)
- Align with NCBI's incremental subject chunking (`s_GetNextSubjectChunk`).
  - NCBI: `blast_engine.c:478-500`.
- Expected benefit: lower upfront encoding cost and peak memory.
- Status: completed for single-thread via streamed subjects + per-subject encoding in `LOSAT/src/algorithm/blastn/blast_engine/run.rs`.

### G) Remove per-hit `String` assembly in outfmt6
- Write directly to the output buffer with numeric formatting and without intermediate strings.
  - NCBI: `align_format_util.cpp:986-1011`.
- Expected benefit: large speedups when output volume is high.
- LOSAT targets:
  - `LOSAT/src/report/outfmt6.rs` (`format_hit` + `write_outfmt6`).

## Practical Next Steps
1. Address A + B first (hit storage/pruning): largest expected impact on memory and sorting overhead.
2. Address C + D next (scan/offset buffer reuse): reduces hot-loop allocations/dispatch cost.
3. Address G (output formatting) if output volume is large in the benchmark.
4. Address F for large database runs or multi-subject FASTA inputs.

## Notes
- All suggestions above align with existing NCBI BLAST+ behavior and should preserve output parity if implemented carefully.
- If you want targeted timing data for LOSAT BLASTN, we can add scoped timers that mirror NCBI phase boundaries (e.g., scan, ungapped, gapped, traceback), but any instrumentation must be guarded to avoid output divergence.
