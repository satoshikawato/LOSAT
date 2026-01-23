# LOSATN Megablast Additional Opportunities (Parity-Safe)

## Purpose
Document additional performance opportunities found in a code scan beyond the
existing acceleration plan, and outline a parity-safe implementation order.

## Scope
- Task: blastn (megablast), single-thread focus.
- Requirement: bit-perfect parity with NCBI BLAST+.
- Source of truth: NCBI BLAST C/C++ code only.

## Findings (Additional Opportunities)

Status legend: [x] done, [ ] pending.

### 1) [x] Index-only hit identifiers (remove per-hit Arc clones)
- Status: done. `Hit` now carries `q_idx`/`s_idx` only; IDs are resolved at output time.
- Observation (pre-fix): `Hit` stored `Arc<str>` for `query_id`/`subject_id`,
  cloned into each hit and re-cloned in `BlastnHsp::into_hit`.
- LOSAT refs: `LOSAT/src/common.rs:80`, `LOSAT/src/common.rs:138`,
  `LOSAT/src/algorithm/blastn/hsp.rs:146`, `LOSAT/src/report/outfmt6.rs:395`.
- NCBI refs: `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/include/algo/blast/core/blast_hits.h:154`,
  `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/include/algo/blast/core/blast_hits.h:155`
  (`BlastHSPList` stores `oid` and `query_index`, not strings).
- Implementation: carry only indices in `Hit`/`BlastnHsp`; resolve IDs at output time.

### 2) [ ] Replace string-keyed caches with index-keyed structures in re-evaluation
- Observation: `filter_hsps` uses `FxHashMap<(Arc<str>, Arc<str>)>` for
  sequences and encoded caches, plus an `id_map` keyed by indices.
- LOSAT refs: `LOSAT/src/algorithm/blastn/blast_engine/mod.rs:43`,
  `LOSAT/src/algorithm/blastn/blast_engine/mod.rs:163`.
- NCBI refs: `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:651`,
  `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:657`
  (uses direct pointers from `query_blk`/`subject` for re-evaluation).
- Plan idea: pass pre-encoded sequences into `filter_hsps`, and use
  `(q_idx, s_idx)` indexing (Vec or small array) instead of hashing strings.

### 3) [ ] Avoid per-chunk allocation in prelim-hit processing
- Observation: `prelim_hits.drain(..).collect()` allocates a new Vec per chunk.
- LOSAT refs: `LOSAT/src/algorithm/blastn/blast_engine/run.rs:6724`.
- NCBI refs: `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:491`
  (`BlastInitHitListReset` resets in place).
- Plan idea: reuse a scratch Vec by swapping/taking, preserving capacity.

### 4) [ ] Reduce subject chunk range allocations
- Observation: each `SubjectChunk` owns `seq_ranges: Vec<(i32, i32)>`; ranges are
  rebuilt for every chunk even when unmasked.
- LOSAT refs: `LOSAT/src/algorithm/blastn/blast_engine/run.rs:157`,
  `LOSAT/src/algorithm/blastn/blast_engine/run.rs:4244`.
- NCBI refs: `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:150`,
  `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:154`
  (`s_BackupSubject` keeps existing ranges).
- Plan idea: reuse shared ranges for unmasked subjects and avoid per-chunk Vecs
  (e.g., borrow slices or use a small fixed-size buffer for the common case).

### 5) [ ] Buffered output for outfmt 0/7 (if output-heavy)
- Observation: outfmt6 is optimized, but pairwise/outfmt7 still use many
  `format!`/`writeln!` calls and temporary strings.
- LOSAT refs: `LOSAT/src/report/pairwise.rs`, `LOSAT/src/report/outfmt6.rs:597`.
- NCBI refs: `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/align_format/showalign.cpp:1414`,
  `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/objtools/align_format/showalign.cpp:1435`
  (buffered formatting to stream).
- Plan idea: add a buffered write path similar to outfmt6 for outfmt0/7.

### 6) [ ] Reuse packed subject encoding when limit_lookup preloads subjects
- Observation: `build_db_word_counts` encodes subjects; the main scan encodes
  again when `load_subjects` is true.
- LOSAT refs: `LOSAT/src/algorithm/blastn/lookup.rs:156`,
  `LOSAT/src/algorithm/blastn/blast_engine/run.rs:4238`.
- NCBI refs: `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:1191`,
  `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:1208`,
  `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_nalookup.c:1214`
  (DB scan is separate, but NCBI reads from seqsrc rather than preloaded records).
- Plan idea: when subjects are already in memory, cache packed `ncbi2na` once and
  share across the word-count scan and main scan.

## Proposed Plan (Suggested Order)
1) [x] Index-only hit identifiers (largest memory/alloc win, matches NCBI layout).
2) [ ] Index-keyed re-evaluation caches and pass pre-encoded sequences.
3) [ ] Reuse prelim-hit scratch Vecs to avoid per-chunk allocations.
4) [ ] Reduce subject chunk range allocations for unmasked subjects.
5) [ ] Buffered output for outfmt0/7 (only if output-heavy workloads).
6) [ ] Optional packed subject reuse for limit_lookup/preload mode.

## Notes / Guardrails
- Every code change must include NCBI reference comments with file path and
  line numbers, per AGENTS.md.
- Do not run tests until known parity divergences are resolved or requested.
