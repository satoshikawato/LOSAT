# LOSATN Megablast Acceleration Plan (EDL933 vs Sakai, Single-Thread)

## Purpose
Define a parity-safe performance improvement plan for LOSATN megablast on the
EDL933 vs Sakai dataset, preserving bitwise output parity with NCBI BLAST+.

This plan is constrained by LOSAT's parity rules: NCBI BLAST C/C++ is the only
source of truth, and all algorithmic behavior must match NCBI exactly. Any
behavioral change requires direct NCBI code verification.

## Scope
- Mode: single-thread (`-num_threads 1`), megablast task.
- Dataset: `EDL933` (query) vs `Sakai` (subject).
  - Query: `LOSAT/tests/fasta/EDL933.fna`
  - Subject: `LOSAT/tests/fasta/Sakai.fna`
  - Baseline NCBI output samples: `LOSAT/tests/blast_out/EDL933.Sakai.blastn.megablast.*`
- Target: reduce runtime gap vs NCBI without altering any output fields.

## Baseline (from existing investigation)
Source: `docs/blastn_performance_investigation.md`
- Megablast EDL933 vs Sakai:
  - NCBI BLAST+: ~0.998s
  - LOSAT: ~3.359s
  - Gap: ~3.37x slower

## Parity Requirements (Non-Negotiable)
- Output must match NCBI BLAST+ bit-for-bit.
- Timing/order of algorithm phases must mirror NCBI.
- Floating-point precision must match NCBI.
- No new features or behavior deviations.
- All code changes must include NCBI reference comments with file path and line
  numbers (captured during implementation).

## Primary NCBI References (Verify Line Numbers During Implementation)
These are the authoritative references for behavior and hot-path design:
- `c++/src/algo/blast/core/blast_engine.c`
- `c++/src/algo/blast/core/blast_nascan.c`
- `c++/src/algo/blast/core/lookup_wrap.c`
- `c++/src/algo/blast/core/blast_hits.c`
- `c++/src/algo/blast/core/blast_hits.h`
- `c++/src/algo/blast/core/blast_extend.c`
- `c++/src/algo/blast/core/blast_traceback.c`
- `c++/src/algo/blast/core/blast_gapalign.c`
- `c++/src/algo/blast/core/blast_parameters.c`
- `c++/src/algo/blast/core/blast_stat.c`

Local quick-reference copies:
- `.ncbi_ref/blast_engine.c`
- `.ncbi_ref/blast_gapalign.c`
- `.ncbi_ref/blast_parameters.c`
- `.ncbi_ref/blast_setup.c`
- `.ncbi_ref/na_ungapped.c`

## Known LOSAT Hotspots and NCBI-Alignment Opportunities
Source: `docs/blastn_performance_investigation.md`

1) Hit storage and pruning
- Issue: LOSAT stores `Arc<str>` IDs per hit and defers pruning to a global pass.
- NCBI: `BlastHSPList` stores indexes, prunes per query/subject during traceback.
- LOSAT files:
  - `LOSAT/src/common.rs`
  - `LOSAT/src/algorithm/blastn/blast_engine/run.rs`
- NCBI references:
  - `blast_hits.h` (Hit/HSP list structure)
  - `blast_hits.c` (trim by `max_hsps_per_subject`)
  - `blast_traceback.c` (hitlist pruning order)

2) Offset-pair buffer reallocation
- Issue: LOSAT allocates `Vec<OffsetPair>` per subject chunk.
- NCBI: one-time sizing via `GetOffsetArraySize`, reused buffer.
- LOSAT file: `LOSAT/src/algorithm/blastn/blast_engine/run.rs`
- NCBI reference: `lookup_wrap.c`

3) Scan dispatch overhead
- Issue: LOSAT uses per-call match/closure dispatch in scan loop.
- NCBI: chooses scan routine once with a function pointer.
- LOSAT file: `LOSAT/src/algorithm/blastn/blast_engine/run.rs`
- NCBI reference: `blast_nascan.c`

4) Diagonal table clearing
- Issue: LOSAT clears full arrays each chunk.
- NCBI: avoids full clears using offset-based approach.
- LOSAT file: `LOSAT/src/algorithm/blastn/blast_engine/run.rs`
- NCBI reference: `blast_extend.c`

5) Output formatting
- Issue: LOSAT allocates per-hit `String` and then writes it.
- NCBI: direct buffer formatting without intermediate allocations.
- LOSAT file: `LOSAT/src/report/outfmt6.rs`
- NCBI reference: `align_format_util.cpp` (formatting helpers)

6) Subject encoding strategy
- Issue: LOSAT encodes all subjects up-front.
- NCBI: streams subject chunks on demand.
- LOSAT file: `LOSAT/src/algorithm/blastn/coordination.rs`
- NCBI reference: `blast_engine.c`

## Implementation Plan (Parity-First, Single-Thread)
The order below aims to maximize impact while preserving NCBI parity.

## Status (Current)
- Phase 0: Completed (timing instrumentation in `LOSAT/src/algorithm/blastn/blast_engine/run.rs`).
- Phase 1: Completed (NCBI-style hitlist update/heap, prelim hitlist sizing, per-subject hitlists, post-traceback filter/sort/prune order; unit tests added in `LOSAT/src/algorithm/blastn/hsp.rs`).
- Phase 2: Completed (offset_pairs preallocated per thread using `GetOffsetArraySize`, reused across chunks).
- Phase 3: Completed (preselected megablast scan routine per subject chunk to avoid per-call dispatch in scan loop).
- Phase 4: Pending.
- Phase 5: Completed (direct buffered outfmt6 writes in `LOSAT/src/common.rs`, NCBI-style scientific formatting + tests in `LOSAT/src/report/outfmt6.rs`).
- Phase 6: Pending.

### Phase 0: Baseline and Instrumentation (No Behavior Change)
Goal: measure and attribute runtime without affecting output.
- Add timing hooks guarded by `LOSAT_TIMING=1` and only when output is unchanged.
- Mirror NCBI phase boundaries:
  - lookup/scan
  - ungapped extension
  - gapped extension
  - traceback and pruning
  - formatting/output
- Ensure instrumentation is side-effect-free and compiled out or disabled by
  default to avoid output divergence.

Deliverables:
- Timing log entry points in `LOSAT/src/algorithm/blastn/blast_engine/run.rs`.
- Confirm output is identical when `LOSAT_TIMING` is disabled.

### Phase 1: Hit Representation and Pruning (Largest Expected Gain)
Goal: reduce per-hit allocation and prune early, matching NCBI.
- Replace per-hit `Arc<str>` ID storage with index-based references.
  - Use query/subject indices, resolve to strings only at output time.
- Introduce NCBI-like `BlastHitList`/`BlastHSPList` semantics for pruning:
  - Apply `max_hsps_per_subject` and `hitlist_size` in the same order as NCBI.
  - Align with NCBI's pruning timing in traceback.

NCBI parity checkpoints:
- `blast_hits.h` and `blast_hits.c` for data structures and trim logic.
- `blast_traceback.c` for pruning order and comparator behavior.

### New Findings (Phase 1)
- `GetPrelimHitlistSize` expands preliminary hitlist size for gapped searches and honors `ADAPTIVE_CBS`; this impacts prune order and must be mirrored exactly.
- `Blast_HitListUpdate` heapifies by e-value after sorting each `BlastHSPList` by e-value; using score order here changes pruning decisions.
- `Blast_HSPListSubjectBestHit` assumes the list is already score-sorted and does not sort internally; calling it after any re-ordering changes results.
- `Blast_TrimHSPListByMaxHsps` only truncates and does not recompute `best_evalue`; NCBI keeps the prior value through post-filter sorting.

### Phase 2: Offset-Pair Buffer Reuse
Goal: reduce allocator overhead in the hot scan loop.
- Move `offset_pairs` into per-thread scratch storage and reuse across chunks.
- Size buffer using NCBI-equivalent logic (`GetOffsetArraySize`).

NCBI parity checkpoints:
- `lookup_wrap.c` for size calculation and buffer lifetime expectations.

### New Findings (Phase 2)
- `GetOffsetArraySize` sizing can be computed once from the chosen lookup table and reused for all subjects/chunks; this mirrors NCBI's aux-struct allocation model and removes per-chunk reserve calls.

### Phase 3: Scan Dispatch Simplification
Goal: avoid per-call dynamic dispatch in scan loop.
- Preselect scan routine once (matching NCBI `s_MBChooseScanSubject` logic).
- Store function pointer or enum outside the hot loop.

NCBI parity checkpoints:
- `blast_nascan.c` scan selection logic and edge cases.

### New Findings (Phase 3)
- Scan routine selection can be cached once per subject chunk; masked subjects fall back to the generic scan path (mirrors `BlastChooseNucleotideScanSubjectAny` behavior).

### Phase 4: Diagonal Table Clearing Optimization
Goal: avoid full clears on large tables.
- Implement NCBI-style offset-based clearing (no changes to hit detection).
- Verify all counters and sentinel behavior are identical.

NCBI parity checkpoints:
- `blast_extend.c` diagonal table layout and clearing method.

### Phase 5: Output Formatting (If Output-Heavy)
Goal: reduce per-hit formatting cost.
- Replace per-hit `String` assembly with direct buffered writes.
- Match NCBI formatting (precision, rounding, order).

NCBI parity checkpoints:
- `align_format_util.cpp` and `outfmt` routines for numeric formatting.

### New Findings (Phase 5)
- NCBI scientific notation uses sign + zero-padded exponents via `snprintf("%.*e")` and `NStr::DoubleToString(..., fDoubleScientific)`; Rust's default `e` formatter omits padding and `+`.
- Fixture output (`EDL933.Sakai.blastn.megablast.out`) includes `e+05` and `e-05`, so tests must assert the padded exponent format.

### Phase 6: Subject Encoding Streaming (Optional, Larger Refactor)
Goal: reduce upfront encoding overhead for large subjects.
- Stream subject chunks (`s_GetNextSubjectChunk` pattern).
- Keep packing and padding identical to NCBI.

NCBI parity checkpoints:
- `blast_engine.c` subject chunk lifecycle, boundary conditions.

## Acceptance Criteria
Performance:
- Reduce LOSAT runtime for EDL933 vs Sakai megablast by at least 2x (target:
  within 1.3x of NCBI on the same machine), without regressions on other cases.

Parity:
- Bitwise-identical output vs NCBI BLAST+ for the same inputs and parameters.
- No changes to hit counts, scores, E-values, or coordinates.

Stability:
- No new memory leaks or large memory spikes on large subjects.

## Test Plan (Parity-Driven)
Testing is only allowed after known NCBI divergences are resolved per AGENTS.md.

### Unit Tests (Targeted)
- Add unit tests for any new NCBI-ported functions or data structures.
- Reference NCBI unit tests when available:
  `ncbi-blast/c++/src/algo/blast/unit_tests/`.

### Integration Tests (Parity)
- Compare LOSAT output with NCBI for EDL933 vs Sakai megablast.
  - Use `LOSAT/tests/blast_out/EDL933.Sakai.blastn.megablast.*` as references.
- Verify full output parity: hit counts, bit scores, E-values, coordinates,
  and ordering.

### Performance Benchmarking
- Re-run EDL933 vs Sakai megablast in single-thread mode.
- Record timings with `LOSAT_TIMING=1` for phase attribution.
- Track regressions or improvements per phase.

## Risks and Mitigations
- Risk: pruning logic changes hit ordering or inclusion.
  - Mitigation: strict NCBI reference mapping and unit tests for tie-breakers.
- Risk: buffer reuse introduces stale data bugs.
  - Mitigation: explicit initialization matching NCBI's clearing semantics.
- Risk: formatting changes introduce floating-point rounding drift.
  - Mitigation: match NCBI formatting routines and precision constants.

## Proposed Work Breakdown (Concrete Steps)
1) Baseline timing instrumentation (guarded).
2) Hit representation refactor + NCBI-aligned pruning.
3) Offset-pair buffer reuse.
4) Scan dispatch preselection.
5) Diagonal table clearing optimization.
6) Output formatting optimization (if output-heavy).
7) Optional subject streaming refactor.

Each step must include NCBI reference comments with file path and line numbers
immediately above any new or modified Rust code.
