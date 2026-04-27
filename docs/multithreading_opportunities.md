# Multithreading Opportunities for LOSAT

## Scope

This note ranks multithreading opportunities for TBLASTX, BLASTN, and BLASTP
while preserving bit-identical single-thread output. NCBI BLAST+ source is the
reference for safe work partitioning and merge order.

Primary NCBI references:

- Preliminary search subject loop:
  `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1409-1498`
- Per-thread preliminary setup and local diagnostics:
  `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1633-1688`
- Multi-thread traceback batching:
  `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:1436-1725`

## Priority Ranking

| Priority | Program | Parallelization point | Expected speedup | Parity risk | LOSAT status |
|---:|---|---|---|---|---|
| 1 | BLASTP | Subject-level preliminary search: scan, ungapped extension, gapped preliminary extension, per-subject provisional HSPList construction | Very high. Current BLASTP subject loop is serial. | Medium | Not implemented |
| 2 | BLASTN | Subject-level scan, ungapped/gapped extension, traceback, pruning | High for multi-subject databases. | Low to medium | Implemented |
| 3 | TBLASTX | Subject-level six-frame translation, scan, ungapped extension, reevaluation, linking, filtering | High for multi-subject databases. | Low | Implemented |
| 4 | BLASTN | Chunk-level processing for very large subjects | High when one subject dominates runtime. | Medium | Implemented |
| 5 | BLASTP | Composition-based statistics / Kappa redo by query or HSPList batch | Medium to high, workload-dependent. | High | Not implemented |
| 6 | All | Record-level preprocessing: FASTA encoding, translation, SEG/DUST masking | Medium. Helps large inputs but less than search-stage parallelism. | Low | Partially not implemented |
| 7 | All | Alignment rendering / formatting precomputation, then deterministic final write | Small to medium, mostly output-heavy runs. | Medium | Not implemented |

## Program Notes

### BLASTP

Best immediate target: subject-level preliminary search.

LOSAT currently builds shared query/lookup/statistical state, then runs a serial
subject loop in `LOSAT/src/algorithm/blastp/blast_engine.rs:2505`. The heavy
work inside that loop includes amino-acid scanning, one-hit or two-hit ungapped
extension, initial HSP sorting, optional chaining, preliminary gapped extension,
and per-query provisional HSPList updates.

Parallelization shape:

- Keep query lookup tables, contexts, cutoffs, search spaces, and encoded
  subjects read-only.
- Give each worker its own `offset_pairs`, diagonal table, interval tree,
  initial HSP vector, preliminary HSP list, and composition/gap scratch where
  needed.
- Return per-subject, per-query `BlastpHspList` results from workers.
- Merge after all workers complete in deterministic NCBI order.

Do not parallelize inside a subject first. The subject-local diagonal table,
`diag_offset`, two-hit state, restricted alignment retry, and
`redo_index`/`redo_query` logic are order-sensitive.

Composition-based stats are a second phase opportunity, but higher risk. NCBI
has dedicated multi-thread redo alignment code in traceback, while LOSAT's
current flow uses one `BlastCompositionWorkspace` and processes queries and
matches in deterministic order. This should wait until subject-level BLASTP
parallelism is stable.

### BLASTN

Subject-level parallelism already exists in
`LOSAT/src/algorithm/blastn/blast_engine/run.rs:8190`. Each worker owns
`GapAlignScratch` and `SubjectScratch`; the writer thread receives subject hits,
updates `BlastnHitList`s, and performs final pruning/output.

Chunk-level parallelism for split long subjects also exists around
`LOSAT/src/algorithm/blastn/blast_engine/run.rs:7436`. This is useful when a
single long subject would otherwise leave most threads idle.

Potential improvement:

- The current parallel writer updates hitlists as messages arrive. For strict
  thread-count-independent behavior, prefer collecting per-subject results with
  their `oid`, sorting by `oid`, and then applying `Blast_HitListUpdate` in that
  deterministic order.

### TBLASTX

Subject-level parallelism already exists in
`LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs:1522`. Each worker owns
offset-pair and diagonal scratch. Hits are collected and final output is sorted
with NCBI-compatible query, subject, and HSP comparators.

Do not split a subject frame loop across threads unless NCBI parity is proven.
The six subject frames share diagonal state and `diag_offset` behavior within a
subject, so frame-level parallelism can alter which seeds trigger extensions.

## Bit-Identical Output Requirements

These rules are mandatory for any parallel implementation:

- Workers must never write final output directly.
- Merge must not depend on channel receive order, Rayon scheduling, or thread
  count.
- Final output order must be reconstructed with NCBI-compatible comparators:
  query order, subject HSPList order, then HSP score/coordinate order.
- Subject-local scan order must remain unchanged inside each subject.
- Diagonal-table state must be worker-local and reset at the same timing as
  NCBI.
- Floating-point calculations must run in the same per-HSP/per-HSPList order
  used by the single-thread implementation where order can affect admission,
  heap replacement, or E-value pruning.
- Hitlist pruning must be deterministic. If `Blast_HitListUpdate` behavior is
  order-sensitive on ties, feed HSPLists in NCBI subject order after parallel
  collection.

## Implementation Plan

### 1. BLASTP subject-level parallelism

Refactor the body of the BLASTP subject loop into a function that returns
per-subject provisional HSPLists instead of mutating global hitlists directly.
Then add a Rayon path analogous to BLASTN/TBLASTX.

Safe merge strategy:

1. Worker returns `(s_idx, Vec<(q_idx, BlastpHspList)>)`.
2. Main thread sorts worker results by `s_idx`.
3. Main thread applies existing `BlastpHitList::update` in that sorted order.
4. Existing composition/Kappa post-processing and final output order remain
   unchanged.

This is the largest unimplemented speedup and should be done before smaller
formatting or preprocessing work.

### 2. Harden BLASTN deterministic merge

BLASTN is already parallel, but the channel writer currently updates hitlists
as messages arrive. Change the parallel path to collect `(s_idx, hits)`, sort by
`s_idx`, then call `update_hitlists_with_subject_hits`. This should remove any
residual scheduling dependency while keeping the same subject-level work split.

### 3. Add preprocessing parallelism

Parallelize independent encoding and masking passes after search-stage
parallelism is stable:

- BLASTP subject protein encoding.
- TBLASTX subject/query translation and SEG masking where data ownership is
  clean.
- BLASTN subject packed/blastna cache construction on paths that still preload
  subjects.

### 4. Output rendering only after search parity is locked

Formatting can be parallelized by building rendered lines or alignment views
per hit, then writing in final NCBI order. This should be last because it adds
buffering and can mask ordering bugs in search-stage changes.

## Open Questions

- For BLASTP, should the first implementation target tabular output only, or
  include pairwise `outfmt 0` in the same pass?
- Should `-num_threads 0` continue to mean all logical CPUs for LOSAT, or should
  it be capped to match NCBI's CLI behavior on the target platform?
- For BLASTN, do we want to prioritize deterministic merge hardening before
  adding any new parallel work?

