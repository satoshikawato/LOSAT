# BLASTN Traceback/Prune Performance Implementation Plan

Date: 2026-05-13

Goal: reduce LOSAT BLASTN cases that are slower than NCBI BLAST+ while
preserving bit-perfect parity. The first target is the `--task blastn`
traceback/prune-heavy path observed in self and high-HSP-count comparisons.

NCBI BLAST remains the only source of truth. NCBI binaries may be used only as
comparison oracles during diagnostics and validation. LOSAT runtime code must
not call NCBI BLAST, link to NCBI BLAST, or use NCBI BLAST as a fallback.

## Current Observation

The timing table shows mixed behavior:

- LOSAT is faster on some task-blastn cases, for example
  `MjPMNV_vs_MlPMNV`.
- LOSAT is slightly slower on several short jobs, where fixed costs can dominate.
- `NZ_CP006932_Self` is a larger regression: BLAST+ 5.215s vs LOSAT 9.986s in
  the supplied table.

Re-running `NZ_CP006932_Self` with `LOSAT_TIMING=1` showed:

| Stage | Time | Calls |
| --- | ---: | ---: |
| `scan_lookup` | 0.222s | 1 |
| `ungapped_extend` | 0.151s | 1,998,965 |
| `gapped_extend` | 0.856s | 66,004 |
| `traceback_prune` | 11.533s | 2 |
| `format_output` | 0.004s | 1 |

The exact wall time increases under timing instrumentation, but the ratio is
clear: the regression is not lookup scanning or ungapped extension. It is the
final HSP materialization and pruning region currently aggregated under
`traceback_prune`.

Primary Rust region:

- `LOSAT/src/algorithm/blastn/blast_engine/run.rs:7912-8703`

This region includes:

- traceback-containing gapped alignment
- ambiguity re-evaluation
- identity/length retesting
- common endpoint purge pass 1 and pass 2
- score re-sort
- final interval-tree containment pass

## Scope

Primary scope:

- `LOSAT blastn --task blastn`
- single-thread `-n 1` first, because the supplied comparison table is
  single-threaded
- outfmt 6/7 parity for hit order, coordinates, bit score, e-value, identity,
  mismatches, gap opens, and alignment length
- performance of self/high-overlap cases without regressing small pairwise cases

Secondary scope after task-blastn is stable:

- default megablast mode
- `-n > 1` behavior and writer-thread overhead
- multi-query/multi-subject databases

Out of scope:

- delegating any behavior to NCBI BLAST from LOSAT runtime code
- dropping parity-sensitive traceback or purge steps
- optimizing by changing cutoffs, x-drop, sorting keys, or final HSP filtering
- accepting hit-count thresholds as success. Hit counts are diagnostics only;
  final acceptance is bit-perfect parity.

## Primary NCBI References

Read these before changing Rust code. Every parity or performance code change
must include the relevant NCBI snippet comment with source path and line
numbers.

- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:350-613`
  - traceback loop, start offset selection, subject-range adjustment, final
    traceback, identity/length test
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:633-692`
  - post-traceback endpoint purge, ambiguity re-evaluation, score re-sort,
    final containment purge
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3810-4105`
  - preliminary gapped alignment, interval tree containment, seed selection,
    HSP initialization
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3323-3390`
  - `BlastGetStartForGappedAlignmentNucl`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:479-642`
  - `Blast_HSPReevaluateWithAmbiguitiesGapped`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2268-2535`
  - common-start/common-endpoint sort keys and purge behavior
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_itree.c:797-847`
  - `s_HSPIsContained`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_itree.c:930-995`
  - `BlastIntervalTreeContainsHSP`

Primary Rust files:

- `LOSAT/src/algorithm/blastn/blast_engine/run.rs`
- `LOSAT/src/algorithm/blastn/alignment/greedy.rs`
- `LOSAT/src/algorithm/blastn/alignment/gapped.rs`
- `LOSAT/src/algorithm/blastn/filtering/purge_endpoints.rs`
- `LOSAT/src/algorithm/blastn/interval_tree.rs`
- `LOSAT/src/algorithm/blastn/hsp.rs`

## Phase 0: Lock Measurement Conditions

Objective: make performance comparisons reproducible and ensure BLAST+ and
LOSAT are running the same task and options.

Tasks:

1. Use NCBI task-blastn oracle files for task-blastn comparisons:
   - NCBI: `LOSAT/tests/blast_out/NZ_CP006932.NZ_CP006932.task_blastn.out`
   - LOSAT: `LOSAT/tests/losat_out/NZ_CP006932.NZ_CP006932.losatn.blastn.out`
2. Keep NCBI default/megablast outputs separate from task-blastn outputs.
3. Record all options in a manifest row for each timing case:
   `task`, `reward`, `penalty`, `gap_open`, `gap_extend`, `word_size`,
   `dust`, `soft_masking`, `strand`, `evalue`, `max_target_seqs`,
   `max_hsps_per_subject`, `num_threads`, `outfmt`, query, and subject.
4. Run each timing case at least three times and report median wall time.
5. Capture both elapsed wall time and `LOSAT_TIMING=1` stage counters for LOSAT.

Acceptance:

- No task-blastn performance result is compared against a megablast oracle.
- A single command can reproduce each table row.
- Current parity status is known before any optimization is applied.

## Phase 1: Split `traceback_prune` Timing

Objective: identify the exact substage that dominates `traceback_prune`.

Current aggregate timing is too coarse. Add timing counters behind
`LOSAT_TIMING=1` only. The default runtime path must remain silent and
behaviorally unchanged.

Add counters for:

- `traceback_sort_prelim`
- `traceback_tree_precheck`
- `traceback_start_offsets`
- `traceback_alignment`
- `traceback_hsp_build`
- `purge_common_endpoint_pass1`
- `reevaluate_ambiguities`
- `identity_length_test`
- `purge_common_endpoint_pass2`
- `traceback_score_resort`
- `traceback_final_tree_contains`
- `traceback_final_tree_add`

Add counts for:

- preliminary HSPs entering traceback
- preliminary HSPs skipped by tree precheck
- HSPs receiving full traceback
- HSPs deleted by e-value/cutoff
- HSPs deleted by re-evaluation
- HSPs deleted by identity/length test
- HSPs removed by endpoint purge pass 1 and pass 2
- HSPs removed by final interval-tree containment
- max and median edit-script length
- max and median HSP alignment length

Acceptance:

- `NZ_CP006932_Self` explains at least 90% of `traceback_prune` time by named
  substages.
- The same counters are available for a fast case, such as
  `MjPMNV_vs_MlPMNV`, so regressions can be compared against a non-self case.
- `LOSAT_TIMING` does not change output bytes.

## Phase 2: Confirm First Parity Divergence Before Optimizing

Objective: avoid optimizing code that is still doing the wrong work.

Known current task-blastn artifact counts:

| Case | LOSAT hits | NCBI task-blastn hits |
| --- | ---: | ---: |
| `NZ_CP006932_Self` | 11,479 | 12,340 |
| `MjPMNV_vs_MlPMNV` | 54,310 | 54,402 |
| `MelaMJNV_vs_PemoMJNVA` | 2,708 | 2,729 |

The self case still has parity deltas. Before changing algorithms for speed,
use stage-level tracing to determine whether LOSAT is spending time on HSPs
that NCBI rejects earlier, or whether both implementations perform equivalent
work but LOSAT has higher constant factors.

Tasks:

1. Select one LOSAT-only HSP and one NCBI-only HSP for `NZ_CP006932_Self`.
2. Trace each through:
   - lookup hit retrieval
   - ungapped extension
   - preliminary gapped score-only extension
   - traceback tree precheck
   - final traceback
   - endpoint purge pass 1
   - ambiguity re-evaluation
   - endpoint purge pass 2
   - final interval-tree containment
3. Compare timing counters for HSP groups that are later deleted.

Acceptance:

- The first divergent stage is known for at least one LOSAT-only and one
  NCBI-only self-case HSP.
- If the first divergence is a parity bug, fix parity before performance work.
- If work is parity-equivalent, proceed to constant-factor optimization.

## Phase 3: Optimize Common Endpoint Purge

Objective: reduce purge cost without changing NCBI ordering or deletion
behavior.

Likely risks:

- repeated `Vec` movement during purge
- repeated sorting of large HSP vectors
- avoidable allocation while grouping by subject/query/endpoints

Tasks:

1. Measure `purge_hsps_with_common_endpoints_ex` input sizes, delete counts,
   trimmed counts, and time per pass.
2. Inspect NCBI `blast_hits.c:2268-2535` and confirm the Rust sort keys and
   deletion/trim behavior match exactly.
3. Replace repeated element removal with mark-and-compact or NCBI-equivalent
   pointer-array compaction only if parity allows.
4. Reuse scratch buffers for per-subject grouping and endpoint sorting.
5. Keep stable behavior for tie cases by matching NCBI comparators exactly.

Acceptance:

- Output is byte-identical to the pre-change LOSAT output on current fixtures.
- Output remains bit-perfect where LOSAT already matches NCBI.
- `purge_common_endpoint_pass1 + purge_common_endpoint_pass2` time drops on
  `NZ_CP006932_Self`.
- Fixed-cost overhead does not increase short cases by more than measurement
  noise.

## Phase 4: Optimize Interval-Tree Containment

Objective: reduce containment cost in dense overlap cases.

Tasks:

1. Measure time and call counts separately for:
   - preliminary tree precheck
   - final containment pass
   - tree insertion
2. Count midpoint-list traversal length and max tree depth in
   `BlastIntervalTree::contains_hsp`.
3. Compare Rust traversal against NCBI `blast_itree.c:930-995` and
   `s_HSPIsContained`.
4. Optimize only parity-preserving details:
   - reduce temporary object construction in hot loops
   - avoid repeated context-offset recomputation
   - reserve tree node capacity from observed HSP count
   - reuse interval tree buffers across subject chunks when NCBI timing/order
     permits

Acceptance:

- Containment decisions match current LOSAT behavior before any parity fix.
- After parity fixes, containment decisions match NCBI-traced decisions.
- Dense self-case final containment time drops without moving output order.

## Phase 5: Reduce Traceback Alignment Overhead

Objective: preserve NCBI traceback behavior while lowering constant factors in
full traceback alignment.

Tasks:

1. Split `traceback_alignment` timing by implementation:
   - greedy traceback path
   - DP traceback path
2. Measure edit-script allocation count, total edit operations, and max edit
   script length.
3. Inspect whether `GapEditOp` vectors are cloned or moved unnecessarily.
4. Reuse traceback scratch buffers where NCBI permits equivalent lifetime.
5. Avoid allocating final `gap_info` for HSPs that will immediately fail
   e-value/cutoff if the NCBI order allows that decision before ownership is
   finalized.

Acceptance:

- Traceback coordinates, score, edit script semantics, identity, mismatch, gap
  opens, and alignment length remain unchanged.
- `traceback_alignment` time drops in self and high-HSP cases.
- No shortcut skips NCBI-required ambiguity re-evaluation or identity testing.

## Phase 6: Optimize Ambiguity Re-Evaluation and Identity Counting

Objective: reduce cost of post-traceback re-evaluation while preserving
`Blast_HSPReevaluateWithAmbiguitiesGapped` and
`Blast_HSPTestIdentityAndLength` behavior.

Tasks:

1. Measure `reevaluate_hsp_with_ambiguities_gapped_ex` by total edit-script
   length, not just HSP count.
2. Confirm score-factor and gap-penalty handling against
   `blast_hits.c:508-608`.
3. Confirm identity counting against `blast_hits.c:745-818`.
4. Add a fast path only when the NCBI conditions prove it is equivalent, for
   example unambiguous contiguous `Sub(n)` spans with no trimming.
5. Keep fallback path identical for ambiguous bases, trims, and mixed edit
   operations.

Acceptance:

- Re-evaluation decisions and recomputed scores/e-values remain unchanged.
- Identity, mismatch, gap open, and length fields remain byte-identical.
- Self-case re-evaluation time drops, or the measurement proves it is not a
  meaningful bottleneck.

## Phase 7: Fixed-Cost Improvements for Short Cases

Objective: address small regressions where LOSAT is slower by 0.02-0.20s.

Only start this phase after traceback/prune counters are available, because
short-case regressions may be fixed by the same HSP-path changes.

Tasks:

1. Split startup/read/filter/build/search/write timing:
   - FASTA read
   - subject metadata scan
   - DUST
   - lookup build
   - output writer setup
2. Check whether subject FASTA is read or scanned more than necessary.
3. Reuse lookup and mask buffers within a run.
4. Avoid progress-bar or writer-thread setup when `num_threads == 1` and output
   is a plain file path.
5. Keep CLI behavior and verbose/debug output unchanged.

Acceptance:

- Short task-blastn cases no longer regress materially against BLAST+ when
  parity is maintained.
- Long cases do not slow down.

## Validation Matrix

Run after each phase that changes code:

| Case | Purpose |
| --- | --- |
| `NZ_CP006932.NZ_CP006932 --task blastn` | dense self-case stress |
| `MjPMNV.MlPMNV --task blastn` | high hit count where LOSAT is already faster |
| `MelaMJNV.PemoMJNVA --task blastn` | short fixed-cost regression |
| `PesePMNV.MjPMNV --task blastn` | very small output case |
| `SiNMV.ChdeNMV --task blastn` | moderate case with current LOSAT speedup |

For each case collect:

- wall time median over at least three runs
- user/sys time
- hit count
- byte-for-byte diff against the correct NCBI task oracle
- coordinate-key diff summary
- stage timing counters

## Implementation Order

1. Add substage timing and counters under `LOSAT_TIMING=1`.
2. Re-run the validation matrix and identify the dominant substage.
3. Trace first parity divergence for `NZ_CP006932_Self`.
4. Fix any confirmed parity bug before speed work.
5. Optimize common endpoint purge if it is a measured bottleneck.
6. Optimize interval-tree containment if it is a measured bottleneck.
7. Optimize traceback allocation/scratch if it is a measured bottleneck.
8. Optimize ambiguity re-evaluation and identity counting if measured.
9. Address fixed startup/read/build costs for short jobs.
10. Update this document with measured before/after tables.

## Stop Conditions

Stop and re-check NCBI source if any of these happen:

- a proposed optimization changes hit order, coordinates, scores, identity, or
  e-values
- a faster path skips a NCBI-required traceback, re-evaluation, or purge step
- a timing improvement depends on comparing LOSAT task-blastn to NCBI megablast
- a change improves `NZ_CP006932_Self` but regresses `MjPMNV_vs_MlPMNV`
  materially
- NCBI source for the corresponding behavior cannot be located

## Expected Outcome

The most likely first win is not lookup scanning. It is reducing unnecessary
work and constant factors inside the final HSP path:

- common endpoint purge
- interval-tree containment
- traceback edit-script allocation
- ambiguity re-evaluation

If Phase 1 shows `traceback_alignment` dominates, optimize scratch reuse and
edit-script ownership first. If endpoint purge or final containment dominates,
optimize those before touching alignment code. Parity remains the gate for all
performance changes.
