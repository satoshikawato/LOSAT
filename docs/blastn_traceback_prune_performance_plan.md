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

## Implementation Update: 2026-05-13

This update focused on actual BLASTN traceback/prune hot-path changes rather
than adding diagnostics or harness-only code. NCBI source comments were added
next to the Rust implementation changes, following the repository porting
requirements.

Changed hot paths:

- `LOSAT/src/algorithm/blastn/filtering/purge_endpoints.rs`
  - Added a per-subject fast path for
    `purge_hsps_with_common_endpoints_ex`, matching the normal
    `BlastHSPList` shape used by NCBI traceback.
  - Replaced repeated `Vec::remove` movement during common-start/common-end
    purge with grouped compaction that preserves the NCBI active-prefix and
    tail-order behavior from `blast_hits.c:2268-2535`.
- `LOSAT/src/algorithm/blastn/interval_tree.rs`
  - Matched NCBI's initial interval-tree node capacity of 100.
  - Added `reserve_nodes_for_hsps` and used it before dense traceback/final
    containment passes to avoid repeated `Vec` reallocations.
  - Made `TreeHsp` `Copy` to remove hot-path clones.
- `LOSAT/src/algorithm/blastn/blast_engine/run.rs`
  - Reserved interval-tree storage from the observed preliminary/final HSP
    counts before containment passes.
- `LOSAT/src/algorithm/blastn/alignment/gapped.rs`
  - Optimized traceback DP cell access by using raw pointer access for
    `score_array` and unchecked score-matrix access after the NCBI sentinel
    checks.
  - Optimized BLASTNA identity counting by processing contiguous substitution
    spans instead of checking bounds per aligned base.

Rejected implementation:

- A `MaybeUninit`-based trace-row allocation path was tested to avoid row
  zero-fill. It either regressed the dense self case or made capacity handling
  fragile, so it was removed from the final patch.

Validation performed:

- `cargo check`
- `cargo build --release`
- `cargo test test_purge_common_start_keeps_longer_end_on_score_tie`
- `cargo test test_interval_tree_basic`
- `PesePMNV.MjPMNV --task blastn -n 1` remained byte-identical to the current
  LOSAT fixture.
- `NZ_CP006932.NZ_CP006932 --task blastn -n 1` remained byte-identical to the
  pre-final-DP-hot-loop LOSAT output, so the final DP pointer optimization did
  not introduce a new output delta.

Observed performance on the instrumented `NZ_CP006932.NZ_CP006932 --task blastn
-n 1` run:

| Measurement point | `traceback_alignment` | `traceback_prune` | Total |
| --- | ---: | ---: | ---: |
| Before final DP pointer optimization | 11.422s | 11.564s | 13.079s |
| After final DP pointer optimization | 8.403s | 8.546s | 9.995s |

The largest retained win came from reducing constant factors inside full
traceback DP, not from changing BLASTN filtering order or skipping any NCBI
required purge/re-evaluation step.

Current caveat:

- The current `NZ_CP006932.NZ_CP006932 --task blastn` LOSAT output still differs
  from the existing fixture/oracle in this working tree. The final DP hot-loop
  change did not create that delta, but parity for this dense self case remains
  a blocker before treating the performance work as complete.

## Next Direction

Immediate priorities:

1. Resolve the remaining `NZ_CP006932.NZ_CP006932 --task blastn` parity delta
   before accepting further risky speedups. Focus on the first divergent stage
   already identified in this plan: preliminary HSP handling, traceback tree
   containment, common endpoint trimming, ambiguity re-evaluation, and final
   containment.
2. Re-run the validation matrix with at least three wall-time runs per case
   after the parity delta is understood. Keep task-blastn outputs separate from
   megablast outputs.
3. If parity is confirmed for the changed regions, continue with hot-path
   improvements in this order:
   - Reuse or pre-size traceback edit-script storage only where the NCBI
     lifetime/order permits identical output.
   - Reduce temporary `TreeHsp` construction and repeated coordinate/context
     recomputation in interval-tree containment.
   - Add narrowly proven fast paths in ambiguity re-evaluation and identity
     counting for unambiguous contiguous substitution spans.
   - Address short-case fixed costs only after the dense traceback path is no
     longer the dominant regression.

Do not pursue optimizations that change hit order, score/e-value rounding,
coordinate adjustment, endpoint trimming order, or the timing of NCBI-required
re-evaluation. Any further unsafe optimization should first be tied to the
corresponding NCBI source block and then validated against byte-for-byte LOSAT
output before comparing against NCBI BLAST+.

## Implementation Update: 2026-05-13, In-Place Endpoint Trim

This update kept the BLASTN traceback/prune behavior on the NCBI path and
focused on reducing hot-path allocation in endpoint trimming.

Changed hot path:

- `LOSAT/src/algorithm/blastn/filtering/purge_endpoints.rs`
  - Updated `cut_off_gap_edit_script` to move the HSP `gap_info` vector out of
    the HSP, mutate it, and put it back instead of cloning the whole edit
    script before every trim.
  - Replaced the temporary `new_ops` allocation with in-place `drain` for
    `cut_begin=true` and `truncate` for `cut_begin=false`.
  - This follows NCBI `s_CutOffGapEditScript`
    (`blast_hits.c:2392-2452`), where `GapEditScript *esp = hsp->gap_info` is
    rewritten in place and `esp->size` is adjusted after moving or truncating
    the active edit-script prefix.

Rejected implementation:

- An explicit fast path for a single exact `Sub(n)` A/C/G/T HSP during
  `Blast_HSPReevaluateWithAmbiguitiesGapped` was removed. NCBI BLAST does not
  contain that branch; it always runs the `blast_hits.c:541-612` per-base
  substitution loop for `eGapAlignSub` and the `blast_hits.c:619-637`
  extension loop. Even though the shortcut can be reasoned about as equivalent
  for a narrow case, it is not an NCBI control-flow equivalent and should not be
  used until there is an explicit source-backed justification accepted for this
  repository.

Validation performed:

- `rustfmt LOSAT/src/algorithm/blastn/filtering/purge_endpoints.rs --edition 2021`
- `cargo check`
- `cargo test test_purge_common_start_keeps_longer_end_on_score_tie`
- `cargo build --release`
- `PesePMNV.MjPMNV --task blastn -n 1` remained byte-identical to the current
  LOSAT fixture.
- `MelaMJNV.PemoMJNVA --task blastn -n 1` remained byte-identical to the
  current LOSAT fixture in both debug `cargo run` and `target/release/LOSAT`.

Current caveat:

- This change removes edit-script clone/allocation work only when common
  endpoint trimming actually calls `s_CutOffGapEditScript`-equivalent logic.
  It is not expected to move cases dominated by full traceback DP unless those
  cases also have substantial trimmed-HSP traffic.

Next direction:

1. Keep re-evaluation and identity-counting control flow literal to NCBI unless
   a proposed shortcut exists in the NCBI source or is accepted as a strictly
   implementation-level Rust optimization with byte-for-byte proof.
2. Measure `purge_common_endpoint_pass1` and `purge_common_endpoint_pass2` on
   `NZ_CP006932.NZ_CP006932 --task blastn -n 1` before and after this endpoint
   trim change to quantify whether clone removal matters in the dense self
   case.
3. Continue with NCBI-shaped allocation reductions:
   - Reuse existing edit-script buffers where NCBI mutates `GapEditScript` in
     place.
   - Avoid temporary `TreeHsp` or coordinate recomputation only where the
     `BlastIntervalTreeContainsHSP/AddHSP` order and inputs remain identical.
   - Keep final containment, endpoint purge, ambiguity re-evaluation, and
     identity/length testing in the same order as `blast_traceback.c:633-692`.

## Implementation Update: 2026-05-13, Direct Traceback Edit-Script Build

This update continued the hot-path work inside BLASTN traceback rather than
adding diagnostics or harness-only code. The changed code remains on the NCBI
traceback control path and keeps the NCBI-required order of traceback,
endpoint purge, ambiguity re-evaluation, identity/length testing, score resort,
and final interval-tree containment.

Changed hot path:

- `LOSAT/src/algorithm/blastn/alignment/gapped.rs`
  - Replaced the traceback edit-script construction path that first accumulated
    `(op_type, count)` tuples and then converted them into `GapEditOp` values.
    Traceback now appends directly into `Vec<GapEditOp>` with the same
    run-length merge semantics as `GapPrelimEditBlockAdd`
    (`gapinfo.c:174-185`).
  - Kept the NCBI forward/reverse traceback ordering from
    `blast_gapalign.c:689-727` and `blast_gapalign.c:2494-2496`: forward
    traceback is reversed after construction, while reverse traceback is already
    forward.
  - Reused the left traceback vector as the final combined edit-script buffer
    and moved right traceback operations into it, preserving the NCBI merge
    behavior from `blast_gapalign.c:2498-2534` while avoiding an extra left-side
    copy and right-side slice copy.
  - Replaced the post-combine terminal-gap trimming copy
    (`combined_edit_ops[start_idx..end_idx].to_vec()`) with in-place
    `truncate`/`drain`, matching NCBI's in-place shortening in
    `blast_gapalign.c:4683-4712`.

Validation performed:

- `cargo check`
- `cargo build --release`
- `PesePMNV.MjPMNV --task blastn -n 1` remained byte-identical to the current
  LOSAT fixture.
- `MelaMJNV.PemoMJNVA --task blastn -n 1` remained byte-identical to the
  current LOSAT fixture.
- `MjPMNV.MlPMNV --task blastn -n 1` remained byte-identical to the current
  LOSAT fixture.

Current caveat:

- No new reliable timing table was recorded for this update. The implementation
  removes allocation and copy work on the traceback edit-script path, but the
  measured wall-clock delta should be collected in a quiet environment with the
  validation matrix before treating this as a quantified performance win.

Next direction:

1. Re-run `NZ_CP006932.NZ_CP006932 --task blastn -n 1` with `LOSAT_TIMING=1`
   in a quiet environment and compare `traceback_alignment`,
   `traceback_hsp_build`, and `reevaluate_ambiguities` against the previous
   baseline. This edit-script change should primarily affect traceback
   allocation/copy overhead, not purge ordering.
2. If byte identity remains stable on the dense self case, inspect the remaining
   traceback allocation sites that still clone or copy `GapEditOp` buffers:
   reevaluation partial-copy (`blast_hits.c:460-470` equivalent) and endpoint
   trimming tails.
3. Do not add shortcuts to ambiguity re-evaluation unless they are either found
   in NCBI source or proven as a Rust-only storage optimization that preserves
   the exact `Blast_HSPReevaluateWithAmbiguitiesGapped` loop decisions.
4. Continue prioritizing allocation/copy reduction in existing NCBI-shaped hot
   paths before touching scoring, sort keys, endpoint selection, or final
   containment behavior.
