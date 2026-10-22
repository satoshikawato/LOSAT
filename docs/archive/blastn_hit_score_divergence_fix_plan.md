# BLASTN Hit Count and Score Distribution Divergence Fix Plan

Date: 2026-05-14

## Goal

Make LOSATN `blastn --task blastn` match NCBI BLASTN bit-for-bit for final
outfmt 6/7 output: hit count, hit order, coordinates, raw score, bit score,
E-value, identity, mismatch count, gap opens, and alignment length.

NCBI BLAST is the source of truth. NCBI binaries may be used only as diagnostic
or comparison oracles. LOSAT runtime code must remain a pure Rust
implementation and must not call, link, wrap, or fall back to NCBI BLAST.

## Fixed Comparison Conditions

Do not compare LOSATN `--task blastn` against NCBI default `blastn`, because
NCBI default `blastn` is megablast. Task-blastn comparisons must use NCBI
`-task blastn`.

Primary manifest:

- `LOSAT/tests/blastn_parity_manifest.tsv`

Primary comparison helper:

- `LOSAT/tests/compare_blastn_parity.py`

Relevant NCBI references:

- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:348-374`
  - task cutoff update and blastn query-length doubling for cutoff search space
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_setup.c:821-847`
  - length adjustment and per-context effective search space assignment
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1887-1890`
  - final HSP E-value uses `query_info->contexts[hsp->context].eff_searchsp`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:633-692`
  - post-traceback common-endpoint purge, ambiguity reevaluation, second purge,
    score resort, and final interval-tree containment
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3323-3390`
  - `BlastGetStartForGappedAlignmentNucl`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:735-962`
  - `Blast_SemiGappedAlign`

## Current Baseline

Command used from crate root:

```bash
python tests/compare_blastn_parity.py --manifest tests/blastn_parity_manifest.tsv --limit 5
```

Current task-blastn artifacts:

| Case | NCBI hits | LOSATN hits | Common coordinate keys | NCBI-only | LOSAT-only | Same-coordinate bit score diffs |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `MelaMJNV.PemoMJNVA.task_blastn` | 2729 | 2708 | 2675 | 54 | 33 | 0 |
| `PemoMJNVA.PeseMJNV.task_blastn` | 2940 | 2931 | 2859 | 81 | 72 | 0 |
| `NZ_CP006932.NZ_CP006932.task_blastn` | 12340 | 11479 | 11144 | 1196 | 335 | 0 |

Current megablast artifacts are not yet at the same parity level:

| Case | NCBI hits | LOSATN hits | Common coordinate keys | Same-coordinate bit score diffs |
| --- | ---: | ---: | ---: | ---: |
| `EDL933.Sakai.megablast` | 5718 | 1424 | 993 | 843 |
| `Sakai.MG1655.megablast` | 6476 | 6476 | 6476 | 0 |

Task-blastn interpretation:

- Same-coordinate HSPs have matching bit scores.
- The primary task-blastn defect is survivor-set divergence: which HSPs are
  created, trimmed, reevaluated, purged, or retained.
- NCBI-only and LOSAT-only rows are concentrated mostly in low score bins
  around 30-49 bits, with a small number of higher-score boundary/merge cases.
- E-values differ for many same-coordinate rows despite identical bit scores,
  so E-value/statistics parity is a separate confirmed defect.

## Likely Root Causes

### 1. Survivor-set divergence after gapped extension

The first task-blastn baseline has zero same-coordinate bit-score differences,
so raw scoring for shared HSPs is probably not the first-order issue.

The leading `NZ_CP006932` self difference shows the same bit score for a
different boundary:

- NCBI-only: `q=625997-626608 s=239620-239018 bit=89.7`
- LOSAT-only: `q=626077-626608 s=239543-239018 bit=89.7`

This points to start/end selection, traceback boundary handling, common-endpoint
trim/delete behavior, or interval-tree containment rather than the bit-score
formula itself.

Primary code to inspect:

- `LOSAT/src/algorithm/blastn/alignment/gapped.rs`
- `LOSAT/src/algorithm/blastn/filtering/purge_endpoints.rs`
- `LOSAT/src/algorithm/blastn/interval_tree.rs`
- `LOSAT/src/algorithm/blastn/blast_engine/run.rs`

### 2. E-value path mismatch

NCBI final HSP E-values use the per-context effective search space:

- NCBI: `blast_hits.c:1887-1890`

LOSAT currently has two visible paths:

- initial final-HSP build uses database-search recalculation in
  `LOSAT/src/algorithm/blastn/blast_engine/run.rs:8645-8651`
- post-reevaluation uses `query_eff_searchsp` in
  `LOSAT/src/algorithm/blastn/blast_engine/run.rs:8804-8817`

This can explain same-coordinate rows with identical bit scores but different
E-values. The fix should make final HSP E-value calculation use NCBI's
per-context `eff_searchsp` path consistently.

### 3. Megablast is a separate track

Megablast has cases with large same-coordinate bit-score and length differences.
Do not use megablast output to diagnose task-blastn DP/purge defects until
task-blastn parity is under control.

## Fix Plan

### Phase 0: Preserve the Baseline

1. Keep `tests/blastn_parity_manifest.tsv` as the canonical task/options
   manifest.
2. Split any new output files by task name:
   - `*.task_blastn.*`
   - `*.megablast.*`
3. Regenerate only the specific case under investigation, not the full suite,
   unless a parity sweep has been completed.
4. Store the comparison command and output in a scratch or docs-linked artifact
   before changing code.

Acceptance:

- Every BLASTN discussion names task, options, query, subject, and oracle file.
- The current task-blastn table above is reproducible.

### Phase 1: Trace One NCBI-only and One LOSAT-only HSP

Use one high-signal pair first:

- NCBI-only target:
  `NZ_CP006932.1 NZ_CP006932.1 q=625997-626608 s=239620-239018 bit=89.7`
- LOSAT-only target:
  `NZ_CP006932.1 NZ_CP006932.1 q=626077-626608 s=239543-239018 bit=89.7`

Trace stages:

1. lookup hit retrieval
2. ungapped extension
3. preliminary gapped score-only extension
4. preliminary common-endpoint purge
5. traceback start selection
6. final traceback
7. common-endpoint trim pass
8. ambiguity reevaluation
9. identity/length test
10. common-endpoint delete pass
11. score resort
12. final interval-tree containment
13. final hit-list pruning/output order

Acceptance:

- For the NCBI-only target, LOSAT reports the first stage where the expected HSP
  is absent, changed, or deleted.
- For the LOSAT-only target, LOSAT reports the first stage where NCBI would
  reject, trim, or merge it.

### Phase 2: Fix E-value Consistency

Make all final HSP E-value assignments follow NCBI:

```c
hsp->evalue = BLAST_KarlinStoE_simple(score, kbp[kbp_context],
                                      query_info->contexts[hsp->context].eff_searchsp);
```

Implementation direction:

1. Add or reuse a helper that takes raw score, Karlin params, and per-context
   `eff_searchsp`.
2. Replace final-HSP creation E-value calculation in
   `blast_engine/run.rs:8645-8651` with the per-context path.
3. Keep bit-score calculation unchanged unless NCBI source proves a mismatch.
4. Add focused tests comparing same-coordinate E-values for task-blastn
   fixtures.

Acceptance:

- Same-coordinate E-value diffs drop to zero for the three task-blastn cases,
  or any remaining diffs are traced to formatting rather than numeric value.
- Same-coordinate bit-score diffs remain zero.

### Phase 3: Align Traceback Boundary Selection

Focus on NCBI:

- `BlastGetStartForGappedAlignmentNucl`
- `Blast_SemiGappedAlign`
- `AdjustSubjectRange`
- `Blast_HSPAdjustSubjectOffset`

Implementation direction:

1. Compare LOSAT traceback start offsets against NCBI for the selected targets.
2. Verify whether LOSAT starts from seed midpoint, HSP gapped start, or
   recalculated offsets at the same timing as NCBI.
3. Verify subject range adjustment and start-shift handling before final
   traceback.
4. Patch only after the exact NCBI line-level behavior is identified.

Acceptance:

- The selected 89.7-bit boundary case has the same final coordinates as NCBI.
- No new same-coordinate bit-score diffs are introduced.

### Phase 4: Align Common-Endpoint and Containment Purges

Focus on NCBI:

- `blast_traceback.c:633-692`
- `blast_hits.c:2268-2535`
- `blast_itree.c:797-847`
- `blast_itree.c:930-995`

Implementation direction:

1. Instrument sort keys and purge decisions for the selected target HSPs.
2. Confirm pass 1 trim behavior with `purge=false`.
3. Confirm pass 2 delete behavior with `purge=true`.
4. Confirm final `ScoreCompareHSPs` ordering before interval containment.
5. Confirm `min_diag_separation` and query-context offset handling.

Acceptance:

- NCBI-only and LOSAT-only counts decrease without changing shared HSP scores.
- The three task-blastn cases reach identical hit counts, then identical final
  output.

### Phase 5: Re-run Parity Sweep

Run only after Phases 2-4 have no known stage-level divergence for selected
targets.

Commands:

```bash
cd LOSAT
python tests/compare_blastn_parity.py --manifest tests/blastn_parity_manifest.tsv --limit 10
```

Then, if task-blastn is clean:

```bash
cd LOSAT/tests
bash run_comparison.sh
```

Acceptance:

- Task-blastn cases have zero coordinate-key differences.
- Same-coordinate bit score, E-value, identity, mismatch, gapopen, and length
  differences are zero.
- Megablast failures are tracked separately after task-blastn is stable.

## Do Not Do

- Do not compare task-blastn against NCBI default megablast output.
- Do not tune thresholds to reduce hit-count deltas.
- Do not accept hit-count deltas as success.
- Do not change score formulas unless the corresponding NCBI source line proves
  LOSAT is wrong.
- Do not delegate unsupported behavior to NCBI BLAST at runtime.
