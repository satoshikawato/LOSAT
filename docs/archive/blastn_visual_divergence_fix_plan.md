# BLASTN Visual Divergence Isolation and Fix Plan

Date: 2026-05-17

Goal: identify and fix the LOSATN vs NCBI BLASTN task divergence visible in the
BLASTN comparison plots, especially the accumulated-length and identity
distribution gaps for long HSPs. NCBI BLAST is the only behavioral oracle.
NCBI binaries may be used only for diagnostics and oracle output generation;
LOSAT runtime code must remain a pure Rust implementation.

This document is a focused execution plan for the plotted BLASTN cases. It
extends, but does not replace, the broader task-blastn remediation history in
`docs/losatn_blastn_bitwise_parity_remediation_plan.md`.

## Plot-Based Diagnosis

The plots show a real parity defect, not just a rendering artifact.

- `PeseMJNV.PemoMJNVB`: BLAST+ has a strong accumulated-length spike around
  multi-kb HSPs that LOSATN does not reproduce. This is high priority because a
  small number of long HSPs can dominate the accumulated-length histogram.
- `PmeNMV.MjPMNV`: bit score vs length largely overlaps, but the identity
  histogram is much narrower in LOSATN. This points away from the raw scoring
  formula and toward HSP survival, trimming, or edit-script/statistic shape.
- `PmeNMV.PesePMNV`: LOSATN and BLAST+ share the same broad score/length curve,
  but some long-HSP bins and identity bands differ. This is a good boundary and
  post-traceback pruning stress case.
- The overall BLASTN row shows broad overlap, so the likely defect is not a
  global score/statistics failure. It is probably in one or more stage-specific
  survivor decisions: preliminary gapped extension, traceback start selection,
  common-endpoint purge, ambiguity re-evaluation, interval-tree containment, or
  final hit-list pruning.

Current high-priority task-blastn notes already show zero shared-coordinate
bit-score and E-value diffs after recent fixes, with survivor-set divergence
remaining. Treat the new plotted cases as broader coverage for the same
survivor-set class until a stage trace proves otherwise.

## Fixed Oracle Conditions

All comparisons in this plan must use NCBI `blastn -task blastn`, not NCBI's
default megablast mode.

Use these fixed options unless a manifest row explicitly records a different
NCBI source-backed option:

- task: `blastn`
- reward: `2`
- penalty: `-3`
- gap open: `5`
- gap extend: `2`
- word size: `11`
- DUST: enabled
- soft masking: enabled
- strand: both
- E-value: `10`
- max target seqs: `500`
- max HSPs per subject: `0`
- threads: `1`
- output: outfmt `7`

The current plotted BLASTN cases are defined in:

- `LOSAT/tests/plot_overall_trend.py:58-66`
- `LOSAT/tests/plot_comparison.py:182-190`
- `LOSAT/tests/run_comparison.sh:31-56` for LOSATN commands
- `LOSAT/tests/run_comparison.sh:281-306` for NCBI BLASTN oracle commands

Promote these visual cases into `LOSAT/tests/blastn_parity_manifest.tsv` before
making behavioral changes:

| Case id | Query | Subject | NCBI output | LOSATN output |
| --- | --- | --- | --- | --- |
| `PesePMNV.MjPMNV.task_blastn` | `tests/fasta/AP027152.fasta` | `tests/fasta/AP027202.fasta` | `tests/blast_out/PesePMNV.MjPMNV.blastn.out` | `tests/losat_out/PesePMNV.MjPMNV.losatn.blastn.out` |
| `PmeNMV.MjPMNV.task_blastn` | `tests/fasta/LC738869.fasta` | `tests/fasta/AP027202.fasta` | `tests/blast_out/PmeNMV.MjPMNV.blastn.out` | `tests/losat_out/PmeNMV.MjPMNV.losatn.blastn.out` |
| `PmeNMV.PesePMNV.task_blastn` | `tests/fasta/LC738869.fasta` | `tests/fasta/AP027152.fasta` | `tests/blast_out/PmeNMV.PesePMNV.blastn.out` | `tests/losat_out/PmeNMV.PesePMNV.losatn.blastn.out` |
| `PeseMJNV.PemoMJNVB.task_blastn` | `tests/fasta/LC738873.fasta` | `tests/fasta/LC738871.fasta` | `tests/blast_out/PeseMJNV.PemoMJNVB.blastn.out` | `tests/losat_out/PeseMJNV.PemoMJNVB.losatn.blastn.out` |
| `PemoMJNVA.PeseMJNV.task_blastn` | `tests/fasta/LC738870.fasta` | `tests/fasta/LC738873.fasta` | `tests/blast_out/PemoMJNVA.PeseMJNV.blastn.out` | `tests/losat_out/PemoMJNVA.PeseMJNV.losatn.blastn.out` |
| `MelaMJNV.PemoMJNVA.task_blastn` | `tests/fasta/LC738874.fasta` | `tests/fasta/LC738870.fasta` | `tests/blast_out/MelaMJNV.PemoMJNVA.blastn.out` | `tests/losat_out/MelaMJNV.PemoMJNVA.losatn.blastn.out` |
| `MjeNMV.MelaMJNV.task_blastn` | `tests/fasta/LC738868.fasta` | `tests/fasta/LC738874.fasta` | `tests/blast_out/MjeNMV.MelaMJNV.blastn.out` | `tests/losat_out/MjeNMV.MelaMJNV.losatn.blastn.out` |
| `MjPMNV.MlPMNV.task_blastn` | `tests/fasta/AP027202.fasta` | `tests/fasta/LC738875.fasta` | `tests/blast_out/MjPMNV.MlPMNV.blastn.out` | `tests/losat_out/MjPMNV.MlPMNV.losatn.blastn.out` |

Keep the existing `NZ_CP006932.NZ_CP006932.task_blastn` self case as a dense
stress case, but do not let it hide visual-case regressions.

## Primary NCBI References

Read these NCBI source ranges before editing any Rust implementation. Every
Rust behavior change must include an immediately adjacent NCBI reference comment
with file path and line numbers.

- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3810-4105`
  - preliminary gapped HSP loop, interval-tree precheck, nucleotide seed shift,
    cutoff check, HSP creation, and immediate interval-tree insertion
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:736-962`
  - `Blast_SemiGappedAlign` score-only DP
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3323-3390`
  - `BlastGetStartForGappedAlignmentNucl`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:364-733`
  - `ALIGN_EX` traceback DP and edit-script row behavior
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:350-613`
  - traceback loop, containment precheck, start selection, subject adjustment,
    final traceback, identity/length test, and tree insertion timing
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:633-692`
  - post-traceback common-endpoint purge, ambiguity re-evaluation, second
    BLASTN purge, score resort, and final interval-tree containment
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1109-1132`
  - `Blast_HSPGetAdjustedOffsets`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1330-1354`
  - `ScoreCompareHSPs`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2268-2535`
  - `s_QueryOffsetCompareHSPs`, `s_QueryEndCompareHSPs`, `s_CutOffGapEditScript`,
    and `Blast_HSPListPurgeHSPsWithCommonEndpoints`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_itree.c:797-847`
  - `s_HSPIsContained`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_itree.c:930-995`
  - `BlastIntervalTreeContainsHSP`

Primary LOSATN files:

- `LOSAT/src/algorithm/blastn/blast_engine/run.rs:8480-9668`
- `LOSAT/src/algorithm/blastn/alignment/gapped.rs:1834-3312`
- `LOSAT/src/algorithm/blastn/filtering/purge_endpoints.rs:132-809`
- `LOSAT/src/algorithm/blastn/interval_tree.rs:895-1108`
- `LOSAT/src/algorithm/blastn/hsp.rs:302-527`
- `LOSAT/src/algorithm/blastn/tracing.rs:1-380`
- `LOSAT/tests/compare_blastn_parity.py`
- `LOSAT/tests/blastn_parity_manifest.tsv`
- `LOSAT/tests/blastn_trace_targets.tsv`

## Phase 0: Refresh and Rank the Visual Cases

Objective: convert the plotted divergence into reproducible coordinate-key
failures, then rank by accumulated-length impact.

Tasks:

1. Add the visual BLASTN cases above to `tests/blastn_parity_manifest.tsv`.
2. Regenerate LOSATN output for those cases only with `--refresh-losat`.
3. Regenerate NCBI oracle output only if the current `tests/blast_out/*.out`
   files do not match the fixed oracle conditions above.
4. Run `tests/compare_blastn_parity.py` for the visual cases.
5. Extend the comparer if necessary to report:
   - accumulated alignment length from NCBI-only rows
   - accumulated alignment length from LOSAT-only rows
   - top N NCBI-only and LOSAT-only rows by `length * bitscore`
   - score bins, identity bins, and length bins for exclusive rows
6. Store the first ranked trace targets in `tests/blastn_trace_targets.tsv`.

Acceptance:

- Each plotted BLASTN case has a manifest row.
- Each case has fresh NCBI-vs-LOSATN counts for common coordinates, NCBI-only,
  LOSAT-only, shared bit-score diffs, shared E-value diffs, shared identity
  diffs, and accumulated-length impact.
- The first patch target is the highest-impact exclusive HSP or HSP cluster, not
  an arbitrary low-score row.

## Phase 1: Find the First Divergent Stage

Objective: determine where one high-impact NCBI-only HSP disappears and where
one high-impact LOSAT-only HSP first survives contrary to NCBI.

Use existing LOSATN trace controls:

- `LOSAT_TRACE_BLASTN_HSP="qstart,qend,sstart,send"`
- `LOSAT_TRACE_BLASTN_SEED="q,s"`
- `LOSAT_TRACE_BLASTN_CONTEXT=<context_idx>`
- `LOSAT_TRACE_BLASTN_SUBJECT=<subject_id_or_index>`
- `LOSAT_TRACE_BLASTN_STAGE=<seed|ungapped|prelim|traceback|purge|hitlist|all>`

Trace these checkpoints in order:

1. lookup hit and ungapped extension
2. preliminary interval-tree containment precheck
3. preliminary score-only gapped extension
4. preliminary cutoff acceptance and immediate tree insertion
5. preliminary common-endpoint purge before traceback
6. traceback interval-tree precheck
7. `BlastGetStartForGappedAlignmentNucl` result
8. `AdjustSubjectRange` start shift and adjusted subject slice
9. left/right `ALIGN_EX` traceback operations
10. final edit-script merge and HSP coordinate construction
11. identity/length test
12. common-endpoint purge pass 1
13. ambiguity re-evaluation
14. BLASTN common-endpoint purge pass 2
15. score resort
16. final interval-tree containment
17. final E-value reap, hit-list pruning, and output ordering

Acceptance:

- For each selected NCBI-only row, record the first stage where LOSATN is
  absent, shortened, shifted, contained, purged, or sorted behind another HSP.
- For each selected LOSAT-only row, record the first stage where NCBI would have
  rejected, trimmed, contained, or deleted it.
- No behavioral patch is made before this first divergent stage is identified
  from NCBI source and LOSATN trace data.

## Phase 2: Patch by Divergence Class

Apply exactly one divergence-class fix at a time.

### Class A: Preliminary Gapped Extension

Patch here if divergence appears before final traceback.

Audit:

- query/subject context offsets and canonical internal subject offsets
- interval-tree initialization bounds and `min_diag_separation`
- nucleotide seed shift by three bases when the 8-base seed condition holds
- `Blast_SemiGappedAlign` x-drop behavior and score-only traceback absence
- immediate tree insertion order after preliminary HSP acceptance
- preliminary common-endpoint purge with BLASTN `purge=TRUE`

Rust target files:

- `LOSAT/src/algorithm/blastn/blast_engine/run.rs`
- `LOSAT/src/algorithm/blastn/alignment/gapped.rs`
- `LOSAT/src/algorithm/blastn/interval_tree.rs`

NCBI anchors:

- `blast_gapalign.c:3810-4105`
- `blast_gapalign.c:736-962`
- `blast_engine.c:545`

### Class B: Traceback Boundary or Edit-Script Shape

Patch here if coordinates match broadly but lengths, identity, mismatch,
gap-open count, or long-HSP boundaries differ.

Audit:

- `BlastGetStartForGappedAlignmentNucl` start movement
- `AdjustSubjectRange` and `Blast_HSPAdjustSubjectOffset`
- left and right `ALIGN_EX` direction-specific tie behavior
- traceback state-array row reuse and unwritten-cell semantics
- `Blast_PrelimEditBlockToGapEditScript` merge order
- stats recomputation from `query_nomask` after final edit script construction

Rust target files:

- `LOSAT/src/algorithm/blastn/alignment/gapped.rs`
- `LOSAT/src/algorithm/blastn/blast_engine/run.rs`

NCBI anchors:

- `blast_gapalign.c:3323-3390`
- `blast_gapalign.c:364-733`
- `blast_traceback.c:436-600`

### Class C: Post-Traceback Purge and Containment

Patch here if LOSATN constructs the same HSP as NCBI but deletes, keeps, trims,
or orders it differently after traceback.

Audit:

- `s_QueryOffsetCompareHSPs` and `s_QueryEndCompareHSPs` exact tie-breakers
- `purge=FALSE` pass 1 trim-vs-delete behavior for BLASTN
- `s_CutOffGapEditScript` boundary cutting and moved-to-end semantics
- ambiguity re-evaluation score and identity/length test timing
- BLASTN-only `purge=TRUE` pass 2
- score resort with `ScoreCompareHSPs`
- final interval-tree reset, contains, and add order

Rust target files:

- `LOSAT/src/algorithm/blastn/filtering/purge_endpoints.rs`
- `LOSAT/src/algorithm/blastn/interval_tree.rs`
- `LOSAT/src/algorithm/blastn/hsp.rs`
- `LOSAT/src/algorithm/blastn/blast_engine/run.rs`

NCBI anchors:

- `blast_traceback.c:633-692`
- `blast_hits.c:2268-2535`
- `blast_hits.c:1330-1354`
- `blast_itree.c:797-995`

### Class D: Final Hit-List Pruning or Output Coordinates

Patch here only if internal survivor sets match but final outfmt rows differ.

Audit:

- `Blast_HSPGetAdjustedOffsets` minus-strand coordinate adjustment
- use of internal contexted offsets for pruning before output conversion
- `max_hsps_per_subject` and `max_target_seqs` order
- final row sort and outfmt 7 formatting

Rust target files:

- `LOSAT/src/algorithm/blastn/hsp.rs`
- `LOSAT/src/algorithm/blastn/blast_engine/run.rs`
- `LOSAT/src/report/`

NCBI anchors:

- `blast_hits.c:1109-1132`
- `blast_hits.c:1330-1354`
- `blast_traceback.c:234-250`

## Phase 3: Verification Gates

Run tests only after a complete patch for the identified divergence class.

Targeted gate:

```bash
cd LOSAT
cargo check
python tests/compare_blastn_parity.py \
  --manifest tests/blastn_parity_manifest.tsv \
  --case-id PeseMJNV.PemoMJNVB.task_blastn \
  --case-id PmeNMV.MjPMNV.task_blastn \
  --case-id PmeNMV.PesePMNV.task_blastn \
  --refresh-losat \
  --limit 10
```

Existing high-priority gate:

```bash
cd LOSAT
python tests/compare_blastn_parity.py \
  --manifest tests/blastn_parity_manifest.tsv \
  --case-id MelaMJNV.PemoMJNVA.task_blastn \
  --case-id PemoMJNVA.PeseMJNV.task_blastn \
  --case-id NZ_CP006932.NZ_CP006932.task_blastn \
  --refresh-losat \
  --limit 10
```

Full task-blastn gate after targeted deltas shrink:

```bash
cd LOSAT
python tests/compare_blastn_parity.py \
  --manifest tests/blastn_parity_manifest.tsv \
  --refresh-losat \
  --fail-on-diff
```

Acceptance:

- Shared-coordinate bit-score and E-value diffs remain zero.
- Shared-coordinate identity, length, mismatch, and gap-open diffs are zero.
- NCBI-only and LOSAT-only coordinate sets shrink to zero for all task-blastn
  manifest cases.
- Final output row order and text match NCBI exactly.
- The plots can be regenerated and no longer show BLASTN accumulated-length or
  identity distribution divergence.

## Implementation Rules

- Do not add features not present in NCBI BLAST.
- Do not call NCBI BLAST from LOSAT runtime, build scripts, feature paths, or
  fallback paths.
- Do not patch from visual intuition alone. Every Rust behavior change must be
  tied to a specific NCBI source line range.
- Keep diagnostic-only tracing behind environment variables.
- If a suspected behavior cannot be found in NCBI source, do not implement it.
- If several differences remain, patch all offenders in the same verified
  divergence class, then fix compile issues and run the relevant gate.

