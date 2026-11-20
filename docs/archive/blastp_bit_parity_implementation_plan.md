# BLASTP Bit-Parity Implementation Plan

Date: 2026-05-12

Goal: make LOSATP output bit-perfect with NCBI BLASTP for the supported pure-Rust
scope, without adding any runtime dependency on NCBI BLAST.

Supported scope for this plan:

- `blastp`
- BLOSUM62
- gap open 11, gap extend 1
- `word_size=3`
- `comp_based_stats=2`
- `use_sw_tback=false`
- default SEG behavior used by LOSATP/NCBI advanced BLASTP options

NCBI BLAST remains the oracle only for source inspection and comparison runs.
Do not call NCBI BLAST from LOSAT runtime code, build scripts, fallback paths, or
unsupported feature paths.

## Current Read

The large apparent hit-count gap is mostly a comparison-condition issue.

`LOSAT/tests/losat_out/*.losatp.out` contains default-style outputs
(`evalue=10`, `max_target_seqs=500`). The active comparison script
`LOSAT/tests/run_blastp_comparison.sh` runs with `-evalue 1e-5` and
`-max_target_seqs 10`. Comparing those two classes makes LOSATP look like it has
thousands of extra hits, but that is not a valid parity comparison.

With default-to-default outputs, remaining differences are small:

- a few missing/extra HSPs near ranking/e-value thresholds
- same-coordinate HSPs with different identity/length/gap counts
- occasional same-coordinate HSPs with different score/e-value

That points away from seed generation as the first target and toward the
composition-based statistics/Kappa postprocess, adjusted-matrix traceback, HSP
re-evaluation, identity recomputation, and final pruning.

## Primary NCBI References

Use these as source of truth before writing Rust changes:

- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:3900-4008`
  - preliminary gapped alignment loop, restricted retry, interval tree rebuild
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4209-4319`
  - `s_BlastProtGappedAlignment`, score-only preliminary gapped extension
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_gapalign.c:4549-4715`
  - `BLAST_GappedAlignmentWithTraceback`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:459-535`
  - `s_ComputeNumIdentities`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_kappa.c:3658-3716`
  - Kappa postprocess order: distinct alignments, contained purge, e-value
    evaluation, score normalization, identity recomputation, heap insert
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:635-688`
  - final traceback re-evaluation, identity/length test, score sort,
    interval-tree containment purge
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:440-642`
  - `Blast_HSPReevaluateWithAmbiguitiesGapped`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:746-824`
  - `s_Blast_HSPGetNumIdentitiesAndPositives`
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1330-1455`
  - score/e-value comparators and HSP sorting
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2455-2535`
  - common-endpoint purge

Rust areas to change only after checking the matching NCBI source:

- `LOSAT/src/algorithm/blastp/blast_engine.rs`
- `LOSAT/src/algorithm/blastp/gapalign.rs`
- `LOSAT/src/algorithm/blastp/kappa.rs`
- `LOSAT/src/algorithm/blastp/hsp.rs`
- `LOSAT/src/algorithm/blastp/alignment.rs`

## Phase 0: Freeze Valid Comparison Sets

Objective: prevent false positives from mixed command-line conditions.

Tasks:

1. Create a small comparison manifest for BLASTP with explicit options:
   `task`, `matrix`, `gap_open`, `gap_extend`, `word_size`, `threshold`,
   `window_size`, `comp_based_stats`, `seg`, `evalue`, `max_target_seqs`,
   `max_hsps_per_subject`, and `num_threads`.
2. Separate two oracle classes:
   - default BLASTP parity: `evalue=10`, `max_target_seqs=500`
   - strict script parity: `evalue=1e-5`, `max_target_seqs=10`
3. For every reported delta, record:
   - NCBI-only coordinate keys
   - LOSAT-only coordinate keys
   - same-coordinate value differences
   - first divergent subject/query pair

Acceptance:

- No parity claim is made unless NCBI and LOSAT were run with identical options.
- Each failing case has one first divergent HSP selected for instrumentation.

## Phase 1: Add Targeted Differential Tracing

Objective: locate the exact first stage where one HSP diverges.

Add tracing gated behind a debug environment variable, for example
`LOSAT_DEBUG_BLASTP_TRACE_HSP="query_id,subject_id,qstart,qend,sstart,send"`.

Trace points:

1. Initial ungapped HSP saved from aa scan.
2. Chained initial HSP list after `chain_blastp_init_hsps`.
3. Preliminary score-only gapped HSP:
   - ungapped seed
   - `gapped_query_start`
   - `gapped_subject_start`
   - restricted vs exact mode
   - raw preliminary score
4. Kappa redo input alignment list.
5. Composition-adjusted matrix metadata:
   - adjustment rule
   - scaled lambda
   - score divisor
   - range begin/end for query and subject
6. Traceback result:
   - edit script
   - query/subject offsets
   - raw score before normalization
7. Postprocess:
   - contained purge decision
   - Spouge e-value before and after composition p-value adjustment
   - normalized score and bit score
   - identity/positives/gaps after recomputation
8. Final hit-list pruning:
   - HSP sort key
   - subject-list sort key
   - heap insertion/replacement decision

Acceptance:

- For one divergent HSP, LOSAT logs show whether the first mismatch occurs in
  preliminary gapped scoring, Kappa redo, identity recomputation, or final
  pruning.
- The tracing code is inactive unless the environment variable is set.
- Every added trace block includes NCBI reference comments with source path and
  line numbers.

## Phase 2: Fix Identity and Alignment-View Differences

Objective: eliminate same-coordinate rows where only `%identity`, mismatch,
gap-open, or alignment length differs.

Primary Rust targets:

- `LOSAT/src/algorithm/blastp/alignment.rs`
- `LOSAT/src/algorithm/blastp/kappa.rs`
- `LOSAT/src/algorithm/blastp/gapalign.rs`

NCBI behavior to mirror:

- `blast_kappa.c:459-535` calls `Blast_HSPGetNumIdentitiesAndPositives` with
  `query_blk->sequence_nomask` and `subject_blk->sequence`.
- `blast_hits.c:746-824` counts identities, positives, alignment length, and
  gap operations by walking `GapEditScript` exactly.
- For protein BLAST, positives use `matrix[*q][*s] > 0`.

Known example to reproduce:

- `WSSV.PajaWSV`
- same coordinate:
  `YP_009220597.1 / BDT63580.1 / q 13..6074 / s 14..6057`
- NCBI `%identity=54.458`
- LOSAT `%identity=54.409`
- same bit score and e-value, so the score path is not the first issue here.

Implementation steps:

1. Add a focused regression test that reconstructs this HSP's final edit script
   and verifies identity/positives/gap counts against the NCBI pairwise output.
2. Compare LOSAT's `GapEditOp` interpretation with NCBI `GapEditScript` op
   semantics for `eGapAlignSub`, `eGapAlignDel`, and `eGapAlignIns`.
3. Verify that query sequence for identity uses unmasked query, while traceback
   scoring uses the correct masked/adjusted sequence according to the NCBI
   postprocess step.
4. Patch only the first confirmed mismatch.

Acceptance:

- Same-coordinate identity-only diffs disappear.
- Existing exact matches remain unchanged.

## Phase 3: Fix CBS/Kappa Redo Traceback Differences

Objective: eliminate same-coordinate rows where score/e-value or alignment span
changes after composition-based postprocess.

Primary Rust targets:

- `LOSAT/src/algorithm/blastp/kappa.rs`
- `LOSAT/src/algorithm/blastp/gapalign.rs`
- `LOSAT/src/core/composition_adjustment/redo_alignment.rs`
- `LOSAT/src/core/composition_adjustment/adjust_scores.rs`

NCBI behavior to mirror:

- `blast_kappa.c:3658-3716` postprocess order is strict:
  `s_HSPListFromDistinctAlignments`, contained purge,
  `s_HitlistEvaluateAndPurge`, score normalization,
  `s_ComputeNumIdentities`, heap insert.
- `blast_gapalign.c:4549-4715` does traceback with:
  - left extension including the seed point
  - right extension excluding the seed point
  - `fwd_prelim_tback` and `rev_prelim_tback`
  - terminal gap pruning after converting preliminary edit blocks
- `blast_gapalign.c:4549-4715` uses the current score matrix in
  `score_params`, which under CBS means the adjusted matrix, not the original
  BLOSUM62 matrix.

Known existing diagnostic:

- `LOSAT/src/algorithm/blastp/gapalign.rs` has an ignored diagnostic test for
  residual `WSSV/PajaWSV CBS=2 traceback diff`. Use that as the first concrete
  debugger entry, but verify the current NCBI source before changing logic.

Implementation steps:

1. Trace one same-coordinate score-different case through Kappa redo.
2. Compare adjusted matrix generation with NCBI for that query/subject range:
   - range begin/end
   - SEG masking decision for subject ranges
   - p-value and lambda ratio
   - score divisor/local scaling factor
3. Compare traceback seeds:
   - incoming preliminary HSP start/end
   - range-relative `q_start` and `s_start`
   - `gap_align->query_start/stop`
   - `gap_align->subject_start/stop`
4. Compare edit scripts before and after terminal-gap pruning.
5. Patch the first mismatch in Rust with NCBI comments immediately above the
   modified logic.

Acceptance:

- Same-coordinate score/e-value diffs disappear for the selected case.
- The fix reduces, not reshuffles, the aggregate divergent-HSP list.

## Phase 4: Fix Borderline Missing/Extra HSPs

Objective: remove the remaining one-to-few HSP count deltas.

Primary Rust targets:

- `LOSAT/src/algorithm/blastp/hsp.rs`
- `LOSAT/src/algorithm/blastp/blast_engine.rs`
- `LOSAT/src/algorithm/blastp/kappa.rs`

NCBI behavior to mirror:

- `blast_hits.c:2455-2535`: common-endpoint purge.
- `blast_traceback.c:635-688`: after re-evaluation, sort by score again and
  remove contained HSPs through interval tree containment.
- `blast_hits.c:1330-1455`: score/e-value comparators and tie-breakers.
- `blast_kappa.c:3658-3716`: heap insert uses best e-value, best score, and
  subject index after Kappa postprocess.

Implementation steps:

1. For each missing/extra HSP, classify:
   - never generated as preliminary
   - generated then lost in restricted/exact preliminary scoring
   - lost in Kappa redo/e-value purge
   - lost in contained/common-endpoint purge
   - lost in hit-list heap or `max_target_seqs` pruning
2. Fix in that order. Do not alter final sort/prune to hide upstream
   divergence.
3. Add comparator regression tests for any tie-breaker mismatch.

Acceptance:

- Per-case line counts match for default BLASTP outputs.
- Coordinate-key sets match before checking formatting.

## Phase 5: Output Formatting Pass

Objective: after semantic parity, make text output bit-perfect.

Targets:

- tabular outfmt 6
- outfmt 7 header/body if enabled
- pairwise outfmt 0 for BLASTP

Checks:

- e-value formatting
- bit-score formatting
- percent identity precision
- HSP order within subject
- subject-list order within query
- empty-query and no-hit sections

Acceptance:

- `diff -u` is empty for every frozen comparison case in the supported scope.

## Phase 6: Verification Discipline

Run tests only after a complete patch sweep for a known divergence.

Minimum verification after each patch sweep:

```bash
cd LOSAT && cargo check
cd LOSAT/tests && bash run_blastp_comparison.sh
```

If the patch touches shared BLASTP/Kappa/gapalign primitives, also run targeted
unit tests around the changed module before full comparison.

Do not run broad "maybe helpful" tests while a known NCBI divergence is still
unclassified.

## Implementation Rules

- Read NCBI C/C++ source before every Rust logic change.
- Add NCBI reference comments with file path and line numbers immediately above
  every changed Rust logic block.
- Do not add unsupported BLASTP modes to make a comparison pass.
- Do not delegate any LOSAT runtime behavior to NCBI binaries or libraries.
- Preserve NCBI timing/order even where Rust could express the logic more
  cleanly.
- Prefer one confirmed discrepancy per patch sweep, then re-run the focused
  comparison.

## Suggested First Patch Sweep

Start with identity-only drift because it is narrower than traceback scoring.

1. Reproduce `WSSV.PajaWSV` same-coordinate identity diff.
2. Add a focused regression around the final `GapEditOp` walk.
3. Compare against NCBI `blast_kappa.c:459-535` and `blast_hits.c:746-824`.
4. Patch `alignment.rs`, `kappa.rs`, or `gapalign.rs` only if the exact NCBI
   mismatch is identified.
5. Re-check the full default BLASTP output set for same-coordinate value diffs.

Then move to CBS traceback score drift using the existing ignored diagnostic in
`gapalign.rs`.
