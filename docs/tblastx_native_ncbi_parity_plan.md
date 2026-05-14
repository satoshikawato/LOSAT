# TBLASTX Native NCBI Full Parity Plan

Date: 2026-05-14

Goal: close the remaining native LOSAT vs NCBI BLAST+ TBLASTX parity gap for
the `LC738874.fasta` vs `LC738875.fasta` fixture, then generalize the fix to
the long-sequence TBLASTX backlog without regressing native/Wasm determinism.

NCBI BLAST is the source of truth. NCBI binaries may be used only as comparison
oracles during diagnostics and validation. LOSAT runtime code must not call,
link, embed, or fall back to NCBI BLAST.

## Current Status

Native/Wasm determinism is already fixed for this fixture. Treat that as a
separate closed issue.

Fixture:

```text
tblastx -q LOSAT/tests/fasta/LC738874.fasta \
        -s LOSAT/tests/fasta/LC738875.fasta \
        --outfmt 6 -n 1
```

Latest artifact directory:

```text
.tmp/tblastx_wasm_parity_codex_run_abs
```

Observed counts from that artifact:

```text
NCBI BLAST+ oracle: 4319 hits
LOSAT native:       4350 hits
LOSAT Wasm SIMD:    4350 hits
LOSAT Wasm scalar:  4350 hits
native/Wasm sorted diff lines: 0
NCBI/native sorted diff lines: 3801
```

Additional classification from the same raw outputs:

```text
exact full-line matches after sorting:                 2434
matches ignoring only e-value text:                    4276
NCBI-only HSP keys after ignoring e-value:               43
LOSAT-only HSP keys after ignoring e-value:              74
all NCBI-only and LOSAT-only HSP-key deltas:        <23 bit
```

Interpretation:

- The remaining defect is native LOSAT vs NCBI TBLASTX parity, not Wasm SIMD.
- Most shared HSPs have the same coordinates, identity fields, raw-derived bit
  score text, and differ only in printed E-values.
- The true hit-set difference is small and concentrated at the low-score
  reporting boundary.
- The first priority is therefore E-value and preliminary reaping parity around
  TBLASTX sum statistics, not another SIMD investigation.

## Acceptance Criteria

For the fixture above:

1. Raw outfmt 6 output is byte-identical to the NCBI oracle, including order.
2. Sorted outfmt 6 output is byte-identical to the NCBI oracle.
3. Hit count is exactly `4319`.
4. Shared-HSP E-values match NCBI at the printed text level and, where traced,
   at the same underlying double-producing formula path.
5. Native, Wasm SIMD, and Wasm scalar LOSAT remain byte-identical to each other.

For follow-up TBLASTX backlog cases:

1. The fix improves or preserves hit-count parity on long-sequence gencode 4
   cases.
2. No performance optimization is mixed into the parity patch.
3. Every Rust code change includes an adjacent NCBI source comment with file
   path and line numbers.

## Primary NCBI References

Read these before changing code in the corresponding area.

- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1482-1541`
  - Ungapped TBLASTX post-search order:
    `Blast_HSPListReevaluateUngapped`,
    `BLAST_LinkHsps` or `Blast_HSPListGetEvalues`,
    `s_Blast_HSPListReapByPrelimEvalue`,
    `Blast_HSPListReapByQueryCoverage`,
    `Blast_HSPListGetBitScores`.
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:643-674`
  - Preliminary E-value reaping uses `hit_params->prelim_evalue`.
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:675-733`
  - `Blast_HSPReevaluateWithAmbiguitiesUngapped`.
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1811-1905`
  - `Blast_HSPListGetEvalues`, including context K/Lambda selection and
    effective search space.
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1907-1927`
  - `Blast_HSPListGetBitScores`.
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2010-2046`
  - Query coverage reaping.
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2455-2535`
  - Common-endpoint purge.
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2609-2673`
  - `Blast_HSPListReevaluateUngapped` list-level translated-subject setup.
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/link_hsps.c:331-374`
  - `s_RevCompareHSPsTbx` translated-query reverse-position comparator.
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/link_hsps.c:476-558`
  - Link wrapper allocation, reverse sort, frame-list construction.
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/link_hsps.c:589-982`
  - Even-gap dynamic programming, small/large gap selection, chain E-values,
    `linked_set`, and `ordering_method`.
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/link_hsps.c:990-1085`
  - Output replay order and chain traversal.
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_stat.c:4269-4567`
  - `s_BlastSumP`, `BLAST_SmallGapSumE`, `BLAST_LargeGapSumE`.

## Primary LOSAT Files

- `LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs`
  - Subject processing, linking call, E-value filter timing, output conversion.
- `LOSAT/src/algorithm/tblastx/blast_gapalign.rs`
  - HSP save order and `ScoreCompareHSPs` emulation.
- `LOSAT/src/algorithm/tblastx/reevaluate.rs`
  - Translated ungapped reevaluation.
- `LOSAT/src/algorithm/tblastx/filtering/purge_endpoints.rs`
  - Common-endpoint purge.
- `LOSAT/src/algorithm/tblastx/sum_stats_linking/cutoffs.rs`
  - `BlastLinkHSPParameters` cutoff calculation.
- `LOSAT/src/algorithm/tblastx/sum_stats_linking/linking.rs`
  - `BLAST_LinkHsps` / `s_BlastEvenGapLinkHSPs` port.
- `LOSAT/src/algorithm/tblastx/stage_dump.rs`
  - Diagnostic stage snapshots.
- `LOSAT/src/report/outfmt6.rs`
  - Final E-value and bit-score formatting.

## Phase 0: Freeze The Native-vs-NCBI Fixture

Objective: make the failing native-vs-NCBI comparison reproducible without
depending on the Wasm comparison script.

Tasks:

1. Add a native-focused comparison script or documented command that captures:
   - NCBI raw and sorted outfmt 6.
   - LOSAT native raw and sorted outfmt 6.
   - exact line counts and SHA-256.
   - full sorted `comm -3`.
   - HSP-key diff with E-value excluded.
   - score-bin summaries for NCBI-only and LOSAT-only HSP keys.
2. Record the NCBI executable path and `tblastx -version` in the scratch
   artifact.
3. Keep scratch artifacts under `.tmp/`; do not commit generated output.
4. Keep the existing Wasm parity script as a native/Wasm regression check, not
   as the main NCBI parity fixture.

Acceptance:

- One command reproduces `NCBI=4319`, `LOSAT=4350`.
- The report separates shared-HSP E-value deltas from true HSP-key deltas.
- The report can name the first score bin where HSP-key sets diverge.

## Phase 1: Align Diagnostic Stage Boundaries With NCBI

Objective: know exactly which NCBI-equivalent stage first creates the remaining
NCBI/LOSAT difference.

Tasks:

1. Extend LOSAT stage dumps, behind diagnostics only, around these boundaries:
   - after initial HSP save and `Blast_HSPListSortByScore` equivalent.
   - after common-endpoint purge if invoked for this path.
   - after `Blast_HSPListReevaluateUngapped`.
   - after `BLAST_LinkHsps`, before output coordinate conversion.
   - after preliminary E-value reaping.
   - after query coverage reaping.
   - after bit-score calculation.
   - after final output-order sorting.
2. For each stage, dump internal keys that match NCBI structures:
   - `context`
   - query offset/end
   - subject offset/end
   - subject frame
   - raw score
   - E-value double bits where available
   - chain flags and ordering method where available
3. If LOSAT stage dumps cannot identify the NCBI boundary, instrument a local
   NCBI checkout for diagnostics only. Do not copy NCBI implementation code
   into LOSAT.
4. Keep the NCBI call order exactly as in `blast_engine.c:1482-1541`. Do not
   move filters based on convenience.

Acceptance:

- The first divergent stage is identified as one of:
  - pre-reevaluation HSP content
  - reevaluation content
  - linking E-value assignment
  - preliminary E-value reaping
  - output ordering/formatting
- No code fix is attempted until this classification is recorded.

## Phase 2: Fix Shared-HSP E-value Parity First

Objective: reduce the large `3801` sorted diff by making E-values for shared
HSPs match NCBI before chasing low-score hit-set deltas.

Rationale: `4276` of `4319` NCBI HSP keys match LOSAT when only E-value text is
ignored, so an E-value calculation mismatch is the dominant defect.

Tasks:

1. Build a shared-HSP E-value comparator keyed by all non-E-value outfmt 6
   fields. For each shared key, report:
   - NCBI E-value text.
   - LOSAT E-value text.
   - numeric ratio where both parse as finite nonzero doubles.
   - raw score and bit score.
   - query frame and subject frame if recoverable from internal dump.
2. Classify E-value deltas:
   - uniform scale factor across many HSPs.
   - context-specific scale factor.
   - chain-length-specific deltas.
   - singletons vs linked sets.
   - small-gap vs large-gap ordering method.
3. Audit scalar single-HSP E-value parity against NCBI
   `Blast_HSPListGetEvalues`:
   - `kbp = sbp->kbp` for ungapped TBLASTX.
   - context index used for K/Lambda.
   - `query_info->contexts[hsp->context].eff_searchsp`.
   - subject length passed into the function.
   - no gapped `round_down` behavior for this path.
4. Audit linked-set E-value parity against NCBI `BLAST_LinkHsps`:
   - `xsum = prior_xsum + score * Lambda - logK`.
   - `query_length = MAX(context.query_length - length_adjustment, 1)`.
   - translated subject length and subject length adjustment division by three.
   - `BLAST_SmallGapSumE` and `BLAST_LargeGapSumE` formulas.
   - `BLAST_GapDecayDivisor` use for chain length.
   - `gap_prob` division for multi-HSP small/large gap chains.
5. Only after a source-backed mismatch is found, patch the smallest Rust area
   that produces the wrong E-value.

Acceptance:

- Shared-HSP E-value text matches NCBI for all HSP keys present in both outputs,
  or every remaining mismatch is explained by a later hit-set/order difference.
- Sorted diff falls from thousands of lines to approximately the true HSP-key
  difference.

## Phase 3: Fix Low-Score HSP-Set Parity

Objective: eliminate the remaining NCBI-only and LOSAT-only HSP-key deltas,
which are currently all below `23` bit.

Tasks:

1. Re-run the HSP-key classifier after Phase 2; do not use the old `43/74`
   numbers if E-value changes alter preliminary reaping.
2. For each remaining NCBI-only/LOSAT-only key, find the last LOSAT stage where
   its internal equivalent exists or disappears.
3. Compare disappearance timing against NCBI:
   - reevaluation deletion by cutoff.
   - common-endpoint purge deletion.
   - chain E-value assignment then preliminary E-value reaping.
   - query coverage reaping.
4. For reevaluation candidates, trace against
   `Blast_HSPReevaluateWithAmbiguitiesUngapped` exactly:
   - score reset when running sum drops below zero.
   - front-part discard when previous top score never reached cutoff.
   - update of best start/end offsets.
   - cutoff score by context.
5. For link/reap candidates, trace exact `hsp->evalue > prelim_evalue`
   decisions. The boundary must use the same double and comparison direction as
   NCBI.
6. For purge candidates, verify the exact common-endpoint sort keys and
   same-context/same-frame grouping.

Acceptance:

- NCBI-only HSP keys: `0`.
- LOSAT-only HSP keys: `0`.
- Hit count: `4319`.

## Phase 4: Restore Raw Output Order Parity

Objective: once the HSP set and E-values match, make raw output order match
NCBI without hiding differences by sorting.

Tasks:

1. Compare raw output after sorted output is byte-identical.
2. Audit final HSP-list and hit-list ordering against:
   - `ScoreCompareHSPs` in `blast_hits.c:1330-1354`.
   - TBLASTX translated replay in `link_hsps.c:990-1085`.
   - hit-list ordering by best E-value and subject OID.
3. Preserve NCBI comparator partial-order behavior. Use the established
   hidden list-order key only when needed to make deterministic emulation of
   comparator-equal rows match observed NCBI behavior.

Acceptance:

- Raw `diff` against NCBI is empty.
- Sorted `diff` remains empty.
- Native/Wasm raw outputs remain byte-identical.

## Phase 5: Generalize To Backlog Cases

Objective: ensure the fixture fix also addresses the known high-priority
TBLASTX parity backlog.

Cases:

- `LC738874.fasta` vs `LC738875.fasta`, `--outfmt 6 -n 1`.
- Long-sequence gencode 4 cases listed in `AGENTS.md`.
- Existing TBLASTX comparison scripts:
  - `bash tests/compare_self_tblastx.sh`
  - `bash tests/compare_long_sequences_debug.sh`
  - `bash compare_tblastx_results.sh`

Tasks:

1. Re-run native/Wasm determinism after every NCBI parity patch.
2. Re-run native-vs-NCBI fixture after every NCBI parity patch.
3. Re-run long-sequence diagnostics only after the small fixture reaches parity
   or a source-backed reason shows it exercises a different NCBI path.
4. If a long-sequence case diverges in a different stage, create a separate
   stage-classified issue instead of folding unrelated work into the fixture
   patch.

Acceptance:

- Small fixture has full raw parity.
- Long-sequence hit inflation is reduced or classified with a separate
  first-divergence stage.
- No NCBI executable is introduced into LOSAT runtime, build, fallback, or
  unsupported-feature paths.

## Implementation Rules For The Fix

- Do not patch based on numeric similarity alone. Every behavioral change must
  be tied to a specific NCBI source path and line range.
- Add NCBI source comments immediately above modified Rust code.
- Do not introduce new behavior absent from NCBI.
- Do not simplify sum-statistics formulas or substitute approximations.
- Preserve `f64`/C `double` operation order where it affects printed E-values.
- Run tests only after the current parity sweep target is patched, or when a
  specific diagnostic test is needed to confirm a source-backed hypothesis.
- Keep performance changes out of the parity PR.

## First Concrete Work Item

Create the native-vs-NCBI fixture classifier from Phase 0 and commit only the
script/documentation, not generated artifacts. The next engineering decision
should be based on that report:

```text
shared E-value mismatch first, or HSP-set divergence first
```

Given the current artifact, the expected next target is the shared-HSP E-value
path in `BLAST_LinkHsps` / sum statistics.

## Progress Log

### 2026-05-14 Phase 0 native-vs-NCBI fixture classifier

Added:

- `tests/compare_tblastx_native_ncbi_parity.sh`
  - One-command native LOSAT vs NCBI BLAST+ fixture capture.
  - Writes raw/sorted outputs, SHA-256 summaries, full sorted `comm -3`,
    NCBI executable path, and `tblastx -version`.
- `tests/classify_tblastx_outfmt6.py`
  - Classifies default outfmt 6 rows by exact full line, HSP key excluding only
    E-value, coordinate/identity key excluding E-value and bitscore, shared
    E-value text mismatches, shared bitscore text mismatches, and score bins.

Validation run:

```bash
bash tests/compare_tblastx_native_ncbi_parity.sh \
  .tmp/tblastx_native_ncbi_parity_codex_phase0
```

Oracle:

```text
/home/kawato/micromamba/bin/tblastx
tblastx: 2.17.0+
Package: blast 2.17.0, build Aug 11 2025 09:46:06
```

Observed current-worktree results:

```text
NCBI BLAST+ oracle:                         4319 hits
LOSAT native:                               4350 hits
sorted comm -3 lines:                       3739
exact full-line matches after sorting:      2465
matches ignoring only E-value:              4293
matches ignoring E-value and bitscore:      4319
shared HSP E-value mismatch rows:           1828
shared HSP bitscore mismatch rows:            26
NCBI-only keys ignoring E-value:              26
LOSAT-only keys ignoring E-value:             57
NCBI-only coordinate keys without scores:      0
LOSAT-only coordinate keys without scores:    31
LOSAT-only coordinate score bins:          22-29 bit floors
```

Interpretation:

- The current native output contains every NCBI coordinate/identity row.
- The remaining hit-count excess is 31 LOSAT-only low-score coordinate rows.
- The 26 NCBI-only rows from the E-value-only key comparison are not missing
  HSP coordinates; they are shared rows with bitscore text drift
  (`30.9` vs `30.8`, `61.6` vs `61.5`, `92.3` vs `92.2`).
- The dominant diff remains shared-HSP E-value text drift, so the next parity
  work should audit `BLAST_LinkHsps` / `Blast_HSPListGetEvalues` and bit-score
  calculation/formatting before chasing the 31 low-score LOSAT-only rows.

### 2026-05-14 Phase 2 first parity patch: context Karlin/ideal Kbp

Source-backed changes:

- `LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs`
  - Changed final outfmt bit-score calculation to use
    `contexts_ref[h.ctx_idx].karlin_params`, matching
    `blast_hits.c:1918-1928` (`kbp[hsp->context]`).
- `LOSAT/src/algorithm/tblastx/lookup/backbone.rs`
  - Changed translated-query `kbp_ideal` construction from static table lookup
    to NCBI-style standard-composition calculation, matching
    `blast_stat.c:2833-2848` (`Blast_ScoreBlkKbpIdealCalc`).

Validation run:

```bash
bash tests/compare_tblastx_native_ncbi_parity.sh \
  .tmp/tblastx_native_ncbi_parity_codex_ideal_kbp
```

Observed after patch:

```text
NCBI BLAST+ oracle:                         4319 hits
LOSAT native:                               4350 hits
sorted comm -3 lines:                       1159
exact full-line matches after sorting:      3755
matches ignoring only E-value:              4319
matches ignoring E-value and bitscore:      4319
shared HSP E-value mismatch rows:            564
shared HSP bitscore mismatch rows:             0
NCBI-only keys ignoring E-value:               0
LOSAT-only keys ignoring E-value:             31
NCBI-only coordinate keys without scores:      0
LOSAT-only coordinate keys without scores:    31
LOSAT-only coordinate score bins:          22-29 bit floors
```

Impact:

- Eliminated all shared-HSP bitscore text drift.
- Reduced shared-HSP E-value mismatch rows from `1828` to `564`.
- Reduced sorted diff from `3739` to `1159`.
- Hit-count discrepancy is unchanged and remains exactly the 31 low-score
  LOSAT-only coordinate rows.

Residual E-value classification:

```text
near 1% ratio window: 124 rows
far ratio rows:      440 rows
```

Next direction:

- Continue in `BLAST_LinkHsps` / sum-statistics chain formation and chain
  E-value assignment. The far-ratio rows show linked-set E-values assigned in
  different chains/orderings, not simple scalar `Blast_HSPListGetEvalues`
  drift.

### 2026-05-14 Phase 2 chain-ordering trace: small/large selection split

Source-backed diagnostic change:

- `LOSAT/src/algorithm/tblastx/sum_stats_linking/linking.rs`
  - Extended the existing `LOSAT_TRACE_HSP` path to print the two NCBI
    `link_hsps.c:901-952` ordering candidates:
    `prob[0] = BLAST_SmallGapSumE(...)`, `prob[1] = BLAST_LargeGapSumE(...)`,
    and the selected `prob[0] <= prob[1] ? small : large` method.
  - This is trace-only diagnostic output; normal runtime output and algorithms
    are unchanged.

Targeted traces:

```bash
LOSAT_TRACE_HSP=187401,187733,241784,242116 \
  LOSAT/target/release/LOSAT tblastx \
  -q LOSAT/tests/fasta/LC738874.fasta \
  -s LOSAT/tests/fasta/LC738875.fasta \
  --outfmt 6 -n 1 -o .tmp/trace_hsp_187401_241784.out

LOSAT_TRACE_HSP=187429,187680,241812,242063 \
  LOSAT/target/release/LOSAT tblastx \
  -q LOSAT/tests/fasta/LC738874.fasta \
  -s LOSAT/tests/fasta/LC738875.fasta \
  --outfmt 6 -n 1 -o .tmp/trace_hsp_187429_241812.out
```

Observed chain split in the representative positive/positive region:

- NCBI assigns E-value `3.85e-49` to the 5-HSP chain:
  `186667-186831`, `186839-186919`, `186987-187223`,
  `187215-187379`, `187401-187733`.
- LOSAT assigns `2.42e-52` to a different 5-HSP small-gap chain:
  `186667-186831`, `186839-186919`, `186987-187223`,
  `187225-187332`, `187429-187680`.
- NCBI assigns `1.64e-57` to the large-gap chain:
  `187225-187332`, `187429-187680`, `187847-187906`.
- LOSAT keeps `187847-187906` in a later/smaller chain with
  `187401-187733` (`1.74e-32`).

New LOSAT ordering trace for `187429-187680 / 241812-242063`:

```text
prob_small=2.420126e-52
prob_large=2.174866e-39
selected=0
best_small=head q_aa=62222..62277 s_aa=80348..80403 num0=5
best_large=head q_aa=2421..2451 s_aa=662..692 num1=60
```

Interpretation:

- The local discrepancy is not caused by raw-score or bit-score drift; NCBI
  `score` outfmt confirms the same raw scores for this region.
- LOSAT selects the local 5-HSP small-gap chain because the current global
  best large-gap candidate has a worse E-value (`2.17e-39`), even though it has
  a larger adjusted large-gap sum.
- NCBI's final output implies earlier large-gap ordering/removal differs before
  this local region is selected. The next audit target is therefore the
  `BLAST_LargeGapSumE` candidate ranking/removal path across the full
  positive/positive group, not the local coordinate conversion for these HSPs.

Follow-up group check:

- NCBI has 35 final rows at E-value `1.64e-57`, including
  `7264-7353 / 1987-2076`, `187225-187332 / 241593-241700`,
  `187429-187680 / 241812-242063`, and
  `187847-187906 / 242203-242262`.
- LOSAT places `7264-7353 / 1987-2076` in a 49-row `2.36e-35` group, while
  placing `187225-187332` and `187429-187680` in the small-gap `2.42e-52`
  group and `187847-187906` in `1.74e-32`.
- This points to a large-gap chain membership/removal-order divergence before
  the local 187k region is processed. The next concrete diagnostic should dump
  selected large-gap chain members for the `7264-7353 / 1987-2076` head and
  compare them against the NCBI `1.64e-57` final group.

Validation notes:

- `rustfmt` on the single edited file completed.
- `cargo build --release` completed successfully with existing warnings.
- Full parity comparison was not rerun after this trace-only change because no
  normal execution path changed.

### 2026-05-14 Phase 2 continuation: chain-selection removal-order trace

Source-backed diagnostic change:

- `LOSAT/src/algorithm/tblastx/sum_stats_linking/linking.rs`
  - Added `LOSAT_TRACE_CHAIN_HSP="qstart,qend,sstart,send"` to dump the selected
    chain that contains a target HSP without enabling the heavy per-candidate
    DP trace.
  - Added `LOSAT_TRACE_LINK_SELECTIONS=1` to dump every selected chain summary
    at the NCBI `link_hsps.c:901-982` selection/removal boundary.
  - These are diagnostic-only paths. With the environment variables unset,
    normal runtime output and linking behavior are unchanged.

Targeted trace:

```bash
LOSAT_TRACE_CHAIN_HSP=7264,7353,1987,2076 \
  LOSAT_TRACE_LINK_SELECTIONS=1 \
  LOSAT/target/release/LOSAT tblastx \
  -q LOSAT/tests/fasta/LC738874.fasta \
  -s LOSAT/tests/fasta/LC738875.fasta \
  --outfmt 6 -n 1 \
  -o .tmp/trace_chain_selection_7264.out \
  > .tmp/trace_chain_selection_7264.stdout \
  2> .tmp/trace_chain_selection_7264.log
```

Observed selection order:

- The traced HSP `7264-7353 / 1987-2076` is selected in
  `context=1 subject_sign=1` at LOSAT round `28`.
- At rounds `8` and `9`, LOSAT selects small-gap chains around
  `182012-182188 / 244446-244622` and
  `186681-186848 / 244455-244622`.
- At those earlier rounds, the best large-gap candidate is still headed by the
  traced `7264-7353 / 1987-2076` HSP, with a larger large-gap chain
  (`num1=61`, `sum1=1110`) than the chain finally selected at round `28`
  (`num1=49`, `sum1=929`).
- By the time the traced HSP is selected, earlier removals have reduced it to a
  49-member large-gap chain and LOSAT assigns E-value `2.357669e-35`.
- NCBI final output instead places this HSP in the `1.64e-57` group, while the
  high-scoring low-branch HSPs that LOSAT links into the 49-member target chain
  appear in NCBI's separate `1.47e-38` group.

Interpretation:

- The remaining divergence is not explained by raw-score, bit-score, or local
  HSP coordinate drift.
- The concrete discrepancy is now at the global chain selection/removal order:
  LOSAT removes one or more small-gap chains before NCBI's large-gap target
  chain is emitted, which changes the later large-gap membership and assigned
  E-values.

Next direction:

- Audit the Rust implementation of NCBI `link_hsps.c:560-982`, especially the
  candidate dominance tests, `xsum` updates, `best[0]/best[1]` selection, and
  `path_changed`/validation timing.
- Use the new trace to compare round `8`/`9` small-gap candidates against the
  large-gap target candidate and identify which NCBI ordering rule should have
  made the large-gap chain win earlier or survive unchanged.

Validation notes:

- `rustfmt --edition 2021 LOSAT/src/algorithm/tblastx/sum_stats_linking/linking.rs`
  completed.
- `cargo check --manifest-path LOSAT/Cargo.toml` completed successfully with
  existing warnings.
- `cargo build --release --manifest-path LOSAT/Cargo.toml` completed
  successfully with existing warnings.

### 2026-05-14 Phase 2 patch: preliminary link/reap before reevaluation

Source-backed behavior change:

- `LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs`
  - Added the missing NCBI preliminary `BLAST_LinkHsps` plus
    `s_Blast_HSPListReapByPrelimEvalue` pass after the initial ungapped HSP
    list is materialized and before `Blast_HSPListReevaluateUngapped`.
  - NCBI source anchor:
    `ncbi-blast/c++/src/algo/blast/core/blast_engine.c:870-899` links and
    prelim-evalue-reaps the raw `hsp_list_out` inside
    `s_BlastSearchEngineCore`.
  - The existing reevaluation and final relink path remains aligned to
    `blast_engine.c:1492-1535`, where NCBI reevaluates translated-subject HSPs,
    calls `BLAST_LinkHsps` again, then reapplies the preliminary e-value
    filter.
  - Added diagnostic-only `LOSAT_DEBUG_CUTOFFS_ALL=1` output for per-context
    `cutoff_score_max`, `gap_trigger`, and final `word_params` cutoff values,
    anchored to `blast_parameters.c:943-946`, `360-374`, and `401-416`.

Why this fixed the main discrepancy:

- NCBI `-evalue` was not just a final output threshold for TBLASTX. It also
  removes HSPs after the first sum-statistics linking pass, before ambiguity
  reevaluation and final relinking.
- LOSAT previously skipped this preliminary link/reap boundary, so high
  preliminary-evalue HSPs survived into the final linker. This made LOSAT
  default behavior match NCBI high-threshold behavior (`-evalue 10000`) for
  the traced `7264-7353 / 1987-2076` HSP.

Targeted validation:

```bash
RAYON_NUM_THREADS=1 LOSAT/target/release/LOSAT tblastx \
  -q LOSAT/tests/fasta/LC738874.fasta \
  -s LOSAT/tests/fasta/LC738875.fasta \
  -e 10 --outfmt 6 -n 1 -o .tmp/losat_prelink_e10.out

RAYON_NUM_THREADS=1 LOSAT/target/release/LOSAT tblastx \
  -q LOSAT/tests/fasta/LC738874.fasta \
  -s LOSAT/tests/fasta/LC738875.fasta \
  -e 100 --outfmt 6 -n 1 -o .tmp/losat_prelink_e100.out

RAYON_NUM_THREADS=1 LOSAT/target/release/LOSAT tblastx \
  -q LOSAT/tests/fasta/LC738874.fasta \
  -s LOSAT/tests/fasta/LC738875.fasta \
  -e 10000 --outfmt 6 -n 1 -o .tmp/losat_prelink_e10000.out
```

Results against the matching NCBI oracle files:

- `-evalue 10`: LOSAT row count `4319`, NCBI row count `4319`.
  The traced HSP now has E-value `1.64e-57`, matching NCBI default/e10.
  Full key match on `qstart/qend/sstart/send/evalue/bitscore` is `4317/4319`;
  the two remaining rows are an E-value swap between
  `197216-197248 / 8251-8219` and `197216-197248 / 8249-8217`.
- `-evalue 100`: LOSAT row count `6959`, NCBI row count `6959`.
  The traced HSP now has E-value `8.13e-36`, matching NCBI.
  Full key match is `6957/6959`.
- `-evalue 10000`: LOSAT row count `49353`, NCBI row count `49353`.
  The traced HSP remains at E-value `2.36e-35`, matching NCBI high-threshold
  behavior. Full key match is `49349/49353`.

Validation notes:

- `rustfmt --edition 2021 LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs`
  completed.
- `cargo build --release --manifest-path LOSAT/Cargo.toml` completed
  successfully with existing warnings.

Next direction:

- Investigate the remaining E-value swaps among near-duplicate short HSPs.
  The first concrete e10 pair is:
  - NCBI: `197216 197248 8251 8219 frame 2/-1 score 47 e=0.93`
  - NCBI: `197216 197248 8249 8217 frame 2/-3 score 47 e=2.3`
  - LOSAT has the same coordinates and bit scores, but the two E-values are
    assigned to the opposite output rows.
- Focus on NCBI `link_hsps.c` frame grouping/sort tie behavior and subject
  translation frame ordering for duplicate short HSPs, not on global hit count
  or cutoff policy.

### 2026-05-14 Phase 2 patch: per-frame append score sort and qsort parity

Source-backed behavior change:

- `LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs`
  - After each subject-frame ungapped HSP list is converted and appended to the
    combined translated-subject HSP list, LOSAT now calls the same NCBI-style
    score sort used by `Blast_HSPListAppend`.
  - NCBI source anchors:
    - `blast_engine.c:804-845` loops translated subject frames and calls
      `Blast_HSPListAppend(&hsp_list_for_chunks, &hsp_list_out, kHspNumMax)`.
    - `blast_hits.c:2758-2766` appends the new HSP pointers to the existing
      combined list and immediately calls `Blast_HSPListSortByScore`.
    - `blast_hits.c:1374-1381` performs the guarded `qsort(...,
      ScoreCompareHSPs)` used by that sort.
  - This is the ordering boundary that controls comparator-equal duplicate HSPs
    before the first `BLAST_LinkHsps` call.
- `LOSAT/src/algorithm/tblastx/ncbi_qsort.rs`
  - Added a native qsort wrapper for `UngappedHit` pointer-array parity. It
    sorts an index array with C `qsort`, then replays that order onto Rust HSP
    values, matching NCBI's `BlastHSP*` / `LinkHSPStruct*` qsort calls without
    moving Rust structs through C.
  - Source anchors: `blast_hits.c:1379-1381` for `ScoreCompareHSPs` qsort and
    `link_hsps.c:483-486` for `s_RevCompareHSPsTbx` qsort.
- `LOSAT/src/algorithm/tblastx/blast_gapalign.rs` and
  `LOSAT/src/algorithm/tblastx/sum_stats_linking/linking.rs`
  now route NCBI score/link sort calls through the qsort wrapper on native
  builds. Wasm keeps the Rust comparator fallback because NCBI BLAST+ native
  parity is the qsort tie-behavior target.

Why this fixed the residual E-value swap:

- The remaining e10 mismatch was not a cutoff or count issue. It was a chain
  member assignment swap between two duplicate short HSPs with the same
  frame-relative AA coordinates and raw score but different subject frames.
- NCBI sorts the combined subject-frame HSP list after every frame append before
  linking. LOSAT previously extended the combined list and deferred the sort,
  so the `link_hsps.c` DP saw duplicate candidates in the wrong order.

Targeted validation:

```bash
cargo build --release --manifest-path LOSAT/Cargo.toml

RAYON_NUM_THREADS=1 LOSAT/target/release/LOSAT tblastx \
  -q LOSAT/tests/fasta/LC738874.fasta \
  -s LOSAT/tests/fasta/LC738875.fasta \
  -e 10 --outfmt 6 -n 1 -o .tmp/losat_framesort_e10.out

RAYON_NUM_THREADS=1 LOSAT/target/release/LOSAT tblastx \
  -q LOSAT/tests/fasta/LC738874.fasta \
  -s LOSAT/tests/fasta/LC738875.fasta \
  -e 100 --outfmt 6 -n 1 -o .tmp/losat_framesort_e100.out

RAYON_NUM_THREADS=1 LOSAT/target/release/LOSAT tblastx \
  -q LOSAT/tests/fasta/LC738874.fasta \
  -s LOSAT/tests/fasta/LC738875.fasta \
  -e 10000 --outfmt 6 -n 1 -o .tmp/losat_framesort_e10000.out
```

Results against the matching NCBI oracle files:

- `-evalue 10`: row counts `4319/4319`; full-key match on
  `qstart/qend/sstart/send/evalue/bitscore` is `4319/4319`.
  The former swap now matches NCBI:
  - `197216-197248 / 8251-8219`: `e=0.93`
  - `197216-197248 / 8249-8217`: `e=2.3`
- `-evalue 100`: row counts `6959/6959`; full-key match `6959/6959`.
- `-evalue 10000`: row counts `49353/49353`; full-key match `49353/49353`.

Validation notes:

- `rustfmt --edition 2021` completed for the touched TBLASTX files.
- `cargo build --release --manifest-path LOSAT/Cargo.toml` completed
  successfully with existing warnings.

Next direction:

- Re-run the larger long-sequence TBLASTX comparison set before changing
  linking heuristics again. The current short LC738874/LC738875 threshold sweep
  is bit-perfect for counts, coordinates, E-values, and bit scores at
  `-evalue 10`, `100`, and `10000`.

### 2026-05-14 Phase 5 long-sequence rerun: AP027131 vs AP027133 gencode 4

Validation baseline:

```bash
bash tests/compare_tblastx_native_ncbi_parity.sh \
  .tmp/tblastx_native_ncbi_parity_codex_current
```

Result:

- `LC738874.fasta` vs `LC738875.fasta`, default `-evalue 10`,
  `--outfmt 6`, `-n 1` is now fully byte-identical after sorting:
  `4319/4319`, `ncbi_vs_native_sorted_comm3_lines=0`.

Long-sequence rerun:

```bash
RAYON_NUM_THREADS=1 LOSAT/target/release/LOSAT tblastx \
  -q LOSAT/tests/fasta/AP027131.fasta \
  -s LOSAT/tests/fasta/AP027133.fasta \
  --query-gencode 4 --db-gencode 4 \
  --outfmt 6 -n 1 \
  -o .tmp/tblastx_long_ap027131_ap027133_current/native.raw.out

tblastx \
  -query LOSAT/tests/fasta/AP027131.fasta \
  -subject LOSAT/tests/fasta/AP027133.fasta \
  -query_gencode 4 -db_gencode 4 \
  -outfmt 6 -num_threads 1 \
  -out .tmp/tblastx_long_ap027131_ap027133_current/ncbi.raw.out
```

Observed classification:

```text
NCBI BLAST+ oracle:                         14871 hits
LOSAT native:                               29926 hits
ratio:                                      2.01x
sorted comm -3 lines:                       37015
exact full-line matches after sorting:       3891
matches ignoring only e-value:               8667
matches ignoring e-value and bitscore:       11248
shared HSP E-value mismatch rows:            4776
shared HSP bitscore mismatch rows:           2581
NCBI-only HSP keys ignoring e-value:         6204
LOSAT-only HSP keys ignoring e-value:       21259
NCBI-only alignment keys without scores:     3623
LOSAT-only alignment keys without scores:   18678
alignment-key diff bit-score floor range:   22-591
```

Representative high-score trace:

- Coordinate shared between outputs:
  `491844-489349 / 366336-363841`, frame `-3/-3`.
- NCBI custom outfmt with raw score:
  `bitscore=1533`, `score=3340`.
- LOSAT trace:
  `init_hsp_saved raw_score=3520`,
  `after_reevaluate raw_score=3520`,
  `output_hit raw_score=3520 bit=1615.8`.
- A direct unmasked BLOSUM62 calculation over the LOSAT translated frame
  interval `q_aa=56754..57586`, `s_aa=79952..80784` gives raw score `3581`
  with `665/832` identities.

Interpretation:

- The long-sequence 2x case is still open; the short-fixture linking/reap/qsort
  parity patch did not generalize to this gencode 4 case.
- At least one top shared-coordinate HSP has a raw-score mismatch before output
  formatting, so the immediate offender is not final bit-score calculation.
- LOSAT is already scoring below the unmasked translated-sequence score, but
  NCBI scores substantially lower. This points to SEG mask coverage or the
  masked sequence used by ungapped extension / `Blast_HSPListReevaluateUngapped`,
  not global link ordering as the first next target.
- Do not change sum-statistics linking again until SEG/reevaluation input
  parity is audited against NCBI source and a source-backed mismatch is found.

Next direction:

- Audit NCBI query SEG setup and translated ungapped reevaluation for this path:
  `blast_filter.c` query masking setup, `blast_seg.c`/SEG interval semantics,
  `blast_hits.c:675-733` and `2609-2737`
  (`Blast_HSPReevaluateWithAmbiguitiesUngapped` /
  `Blast_HSPListReevaluateUngapped`).
- Add a targeted diagnostic, behind an environment variable, to dump SEG mask
  intervals and masked-residue counts for a traced TBLASTX HSP. The first target
  should be `LOSAT_TRACE_HSP=491844,489349,366336,363841` on
  `AP027131.fasta` vs `AP027133.fasta`, gencode 4.
- Use NCBI only as a comparison oracle/source reference. Do not add any NCBI
  runtime dependency.

### 2026-05-14 Phase 5 diagnostic update: SEG ruled out for the top mismatch

Added a targeted diagnostic gate, `LOSAT_TRACE_HSP_MASKS=1`, for traced TBLASTX
HSPs. For the representative AP027131/AP027133 mismatch
`491844-489349 / 366336-363841`, frame `-3/-3`, LOSAT reports:

```text
raw_score=3520 masked_score=3520 unmasked_score=3581
masked_residues=12 mask_runs=209..221
```

NCBI oracle was re-run on a restricted region to avoid the full long-sequence
`-seg no` runtime:

```bash
tblastx \
  -query LOSAT/tests/fasta/AP027131.fasta \
  -subject LOSAT/tests/fasta/AP027133.fasta \
  -query_loc 480000-500000 \
  -subject_loc 350000-380000 \
  -query_gencode 4 -db_gencode 4 \
  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore score qframe sframe' \
  -num_threads 1 \
  -out .tmp/tblastx_long_ap027131_ap027133_current/ncbi.region.score_frames.seg_default.out

tblastx \
  -query LOSAT/tests/fasta/AP027131.fasta \
  -subject LOSAT/tests/fasta/AP027133.fasta \
  -query_loc 480000-500000 \
  -subject_loc 350000-380000 \
  -query_gencode 4 -db_gencode 4 \
  -seg no \
  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore score qframe sframe' \
  -num_threads 1 \
  -out .tmp/tblastx_long_ap027131_ap027133_current/ncbi.region.score_frames.seg_no.out
```

NCBI scores for the same HSP:

- default SEG: `score=3340`, `bitscore=1533`
- `-seg no`: `score=3401`, `bitscore=1561`

The NCBI SEG delta is `61`, exactly matching the LOSAT masked/unmasked delta
`3581 - 3520 = 61`. Therefore the remaining `180` raw-score overrun
(`3520 - 3340`, or `3581 - 3401` without SEG) is not a SEG interval or
stop-codon matrix problem. A small synthetic stop-codon oracle confirmed NCBI
scores a displayed `*/*` aligned residue pair as `+1`, matching the static
BLOSUM62 table.

Current interpretation:

- Query masking for this top shared-coordinate HSP is now likely in parity.
- `Blast_HSPReevaluateWithAmbiguitiesUngapped` alone does not explain the
  mismatch, because the printed NCBI `qseq/sseq` full BLOSUM62 sum is `3581`
  while the NCBI `-seg no` raw score remains `3401`.
- The next source-backed target is the initial ungapped extension score path,
  especially `aa_ungapped.c:s_BlastAaExtendTwoHit` /
  `s_BlastAaExtendLeft` / `s_BlastAaExtendRight`, plus the timing of when that
  initial score is preserved or replaced before linking.

Validation after adding the diagnostic:

- `cargo check --manifest-path LOSAT/Cargo.toml` passed with existing warnings.
- `cargo build --release --manifest-path LOSAT/Cargo.toml` passed with existing
  warnings.

Next direction:

- Instrument or audit the LOSAT two-hit extension result for the traced HSP
  against the NCBI `aa_ungapped.c` loop structure and score-preservation timing.
- Avoid further linking or sum-statistics changes until the `3520` vs `3340`
  initial-score discrepancy is explained by a concrete NCBI source reference.

### 2026-05-14 Phase 5 resolution: local-subject DB genetic-code split

Root cause:

- The remaining `180` raw-score overrun was not caused by two-hit extension
  loop bounds. It came from local `-subject` TBLASTX genetic-code handling for
  translated subjects.
- NCBI BLAST+ with local `-subject` and `-query_gencode 4 -db_gencode 4`
  displays TGA as `W` in subject `qseq/sseq`, but the raw search score is
  computed as though the local subject search translation still used the
  standard code for TGA. A synthetic `ATGTGA...` oracle showed the key behavior:
  displayed `MWMW.../MWMW...` has raw score `11`, exactly the score for
  `M/W` positions against standard-code subject stops (`W/* = -4`), not the
  BLOSUM62 sum for displayed `W/W`.
- Source boundary used for the port:
  - `blast_engine.c:772-775`: search-time translated-subject frames are built
    from `subject->gen_code_string`.
  - `blast_engine.c:1458-1465`: the `db_options->genetic_code` fallback is on
    the database `seq_src` path.
  - `blast_hits.c:2708-2717`: identities are recomputed after reevaluation,
    using the reporting sequence buffers.

Implementation:

- LOSAT now keeps two translated subject frame sets in the native `-s/--subject`
  path:
  - search/scoring frames use standard genetic code 1, matching NCBI local
    `-subject` raw-score behavior;
  - report/identity frames use `--db-gencode`, preserving displayed coordinate
    and identity behavior for code 4.
- `LOSAT_TRACE_HSP_MASKS=1` now prints both the score against the search subject
  translation and the displayed/report translation sum.

Representative trace after the fix:

```text
[TRACE_HSP] seed/extend ... score=3340 ...
[TRACE_HSP] stage=after_reevaluate raw_score=3340 ...
[TRACE_HSP_MASKS] ... raw_score=3340 masked_score=3340 unmasked_score=3401 report_unmasked_score=3581 masked_residues=12 mask_runs=209..221
[TRACE_HSP] stage=output_hit raw_score=3340 bit=1533.3 ...
```

Validation:

- Synthetic TGA code-4 oracle:
  - NCBI: displayed `MWMWMWMWMWMWM/MWMWMWMWMWMWM`, `score=11`,
    `bitscore=7.9`.
  - LOSAT trace now reports `raw_score=11`, `bit=7.9` for the same
    `q=1-39 / s=1-39`, frame `1/1`.
- Long AP027131/AP027133 gencode-4 rerun:

```text
ncbi_line_count                                  14871
native_line_count                                14871
exact_full_line_matches_after_sorting            14871
matches_ignoring_only_evalue                     14871
matches_ignoring_evalue_and_bitscore             14871
shared_hsp_evalue_mismatch_rows                  0
shared_hsp_bitscore_mismatch_rows                0
ncbi_only_hsp_keys_ignoring_evalue               0
native_only_hsp_keys_ignoring_evalue             0
ncbi_only_alignment_keys_ignoring_evalue+score   0
native_only_alignment_keys_ignoring_evalue+score 0
```

- Short LC738874/LC738875 default parity fixture still matches exactly:
  `4319/4319`, with zero E-value, bit-score, HSP-key, or alignment-key
  mismatches.
- `cargo check --manifest-path LOSAT/Cargo.toml` passed with existing warnings.
- `cargo build --release --manifest-path LOSAT/Cargo.toml` passed with existing
  warnings.

Next direction:

- The long-sequence AP027131/AP027133 gencode-4 2x-hit issue is resolved for
  local `-subject` mode. Keep this case as a regression fixture.
- Do not remove the search/report subject translation split unless a new NCBI
  source-backed oracle demonstrates different behavior for a non-local database
  path.
