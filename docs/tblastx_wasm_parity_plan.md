# TBLASTX Wasm Parity Implementation Plan

Date: 2026-05-13

Goal: make LOSAT TBLASTX Wasm output bit-identical to the authoritative NCBI
BLAST behavior, and therefore identical to the validated native LOSAT path for
the same inputs and options. Performance work must wait until Wasm parity is
restored.

NCBI BLAST remains the only source of truth. NCBI binaries may be used only as
comparison oracles during diagnostics and validation. LOSAT runtime code must
not call NCBI BLAST, link to NCBI BLAST, or use NCBI BLAST as a fallback.

## Current Observation

The current release builds show matching hit counts but byte-different TBLASTX
outfmt 6 output:

```text
native: target/release/LOSAT
wasm:   target/wasm32-wasip1/release/LOSAT.wasm, run through Node WASI
case:   tblastx -q tests/fasta/LC738874.fasta -s tests/fasta/LC738875.fasta --outfmt 6 -n 1
hits:   native 4350, wasm 4350
```

Sorted output hashes differ, so this is not only final line ordering. Examples
from the observed diff:

- Same coordinates with E-values exchanged between neighboring HSPs:
  `267767..267735` and `267769..267737` swap `2.3` and `2.9`.
- Similar low-score HSPs shift by one codon-like step on subject coordinates:
  native `250987..251037`, Wasm `250989..251039`.
- Differences cluster around low bit scores (`22.1` to `24.4`) and E-values near
  the reporting threshold.

Primary suspicion order:

1. Wasm SIMD TBLASTX extension or reevaluation deviates from scalar/NCBI in
   endpoint selection.
2. Sum-statistics linking assigns chain E-values to different neighboring HSPs
   because an earlier sort/link order diverged.
3. Final output sort reconstructs internal coordinates from formatted output
   coordinates and loses a TBLASTX-specific tie-break, especially on negative
   subject frames.

## Primary NCBI References

Read these before changing Rust code. Every parity code change must include the
relevant NCBI snippet comment with file path and line numbers.

- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/aa_ungapped.c:846-921`
  - `s_BlastAaExtendRight`, `s_BlastAaExtendLeft`, X-drop termination, best
    endpoint selection.
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:675-733`
  - `Blast_HSPReevaluateWithAmbiguitiesUngapped`, score reset and best segment
    selection.
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1330-1377`
  - `ScoreCompareHSPs`, HSP sort order by raw score and internal coordinates.
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:2268-2535`
  - common endpoint purge sort keys and duplicate handling.
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/link_hsps.c:483-558`
  - TBLASTX reverse-position sort and frame-group construction.
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/link_hsps.c:589-982`
  - chain dynamic programming, small/large gap E-value selection, and
    `linked_set`/`start_of_chain` assignment.
- `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/link_hsps.c:990-1059`
  - translated-query output replay order and chain-member traversal.

Primary Rust files:

- `LOSAT/src/algorithm/tblastx/extension/ungapped.rs`
- `LOSAT/src/algorithm/tblastx/extension/two_hit.rs`
- `LOSAT/src/algorithm/tblastx/reevaluate.rs`
- `LOSAT/src/algorithm/tblastx/blast_gapalign.rs`
- `LOSAT/src/algorithm/tblastx/filtering/purge_endpoints.rs`
- `LOSAT/src/algorithm/tblastx/sum_stats_linking/linking.rs`
- `LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs`
- `LOSAT/src/common.rs`

## Phase 0: Reproduce and Freeze the Failing Case

Objective: create a small, repeatable parity fixture before modifying logic.

Tasks:

1. Add a documented comparison command for the current failing case:
   `LC738874.fasta` vs `LC738875.fasta`, TBLASTX, `--outfmt 6`, `-n 1`.
2. Capture four outputs for the same command:
   - NCBI BLAST+ TBLASTX oracle.
   - Native LOSAT release.
   - Wasm LOSAT release with SIMD enabled.
   - Wasm LOSAT release with Wasm SIMD disabled.
3. Store comparison artifacts under a scratch directory, not committed output
   files:
   - raw outfmt 6
   - sorted outfmt 6
   - `comm -3` of sorted outputs
   - line count
   - SHA-256 of raw and sorted outputs
4. Treat hit-count equality as diagnostic only. Acceptance requires byte parity
   against NCBI output.

Acceptance:

- A single command reproduces the native/Wasm mismatch.
- The diff is classified as one of:
  - Wasm SIMD only
  - all Wasm builds
  - native/Wasm both differ from NCBI
  - final formatting/order only

## Phase 1: Isolate SIMD vs Scalar Behavior

Objective: determine whether the first divergence is caused by Wasm SIMD.

Current `.cargo/config.toml` enables `+simd128` for `wasm32-wasip1`, and the
TBLASTX code has Wasm SIMD branches in extension and reevaluation.

Tasks:

1. Build a scalar Wasm diagnostic binary by temporarily disabling `simd128` or
   by adding a parity-only feature gate around these branches:
   - `reevaluate_ungapped_hit_ncbi_translated_wasm_simd128`
   - `extend_left_ungapped_simd128`
   - `extend_right_ungapped_simd128`
   - both one-hit and two-hit call sites
2. Compare:
   - native release vs native scalar path if available
   - native release vs Wasm scalar
   - Wasm scalar vs Wasm SIMD
3. If scalar Wasm matches native/NCBI but SIMD does not, fix only the offending
   SIMD routine and keep the scalar routine as the reference implementation.
4. If scalar Wasm still diverges, move to stage snapshots before changing SIMD.

Acceptance:

- The plan identifies whether Wasm SIMD is necessary and sufficient to
  reproduce the parity failure.
- No performance optimization is made until scalar/SIMD parity is understood.

## Phase 2: Add Stage Snapshot Diagnostics

Objective: find the first pipeline stage where native and Wasm differ.

Add diagnostics behind an explicit environment variable, for example
`LOSAT_DUMP_TBLASTX_STAGE=/tmp/losat-stage`, with no default runtime output
change. Dumps must use internal values, not formatted output coordinates.

Snapshot rows should include:

- stage name
- `q_idx`, `s_idx`, `ctx_idx`, `s_f_idx`
- query frame and subject frame
- internal amino-acid starts/ends: `q_aa_start`, `q_aa_end`, `s_aa_start`,
  `s_aa_end`
- seed offsets: `q_seed_off`, `s_seed_off`
- `raw_score`
- `e_value` and `e_value.to_bits()` in hex
- `ordering_method`
- `linked_set`, `start_of_chain`, `link_id`, `chain_next_link_id`

Dump stages:

1. After initial ungapped extension, before common-endpoint purge.
2. After common-endpoint purge.
3. After `Blast_HSPListReevaluateUngapped` equivalent.
4. Before `link_hsps`.
5. After `link_hsps`, before output conversion.
6. After output conversion, before final output sort.
7. After final output sort.

Acceptance:

- Native and Wasm snapshots can be compared with `sort`/`comm` using stable keys.
- The first divergent stage is known.
- If differences are only floating-point, the raw `f64::to_bits()` divergence is
  visible before decimal formatting.

## Phase 3: Fix Extension/Reevaluation If It Is the First Divergence

Objective: make Wasm endpoint and score selection match NCBI exactly.

Tasks:

1. For each Wasm SIMD routine, add unit-style parity tests that compare SIMD to
   scalar over boundary-heavy synthetic sequences:
   - best segment starts at lane 0 and lane 15
   - X-drop triggers inside a SIMD block
   - X-drop triggers in the scalar tail
   - score drops below zero inside a word
   - equal max-score ties, where NCBI keeps the first best endpoint because the
     update is strictly `>`, not `>=`
   - stop/sentinel-like residues and masked query residues
2. For `reevaluate_ungapped_hit_ncbi_translated_wasm_simd128`, verify that:
   - `current_start = pos + 1` after negative reset
   - `best_start` and `best_end` are updated only when `sum > score`
   - score below cutoff returns `None`
   - `new_s_off = s_off_raw + best_start`, not any subject-derived local offset
3. For extension routines, verify that returned scanned positions match scalar:
   - right extension returns the current zero-based position when termination
     fires
   - left extension applies lanes near-to-far
4. Prefer deleting a SIMD shortcut over preserving a fast but non-parity path.

Acceptance:

- SIMD and scalar return identical tuples for all added parity cases.
- The failing real TBLASTX stage no longer diverges at extension/reevaluation.

## Phase 4: Fix Linking/E-Value Assignment If It Is the First Divergence

Objective: ensure `link_hsps.c` ordering and chain E-value replay match NCBI.

Tasks:

1. Verify `rev_compare_hsps_tbx`, `rev_compare_hsps_transl`, and
   `fwd_compare_hsps_transl` against NCBI comparators, including complete
   tie-breakers. Avoid `sort_unstable_by` unless the comparator is total for all
   possible HSP pairs.
2. Snapshot frame-group boundaries from native and Wasm. The group count,
   group order, and HSP membership must match before link DP runs.
3. Snapshot selected `best_i`, `ordering`, `prob[0]`, `prob[1]`, chain length,
   and chain member link IDs for the first mismatching HSP cluster.
4. Compare floating-point inputs to sum statistics:
   - effective query length
   - effective subject length
   - effective search space
   - `logK`
   - gap decay divisor
5. If E-values differ only because equal-score neighboring HSPs are replayed in
   different order, fix the ordering comparator, not the E-value math.

Acceptance:

- Native and Wasm produce the same `linked_set`, `start_of_chain`,
  `ordering_method`, chain links, and `e_value.to_bits()` after linking.

## Phase 5: Fix Final Output Sorting If Internal HSPs Match

Objective: avoid final output differences caused by losing internal coordinates
during conversion to `Hit`.

Risk in current code:

- `common.rs::score_compare_hsps` reconstructs subject/query sort keys from
  output coordinates.
- NCBI `ScoreCompareHSPs` sorts by internal `BlastHSP` offsets/ends before
  output coordinate adjustment.
- For TBLASTX negative frames, `min/max` output coordinates may not be a
  lossless substitute for internal offsets.

Tasks:

1. If snapshots show internal HSPs match but final order differs, preserve
   nonprinted internal sort keys in `Hit` or perform the NCBI HSP sort before
   converting coordinates.
2. Include internal keys needed by `ScoreCompareHSPs`:
   - internal query offset/end
   - internal subject offset/end
   - context
   - subject frame
3. Keep outfmt 6/7 printed columns unchanged.
4. Add tests covering positive and negative query/subject frame combinations.

Acceptance:

- Final outfmt 6 order matches NCBI without changing printed hit content.
- Sorting does not depend on `HashMap` or `HashSet` iteration order where NCBI
  uses query/subject input order.

## Phase 6: Validation Matrix

Objective: prove the fix is not narrow to one fixture.

Run at least:

- `LC738874.fasta` vs `LC738875.fasta`
- `small_test.fasta` vs `small_test.fasta`
- `NZ_CP006932.fasta` vs `NZ_CP006932.fasta`
- one gencode 4 long-sequence case from the existing TBLASTX backlog
- one no-SEG diagnostic run if SEG masking affects the first divergence

For each:

- native release vs Wasm release
- Wasm scalar diagnostic vs Wasm SIMD release, until SIMD parity is trusted
- LOSAT vs NCBI BLAST+ oracle
- `--outfmt 6` and, if practical, `--outfmt 7`

Acceptance:

- Byte-identical LOSAT native/Wasm output for all validation cases.
- Byte-identical LOSAT/NCBI output where the existing native parity baseline is
  expected to pass.
- No new runtime dependency on NCBI BLAST.

## Phase 7: Only Then Resume Performance Work

After parity is fixed:

1. Re-run the native/Wasm timing table.
2. Profile only the routines confirmed parity-safe.
3. Keep scalar fallbacks as the reference in tests.
4. Optimize Wasm SIMD only behind tests that compare each SIMD helper to scalar
   and the full TBLASTX output to NCBI.

## Proposed First PR Scope

Keep the first implementation PR narrow:

1. Add the failing TBLASTX Wasm parity fixture command and comparison script.
2. Add scalar-vs-SIMD diagnostic build path or feature gate.
3. Add stage snapshot diagnostics behind `LOSAT_DUMP_TBLASTX_STAGE`.
4. Use those diagnostics to identify the first divergent stage.

Do not combine this with a performance optimization PR. The first PR should
either include the minimal parity fix if the culprit is obvious, or stop after
diagnostics with a documented first-divergence result.

## 2026-05-13 Investigation Update: First Divergence

Diagnostic run:

- Comparison artifacts: `.tmp/tblastx_wasm_parity_20260513_213603`
- Case:
  `tblastx -q LOSAT/tests/fasta/LC738874.fasta -s LOSAT/tests/fasta/LC738875.fasta --outfmt 6 -n 1`
- Output line counts:
  - native LOSAT: 4350
  - Wasm SIMD LOSAT: 4350
  - Wasm scalar LOSAT: 4350
  - NCBI BLAST+ oracle: 4319
- Sorted diff line counts:
  - native vs Wasm SIMD: 36
  - native vs Wasm scalar: 36
  - Wasm scalar vs Wasm SIMD: 0

Stage snapshots were captured with `LOSAT_DUMP_TBLASTX_STAGE` under
`/tmp/losat_tblastx_stage_20260513`. The first divergent stage is not the
initial extension:

- `after_initial_ungapped_extension.tsv` is byte-identical between native and
  Wasm.
- `after_common_endpoint_purge.tsv` is byte-identical between native and Wasm.
- `after_reevaluate.tsv` has different raw hashes, but sorted rows are
  identical (`comm -3` is empty).
- `before_link_hsps.tsv` also has identical sorted rows.

Diagnosis:

- This is not a Wasm SIMD-specific defect. Wasm scalar and Wasm SIMD are
  byte-identical for this fixture.
- Re-evaluation produces the same HSP content on native and Wasm for this
  fixture. The first divergence is row order immediately after the
  `Blast_HSPListReevaluateUngapped` equivalent.
- The current Rust `score_compare_ungapped_hits_ncbi` mirrors NCBI
  `ScoreCompareHSPs` from
  `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1330-1377`.
  That comparator is not total for distinct TBLASTX HSPs because it compares
  score and internal query/subject offsets/ends, but not context, query frame,
  subject frame, original list position, or seed identity.
- `ncbi_qsort_ungapped_hits_by_score` currently uses Rust
  `sort_unstable_by`. Comparator-equal HSPs can therefore be permuted
  differently by native and Wasm builds even when all HSP fields are otherwise
  identical.
- That order then feeds NCBI-style linking. NCBI `BLAST_LinkHsps` builds the
  temporary link array from the current HSP list order before sorting with
  `s_RevCompareHSPsTbx`
  (`/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/link_hsps.c:483-558`).
  The translated replay sorts in
  `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/link_hsps.c:990-1059`
  also use partial comparators. The initial order difference is therefore able
  to change chain links and assign neighboring low-score HSPs different
  E-values.

## Repair Plan: Comparator-Equal HSP Ordering

Objective: remove native/Wasm ordering nondeterminism while preserving NCBI
semantics and the existing NCBI comparator keys.

1. Add a hidden, nonprinted HSP list-order key to the TBLASTX internal HSP
   representation. Assign it at the NCBI-equivalent insertion/list creation
   boundary, before any `Blast_HSPListSortByScore`-equivalent step.
2. Replace `sort_unstable_by` calls that emulate NCBI HSP `qsort` behavior with
   a dedicated helper. The helper must apply the NCBI comparator first and use
   the captured list-order key only when the NCBI comparator returns equality.
   Do not add biological coordinate, frame, hash, or pointer-derived tie-breaks
   that are absent from the NCBI comparator.
3. Preserve the `Blast_HSPListIsSortedByScore` precheck semantics from
   `/mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1369-1377`.
   If the list is already sorted by the NCBI comparator, keep it in its current
   order.
4. Apply the ordering helper consistently to the first observed offender and
   then to later linking sorts only where stage snapshots prove comparator-equal
   HSPs can affect output:
   - `ncbi_qsort_ungapped_hits_by_score`
   - `hits.sort_unstable_by(rev_compare_hsps_tbx)`
   - `results.sort_unstable_by(rev_compare_hsps_transl)`
   - `results.sort_unstable_by(fwd_compare_hsps_transl)`
5. Add focused tests with comparator-equal TBLASTX HSPs that differ only in
   fields omitted by the NCBI comparators, such as context or frame. The tests
   should assert deterministic native/Wasm internal order and unchanged printed
   fields.
6. Re-run `bash tests/compare_tblastx_wasm_parity.sh`. For this fixture, the
   expected acceptance after this fix is:
   - native vs Wasm SIMD sorted diff lines: 0
   - native vs Wasm scalar sorted diff lines: 0
   - Wasm scalar vs Wasm SIMD sorted diff lines: 0
7. Keep NCBI oracle comparison separate in the report. This fixture still shows
   a native LOSAT vs NCBI hit-count gap before this repair (`4350` vs `4319`),
   so the ordering fix should be reported as native/Wasm determinism work, not
   as full TBLASTX NCBI parity.
8. If an NCBI oracle trace shows that the libc `qsort` used by NCBI reorders a
   specific comparator-equal group in a way that differs from input-order
   preservation, model that exact NCBI-observed order with source-backed
   evidence. Until then, preserving list order for comparator-equal HSPs is the
   least invasive deterministic emulation of the NCBI comparator contract.
