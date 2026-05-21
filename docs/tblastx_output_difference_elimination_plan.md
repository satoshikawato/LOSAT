# TBLASTX Output Difference Elimination Plan

Goal: eliminate the remaining observable output differences among NCBI
BLAST+ TBLASTX, LOSAT native TLOSATX, and LOSAT Wasm TLOSATX for the
`TrcuMJNV_vs_MellatMJNV` comparison, without introducing behavior that is not
backed by NCBI BLAST source.

This plan is scoped to output-order and Wasm/native determinism defects. It is
not a license to change scoring, masking, chaining, filtering, or formatting
unless the first divergent stage proves that NCBI does the same thing.

## Current Evidence

Case:

```bash
BLAST+:
  LOSAT/tests/blast_out/TrcuMJNV.MellatMJNV.tblastx.n8.out

LOSAT native:
  LOSAT/tests/losat_out/TrcuMJNV.MellatMJNV.tlosatx.n8.out

LOSAT Wasm:
  LOSAT/tests/losat_out/TrcuMJNV.MellatMJNV.tlosatx.wasm.n8.out
```

Native LOSAT vs BLAST+:

```text
ncbi_line_count                                  5804
native_line_count                                5804
exact_full_line_matches_after_sorting            5804
matches_ignoring_only_evalue                     5804
matches_ignoring_evalue_and_bitscore             5804
shared_hsp_evalue_mismatch_rows                  0
shared_hsp_bitscore_mismatch_rows                0
ncbi_only_hsp_keys_ignoring_evalue               0
native_only_hsp_keys_ignoring_evalue             0
ncbi_only_alignment_keys_ignoring_evalue+score   0
native_only_alignment_keys_ignoring_evalue+score 0
```

Interpretation: native LOSAT and BLAST+ already have identical hit content for
this case after sorting. The visible/native difference is output order, not
HSP set, E-value, bit score, identity, or coordinates.

Wasm LOSAT vs BLAST+ / native:

```text
ncbi_line_count                                  5804
wasm_line_count                                  5749
exact_full_line_matches_after_sorting            3127
matches_ignoring_only_evalue                     5490
matches_ignoring_evalue_and_bitscore             5689
shared_hsp_evalue_mismatch_rows                  2363
shared_hsp_bitscore_mismatch_rows                199
ncbi_only_hsp_keys_ignoring_evalue               314
wasm_only_hsp_keys_ignoring_evalue               259
ncbi_only_alignment_keys_ignoring_evalue+score   115
wasm_only_alignment_keys_ignoring_evalue+score   60
alignment-key diff bit-score floor range         22-30
```

Interpretation: Wasm still has real content differences. The largest bucket is
shared HSPs with different E-values, with a smaller number of bit-score and
HSP-set differences.

## NCBI Source Anchors

Only the NCBI C/C++ implementation is authoritative.

- `ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1330-1356`
  defines `ScoreCompareHSPs`.
- `ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1374-1381`
  performs guarded `qsort(..., ScoreCompareHSPs)` in
  `Blast_HSPListSortByScore`.
- `ncbi-blast/c++/src/algo/blast/core/link_hsps.c:483-486`
  sorts translated-query link HSPs with `s_RevCompareHSPsTbx`.
- `ncbi-blast/c++/src/algo/blast/core/link_hsps.c:990-994`
  sorts translated-query replay order with `s_RevCompareHSPsTransl` and
  `s_FwdCompareHSPsTransl`.
- `ncbi-blast/c++/src/algo/blast/core/link_hsps.c:1018-1059`
  skips non-start chain members and replays linked chains for output.
- `ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3078-3107`
  defines `s_EvalueCompareHSPLists`.
- `ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3331-3338`
  sorts hit lists by E-value order.

Important property: these NCBI comparators are partial for distinct TBLASTX
HSPs. They do not include every LOSAT-visible field, such as query frame,
subject frame, seed identity, or original insertion position. Comparator-equal
rows are therefore delegated to the platform `qsort` behavior in NCBI.

## Phase 0: Preserve Current Diagnostics

Save reproducible classification outputs before changing code:

```bash
python tests/classify_tblastx_outfmt6.py \
  --ncbi LOSAT/tests/blast_out/TrcuMJNV.MellatMJNV.tblastx.n8.out \
  --native LOSAT/tests/losat_out/TrcuMJNV.MellatMJNV.tlosatx.n8.out \
  --out-dir .tmp/trcu_mellat_native_classify

python tests/classify_tblastx_outfmt6.py \
  --ncbi LOSAT/tests/blast_out/TrcuMJNV.MellatMJNV.tblastx.n8.out \
  --native LOSAT/tests/losat_out/TrcuMJNV.MellatMJNV.tlosatx.wasm.n8.out \
  --out-dir .tmp/trcu_mellat_wasm_classify
```

Acceptance:

- Native sorted full-line match remains `5804/5804`.
- Wasm baseline summary is recorded before repair.

## Phase 1: Fix Native Raw Output Order

Native has matching content, so do not touch extension, reevaluation, masking,
linking, or scoring for this phase.

Target files:

- `LOSAT/src/common.rs`
- If needed, a small native qsort helper for final `Hit` / subject-group
  ordering, modeled after `LOSAT/src/algorithm/tblastx/ncbi_qsort.rs`.

Tasks:

1. Remove nondeterministic final-output grouping order.
   - Current final output uses `HashSet` while collecting subject indices in
     `write_output_ncbi_order_with_format_to_writer`.
   - Replace it with deterministic query/subject traversal that follows the
     NCBI hit-list lifecycle.

2. Audit final HSP sorting against `ScoreCompareHSPs`.
   - The comparator must remain exactly score DESC, subject offset ASC,
     subject end DESC, query offset ASC, query end DESC.
   - Do not add frame, chain, insertion-order, or formatted-coordinate
     tie-breakers unless an NCBI source or oracle trace proves such ordering.

3. Audit final subject ordering against `s_EvalueCompareHSPLists`.
   - Best E-value ASC, with the NCBI near-zero comparison.
   - Best score DESC.
   - Subject OID DESC on ties.

4. If raw order still differs after deterministic traversal, reproduce NCBI
   final `qsort` tie behavior for native final output using a C `qsort` index
   wrapper, not Rust-only tie-breakers.

Acceptance:

```bash
cmp -s \
  LOSAT/tests/blast_out/TrcuMJNV.MellatMJNV.tblastx.n8.out \
  LOSAT/tests/losat_out/TrcuMJNV.MellatMJNV.tlosatx.n8.out
```

After regenerating outputs from the same commands, raw native output should be
byte-identical to BLAST+ for this case.

## Phase 2: Fix Wasm Qsort Tie Behavior Before Linking

Wasm content differences must be fixed before final output order is tuned.
The likely offender is that native can route comparator-equal HSP ordering
through platform C `qsort`, while Wasm currently uses Rust `sort_by`.

Target files:

- `LOSAT/src/algorithm/tblastx/ncbi_qsort.rs`
- `LOSAT/src/algorithm/tblastx/blast_gapalign.rs`
- `LOSAT/src/algorithm/tblastx/sum_stats_linking/linking.rs`
- `LOSAT/src/algorithm/tblastx/blast_engine/run_impl.rs`
- `LOSAT/src/algorithm/tblastx/stage_dump.rs`

Tasks:

1. Add stage snapshots for this case if they are not already sufficient:
   - after initial ungapped extension
   - after common-endpoint purge
   - after reevaluation
   - before first link HSP pass
   - after link HSP pass
   - after e-value assignment
   - before final output sort

2. Compare native vs Wasm snapshots using stable internal keys, not formatted
   nucleotide coordinates.

3. Identify the first stage where sorted content diverges.
   - If sorted content first diverges before linking, investigate extension,
     reevaluation, or masking with NCBI source.
   - If sorted content matches before linking but diverges after linking, focus
     on qsort tie behavior and link-order replay.

4. Implement a Wasm-compatible emulation of the native/NCBI qsort tie order for
   the relevant TBLASTX HSP arrays.
   - Do not call NCBI BLAST from LOSAT runtime.
   - Do not add semantic tie-breakers absent from NCBI comparators.
   - If exact platform qsort behavior cannot be reproduced generically, model
     only the observed NCBI/native ordering boundary with source-backed
     evidence and regression tests.

5. Apply the same qsort emulation to all NCBI-equivalent translated-query sort
   boundaries:
   - `ScoreCompareHSPs`
   - `s_RevCompareHSPsTbx`
   - `s_RevCompareHSPsTransl`
   - `s_FwdCompareHSPsTransl`

Acceptance:

```bash
python tests/classify_tblastx_outfmt6.py \
  --ncbi LOSAT/tests/blast_out/TrcuMJNV.MellatMJNV.tblastx.n8.out \
  --native LOSAT/tests/losat_out/TrcuMJNV.MellatMJNV.tlosatx.wasm.n8.out \
  --out-dir .tmp/trcu_mellat_wasm_after
```

Expected summary:

```text
native_line_count                                5804
exact_full_line_matches_after_sorting            5804
matches_ignoring_only_evalue                     5804
matches_ignoring_evalue_and_bitscore             5804
shared_hsp_evalue_mismatch_rows                  0
shared_hsp_bitscore_mismatch_rows                0
ncbi_only_hsp_keys_ignoring_evalue               0
native_only_hsp_keys_ignoring_evalue             0
ncbi_only_alignment_keys_ignoring_evalue+score   0
native_only_alignment_keys_ignoring_evalue+score 0
```

## Phase 3: Final Wasm Raw Output Order

Only after Phase 2 sorted content is identical, tune Wasm raw output order.

Tasks:

1. Compare Wasm raw output to native raw output.
2. If sorted outputs match but raw outputs differ, reuse the final-output order
   repair from Phase 1.
3. Verify that the final order does not depend on `HashMap` or `HashSet`
   iteration order.

Acceptance:

```bash
cmp -s \
  LOSAT/tests/losat_out/TrcuMJNV.MellatMJNV.tlosatx.n8.out \
  LOSAT/tests/losat_out/TrcuMJNV.MellatMJNV.tlosatx.wasm.n8.out
```

After regeneration, native and Wasm output should be byte-identical.

## Regression Matrix

Run these only after the relevant phase is implemented:

```bash
bash tests/compare_tblastx_wasm_parity.sh \
  .tmp/tblastx_wasm_parity_after_trcu_fix

bash tests/compare_tblastx_native_ncbi_parity.sh \
  .tmp/tblastx_native_ncbi_after_trcu_fix
```

Also re-run at least:

- `LC738874.fasta` vs `LC738875.fasta`, default `--outfmt 6 -n 1`.
- `AP027131.fasta` vs `AP027133.fasta`, `--query-gencode 4 --db-gencode 4`.
- `TrcuMJNV.fasta` vs `MellatMJNV.fasta`, `--query-gencode 1 --db-gencode 1`.

## Non-Goals

- Do not delegate runtime behavior to NCBI BLAST executables or libraries.
- Do not change scoring, masking, genetic-code handling, or E-value formulas
  merely to reduce this output-order issue.
- Do not add LOSAT-only features or tie-breakers that are absent from NCBI.
- Do not treat the project-approved local-subject non-default `db_gencode`
  exception as permission for any other BLAST+ deviation.

## Done Criteria

The work is complete when:

1. Native LOSAT and BLAST+ are byte-identical for
   `TrcuMJNV_vs_MellatMJNV` after regenerating outputs.
2. Wasm LOSAT and native LOSAT are byte-identical for the same case.
3. The classifier reports zero HSP-key, alignment-key, E-value, and bit-score
   differences for Wasm vs BLAST+.
4. The short and long TBLASTX regression fixtures listed above still pass.
5. Every code change made for the plan includes NCBI C/C++ reference comments
   with file paths and line numbers immediately above the ported Rust logic.
