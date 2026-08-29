# LOSATN BLASTN Bitwise-Parity Implementation Plan

Status: in progress  
Owner: Codex  
Updated: 2026-08-29

## Objective

Make the supported pure-Rust LOSATN BLASTN path exact for source-defined NCBI
BLASTN behavior and source-compatible for the narrowly documented
source-underdetermined equal-HSP tie class. NCBI C/C++ is the only behavioral
authority; NCBI BLAST+ is used only for comparison and diagnostics.

No runtime, build-time, FFI, or fallback dependency on NCBI BLAST may be added.
The only accepted behavioral deviation is the existing local-subject TBLASTX
non-default `db_gencode` exception; it does not apply to BLASTN. The BLASTN
source-underdetermined tie policy is not a behavioral deviation: it preserves
the set of survivor representations permitted by NCBI source.

## Frozen current evidence

The current native release binary is pinned to LOSAT commit
`92b934d8defd299bdf821055d92053b8f8923fad` and SHA-256
`d87f8c210d3506bd75ac9b0a5d7628a8c64bcb9c26f0b99923c2742fea366582`.
Fresh manifest runs cover all 14 cases from
[`tests/blastn_parity_manifest.tsv`](../LOSAT/tests/blastn_parity_manifest.tsv):

- 13 cases are byte-identical and match NCBI in hit count, coordinate keys,
  identity, alignment length, mismatches, gap opens, E-value, and bit score.
- `Sakai.MG1655.megablast` has 6476 common HSPs and five same-coordinate
  statistic differences. Bit scores and E-values still match.
- The installed binary's first residual differs from LOSAT in the
  edit-script-derived statistics, but a direct current-source trace of all
  five residual HSPs produces the same greedy scripts as LOSAT.
- The active BLASTN writer in [`src/algorithm/blastn/hsp.rs`](../LOSAT/src/algorithm/blastn/hsp.rs)
  emits NCBI-compatible outfmt-6/7 rows and headers, and
  [`tests/compare_blastn_parity.py`](../LOSAT/tests/compare_blastn_parity.py)
  checks both structured fields and raw bytes for fresh paired runs.

The comparison helper retains its structured diagnostic view (comment lines are
ignored and duplicate coordinate keys are reported) and now independently
checks raw bytes for fresh paired runs.

## Compatibility contract for source-underdetermined ties

Source-defined BLASTN cases require exact official NCBI BLAST+ 2.17.0 binary
parity. When distinct HSP edit scripts tie on every field in the NCBI
common-endpoint comparator, the NCBI source does not define which edit script
must survive its `qsort` and first-survivor purge. Under
[`PD-BLASTN-HSP-CANONICALIZATION`](product_decisions/PD-BLASTN-HSP-CANONICALIZATION.md)
Version 1.1, LOSAT preserves that source semantics and does not add a secondary
total order solely to reproduce one precompiled binary's survivor.

For the current 14-case manifest, 13 source-defined cases require exact parity.
The known source-underdetermined case may be accepted only with the same HSP row
count, HSP membership and coordinate keys, E-values, and bit scores, and only
when differences are confined to the observed edit-script-derived `pident`,
alignment-length, mismatch, and gap-open fields. Coordinate, HSP-set, score,
E-value, formatter, and unrelated algorithm differences remain failures.

No cross-platform deterministic survivor is claimed. Supported native/Wasm
executions are to be monitored and tested; a LOSAT-specific canonicalization
may be reconsidered only after evidence of a material supported-platform
incompatibility.

## Scope and non-goals

In scope:

1. Single-query/single-subject `blastn --task blastn` and default megablast
   cases in the current manifest.
2. NCBI-compatible outfmt 6/7 field values, ordering, headers, and newlines.
3. Greedy traceback, preliminary edit-block conversion, gap reduction, and
   post-traceback ambiguity re-evaluation needed by the residual megablast
   cases.
4. Native serial execution first; threaded/Wasm evidence only after the native
   source-compatibility gate is satisfied.

Out of scope:

- Replacing an unresolved Rust path with an NCBI executable or library.
- Generic alignment heuristics, score-equivalent alignment selection, or
  changing `gapopen` statistics independently of the NCBI edit script.
- Unrelated dirty-tree files and the older broad plan in
  `docs/blastn_bit_parity_fix_plan.md`.

## Phase 1 — Make the parity gate truthful

Owners:

- `tests/compare_blastn_parity.py`
- `src/algorithm/blastn/hsp.rs`
- relevant BLASTN output tests

Work:

1. Preserve the existing structured diff for diagnosis.
2. Add a raw-byte comparison gate for each selected outfmt.
3. Route the active BLASTN outfmt-7 path through an NCBI-compatible header
   and data-row writer, preserving query/subject descriptions, field names,
   hit count, ordering, and final newline behavior.
4. Keep outfmt 6 free of comment headers.

NCBI references:

- `c++/src/objtools/align_format/tabular.cpp:1100-1108` — data-row emission.
- `c++/src/objtools/align_format/tabular.cpp:1264-1284` — outfmt-7 header.
- `c++/src/algo/blast/format/blast_format.cpp:770-782` — format dispatch.

Acceptance evidence:

- A focused outfmt-7 fixture has identical SHA-256 and `cmp` exit status 0.
- The same fixture has identical parsed rows and ordering.
- A no-hit outfmt-7 case has identical header and no-hit text.

## Phase 2 — Isolate and correct greedy traceback parity

Owners:

- `src/algorithm/blastn/alignment/greedy.rs`
- `src/algorithm/blastn/blast_engine/run.rs`
- targeted BLASTN alignment tests

NCBI owners:

- `c++/src/algo/blast/core/greedy_align.c:379-751` — non-affine greedy
  traceback and tie order.
- `c++/src/algo/blast/core/blast_gapalign.c:2481-2536` — preliminary block
  conversion.
- `c++/src/algo/blast/core/blast_gapalign.c:2572-2758` — edit-script update,
  rebuild, and `s_ReduceGaps`.
- `c++/src/algo/blast/core/blast_gapalign.c:2762-2936` — forward/reverse
  traceback timing, sequence representation, and script construction.
- `c++/src/algo/blast/core/blast_traceback.c:633-692` — post-traceback
  re-evaluation and purge timing.

Work:

1. Freeze `Sakai/MG1655` row 986 as the first current residual and retain the
   five residual rows as regression targets.
2. Compare, in order, forward prelim block, reverse prelim block, merged edit
   script, post-`s_ReduceGaps` script, and post-re-evaluation script.
3. Verify exact NCBI tie behavior, row origins, sentinel values, pool rewind and
   stale-cell preservation, packed score-only versus uncompressed traceback
   input, and q/s boundary pointers.
4. Fix only the first divergent stage. Do not adjust statistics after the fact.
5. Add small deterministic tests for repeated bases, zero gap penalties,
   minus-strand output, and boundary-length cases.

Acceptance evidence:

- The five residual HSPs match the NCBI 2.17.0 source path through final greedy
  edit-script construction.
- The precompiled-binary residual is classified at the source-underdetermined
  common-endpoint survivor boundary, with row count, membership/coordinate
  keys, E-values, and bit scores unchanged.
- No task-blastn fixture regresses.

## Phase 3 — Full native source-compatibility and regression gate

Work:

1. Run the complete current BLASTN manifest in a fresh temporary directory with
   a release build and fixed `num_threads=1`.
2. Require raw-byte and structured equality for source-defined outfmt 6/7
   cases; apply only the fixed invariants to the known source-underdetermined
   tie.
3. Re-run the existing focused Rust tests and the broader requested Rust gate.
4. Repeat the focused fixture to monitor survivor stability on the tested
   target, without inferring a cross-platform deterministic contract.
5. Record commands, source versions, checksums, first divergence, and residual
   risk in the evidence section below.

Acceptance evidence:

- All 13 source-defined native serial fixtures are byte-identical.
- The known source-underdetermined case satisfies the fixed row-count,
  membership/coordinate, E-value, bit-score, and differing-field invariants.
- No BLASTN deviation from source-defined behavior is used.
- Any unsupported option remains an explicit Rust error rather than an NCBI
  fallback.

## Execution record

### Phase 1

Status: completed  
Behavior implemented: The active BLASTN path now parses outfmt 6/7, emits
NCBI-compatible outfmt-7 headers and footer, preserves query deflines and the
user-supplied subject label, and fails explicitly for unported custom fields.  
Evidence: `env CARGO_TARGET_DIR=/tmp/losatn-parity-test-target cargo test
--release blastn::hsp::tests::test_blastn_outfmt7_header_matches_ncbi_shape`
passed; a fresh `LC738869/AP027202 task=blastn outfmt=7` comparison using the
same relative input paths returned `cmp` exit 0 with both files 16106 bytes,
214 lines, and SHA-256
`30717748c1eac9b994512c1117eb440d1c5a2ec7f05f45dcfe3d50bd3ca13177`.  
Deviation: custom field specifications remain explicitly unsupported and now
fail fast instead of being silently treated as the default field set.  
Remaining risk: the no-hit `AP027078/AvCLPV` fixture was verified fresh but is
not yet part of the permanent manifest; multi-query raw-byte coverage also
remains to be added. The fresh no-hit run returned `cmp` exit 0 with the
expected `# 0 hits found` header.

### Phase 2

Status: current-source traceback parity confirmed; installed-binary discrepancy
remains unresolved
Behavior implemented: The active greedy traceback and post-traceback path was
audited against `greedy_align.c`, `blast_gapalign.c`, `blast_hits.c`, and
`blast_traceback.c`. Temporary internal tracing was used to compare the first
residual HSP, then removed. NCBI-referenced regression tests now cover
repeated-base tie cases, zero-gap-penalty indels, affine traceback, reverse
greedy traversal, empty/one-base boundaries, and minus-query coordinate
conversion. No reducer or score/statistics workaround was introduced because
it would not be supported by the NCBI source.
Evidence: For all five residual HSPs, seed, ungapped extension, preliminary
score/coordinates, forward and reverse preliminary blocks, merged edit script,
`s_ReduceGaps`, and final greedy edit script were compared with the current
NCBI source implementation. A corrected source-level C harness calling
`BLAST_GreedyGappedAlignment` produced these `full_ops` values, exactly
matching LOSAT's traces: row 1
`S152,I1,S1,I1,S8,D2,S68`; row 2
`S30,D2,S1,D1,S3,I1,S1,I2,S66`; row 3
`S30,D2,S4,I1,S3,I1,S19,I1,S41`; row 4 `S139`; row 5
`S34,D1,S65,D1,S10`. The current NCBI source snapshot used for the audit is
commit `598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4` (2025-12-21). The installed
NCBI 2.17.0+ binary (build Aug 11 2025) still emits score-equivalent
alternative gap placements on the same five output rows. Focused isolated
tests pass:
`env CARGO_TARGET_DIR=/tmp/losatn-greedy-test-target cargo test --lib greedy
-- --nocapture` (2 passed) and
`env CARGO_TARGET_DIR=/tmp/losatn-greedy-test-target cargo test --lib
test_adjust_blastn_offsets_minus_query_keeps_internal_subject_order`
(1 passed).
Deviation: five installed-binary residuals remain: length/mismatch/gap-open
differences on 5 rows and pident differences on 2 rows; coordinates, hit
count, bit score, and E-value remain exact. These are not current-source
traceback defects.
Remaining risk: the August 2025 installed binary's exact source/provenance is
not available. A post-hoc score-equivalent alignment preference would violate
the NCBI-only policy and is intentionally not applied.

### Phase 3

Status: native serial source-compatibility contract satisfied; supported-target
monitoring remains
Behavior implemented: The comparison helper now has a raw-byte gate in
addition to its structured diagnostic diff, executes fresh paired NCBI/LOSAT
runs with identical relative input paths, and passes the supported `dust=true`
manifest value to LOSAT's flag-form CLI.
Evidence: release binary rebuilt, then the complete fresh native serial gate
was run:
`python tests/compare_blastn_parity.py --fresh-paired
--paired-output-dir /tmp/losat-blastn-manifest-20260810-full
--losat-bin target/release/LOSAT --ncbi-bin blastn`. NCBI was
`blastn 2.17.0+` (build Aug 11 2025 09:46:06), with binary SHA-256
`0b88d4a00cb7fa579c203653151175b24c83d27fe240d420abe7b5261a3083d1`;
LOSAT was commit `92b934d8defd299bdf821055d92053b8f8923fad` with binary
SHA-256 `d87f8c210d3506bd75ac9b0a5d7628a8c64bcb9c26f0b99923c2742fea366582`,
and Rust was `1.92.0`. Thirteen of the 14 cases have exact raw bytes and
structured rows. Sakai/MG1655 still has 6476/6476 rows and coordinate keys,
but five Phase 2 edit-script statistic residuals. Exact-output SHA-256
examples from the fresh pair are: `PesePMNV.MjPMNV.task_blastn`,
`cadccd5ba65c5632193b2866713328ca3689e684ab8718a773df8decf4d8284c`;
`PmeNMV.MjPMNV.task_blastn`,
`30717748c1eac9b994512c1117eb440d1c5a2ec7f05f45dcfe3d50bd3ca13177`;
`EDL933.Sakai.megablast`,
`4b4b2501f1d87303c3060c5b0c100fbd381604a72e4df130e23701ff81e49fac`;
`compact.multi_query.outfmt7`,
`58b0673dc9379427777880c709b9edb35647806e78de74197899680fac402547`.
`env CARGO_TARGET_DIR=/tmp/losatn-parity-test-target cargo test --release`
passes: 395 library tests passed with 1 ignored, 5 CLI tests passed, 4
compressed-lookup tests passed, and 177 integration tests passed with 2
ignored.  `env CARGO_TARGET_DIR=/tmp/losatn-parity-test-target cargo clippy
--release --all-targets` also passes.
Deviation: the raw-byte diagnostic is truthful and fresh; the only non-identical
output is the installed NCBI binary's five score-equivalent,
edit-script-derived statistic alternatives for Sakai/MG1655. No post-hoc
alignment preference was added because the inspected NCBI source does not
authorize it. `cargo fmt -- --check` remains non-clean because the existing
dirty `src/algorithm/blastn/blast_engine/run.rs` has unrelated formatting drift;
no whole-file formatting rewrite was applied.
Additional fresh Sakai/MG1655 verification on 2026-08-10 returned 6476 common
HSPs and the same five residual statistic differences. Raw output sizes were
536880 bytes each, with first byte difference at offset 82036; SHA-256 was
`2e15963c66e5f552088e468d17a1a86b779be61e7ed6331e57314f6e93c69ffa` for NCBI
and `9edd4883881316976c82c3e3674cab6d4d863e17ac89e3db2052b2ed7996df74` for
LOSAT. Under Product Decision Version 1.1, this is the known
source-underdetermined case: the other 13 manifest cases are byte-identical,
and this case retains identical row count, membership/coordinate keys,
E-values, and bit scores with differences limited to `pident`, alignment
length, mismatch, and gap-open fields. The native serial manifest therefore
satisfies the revised source-compatibility contract. No Wasm/threading survivor
claim is made.
