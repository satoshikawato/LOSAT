# LOSATN BLASTN Bitwise-Parity Implementation Plan

Status: in progress  
Owner: Codex  
Updated: 2026-08-10

## Objective

Make the supported pure-Rust LOSATN BLASTN path byte-identical to the NCBI
BLASTN oracle for the declared fixtures and options. NCBI C/C++ is the only
behavioral authority; NCBI BLAST+ is used only for comparison and diagnostics.

No runtime, build-time, FFI, or fallback dependency on NCBI BLAST may be added.
The only accepted parity exception is the existing local-subject TBLASTX
non-default `db_gencode` exception; it does not apply to BLASTN.

## Frozen current evidence

The current release binary was rebuilt before diagnosis. Fresh manifest runs
cover 11 cases from [`tests/blastn_parity_manifest.tsv`](../LOSAT/tests/blastn_parity_manifest.tsv):

- 10 cases match NCBI in hit count, coordinate keys, identity, alignment
  length, mismatches, gap opens, E-value, and bit score.
- `Sakai.MG1655.megablast` has 6476 common HSPs and five same-coordinate
  statistic differences. Bit scores and E-values still match.
- A traced residual first diverges in the greedy traceback edit script after
  seed, ungapped extension, and preliminary gapped coordinates/raw score
  match.
- The active BLASTN writer in [`src/algorithm/blastn/hsp.rs`](../LOSAT/src/algorithm/blastn/hsp.rs)
  writes data rows only. The generic outfmt-7 header writer is not used by the
  BLASTN engine, so raw outfmt-7 bytes cannot match NCBI.

The comparison helper is diagnostic, not a bitwise gate: it ignores comment
lines, collapses duplicate coordinate keys, and does not validate raw bytes.

## Scope and non-goals

In scope:

1. Single-query/single-subject `blastn --task blastn` and default megablast
   cases in the current manifest.
2. NCBI-compatible outfmt 6/7 field values, ordering, headers, and newlines.
3. Greedy traceback, preliminary edit-block conversion, gap reduction, and
   post-traceback ambiguity re-evaluation needed by the residual megablast
   cases.
4. Native serial execution first; threaded/Wasm parity only after the native
   gate is exact.

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

- The five residual HSPs have identical edit-derived length, identity,
  mismatches, and gap opens.
- The focused megablast output has no structured field differences.
- No task-blastn fixture regresses.

## Phase 3 — Full native parity and regression gate

Work:

1. Run the complete current BLASTN manifest in a fresh temporary directory with
   a release build and fixed `num_threads=1`.
2. Require both raw-byte equality and structured equality for outfmt 6/7.
3. Re-run the existing focused Rust tests and the broader requested Rust gate.
4. Repeat the focused fixture to confirm deterministic output.
5. Record commands, source versions, checksums, first divergence, and residual
   risk in the evidence section below.

Acceptance evidence:

- All declared native serial fixtures are byte-identical.
- No accepted exception is used for BLASTN.
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

Status: source audit complete; residual binary-oracle discrepancy unresolved  
Behavior implemented: The active greedy traceback and post-traceback path was
audited against `greedy_align.c`, `blast_gapalign.c`, `blast_hits.c`, and
`blast_traceback.c`. Temporary internal tracing was used to compare the first
residual HSP, then removed. No reducer or score/statistics workaround was
introduced because it would not be supported by the NCBI source.  
Evidence: For the first residual HSP, seed, ungapped extension, preliminary
score/coordinates, forward and reverse preliminary blocks, merged edit script,
`s_ReduceGaps`, and ambiguity re-evaluation all match the current NCBI
2.17.0+ source path. The installed NCBI 2.17.0+ binary (build Aug 11 2025)
still emits a score-equivalent alternative gap placement on five
Sakai/MG1655 HSPs. The current NCBI 2.17.0+ source snapshot used for the
audit is dated Dec 21 2025; the source-level C harness produced the same edit
script as LOSAT for the traced HSP, including after the NCBI
score-only-then-traceback memory reuse sequence.  
Deviation: five installed-binary residuals remain: length/mismatch/gap-open
differences on 5 rows and pident differences on 2 rows; coordinates, hit
count, bit score, and E-value remain exact.  
Remaining risk: matching the older installed binary would require an
authoritative source or instrumented binary showing its tie choice. A
post-hoc score-equivalent alignment preference would violate the NCBI-only
policy and is intentionally not applied.

### Phase 3

Status: in progress  
Behavior implemented: The comparison helper now has a raw-byte gate in
addition to its structured diagnostic diff.  
Evidence: fresh manifest refresh: all 11 cases have identical hit counts and
coordinate keys; 10 cases have exact structured values; Sakai/MG1655 has the
five residual edit-derived statistic differences documented in Phase 2.
`env CARGO_TARGET_DIR=/tmp/losatn-parity-test-target cargo test --release`
passes (392 passed, 1 ignored in library tests; 5 CLI tests; 4 compressed
lookup tests; 177 integration tests, 2 ignored). `cargo clippy --release
--all-targets` and the release build also pass.  
Deviation: the full raw-byte manifest gate is not yet a valid acceptance claim
because the checked-in NCBI output fixtures were generated with different
relative subject paths; the fresh same-argument focused fixture is byte exact.
 
Remaining risk: run the full Rust test gate and add fresh paired raw-byte
artifacts before declaring native serial parity complete. `cargo fmt --
--check` remains non-clean because the already-dirty `run.rs` contains
unrelated formatting drift; no whole-file formatting rewrite was applied.
Deviation: none  
Remaining risk: native parity must pass before any Wasm/threading claim.
