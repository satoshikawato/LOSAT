# TLOSAN Stage C continuation: traceback start-offset failure and success

Starting branch and commit: `feature/tlosan-tblastn-v0.2.0` at `b0ea2300`.
Read [the preceding bounded checkpoint](CONTINUATION_20260924_CAND_HSP.md)
and its linked records. Pinned NCBI source:
`598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4`; comparison-only
`tblastn` SHA-256:
`e3956f1e107a30439d56c8f72fae4267a7d62ebcc16e4f99d4baf2bdf00402e0`.

**Stage C remains open. D and E have not begun.** No TBLASTN output,
performance, completion, release, or certification claim is made. The public
TBLASTN CLI retains the explicit Stage C-E unimplemented error.

## First difference and NCBI owner

The [failure fixture](start_failure_real_path_20260924/) and
[success fixture](start_success_real_path_20260924/) reuse the pinned
three-query, 6,377-nt local `-subject` search. A comparison-only probe changes
the first preliminary HSP at each query-index 0 traceback entry to query
`0..20`, subject `+1 0..20`, and both gapped starts `0`. It preserves the
pre-sorted score, so this is an artificial function-input state, not a natural
search result or a CLI parity oracle. Unprobed NCBI bytes still match the
previous saved baseline. Each modified run has its exact input, command, source/binary hashes,
complete raw trace in deterministic gzip, ordered event TSV, output, and
SHA-256 manifest retained.

The first reproducible Rust difference was at
`blast_traceback.c:436-445`: NCBI called
`BlastGetOffsetsForGappedAlignment` before traceback when both gapped starts
were zero. Rust's diagnostic entered gapped traceback directly. The pinned
`blast_gapalign.c:3248-3321` computes an 11-residue sliding-window score,
then tests the terminal window. Both are nonpositive for this HSP, so NCBI
returned `FALSE` twice: on the partial translation attempt and after the
full-subject fence retry. It freed the HSP before `Blast_HSPTest` on both
passes. The second successful NCBI traceback call ended with 16 HSPs for
query index 0, after endpoint filtering.

Rust now ports that exact window calculation for the protein matrix and calls
it at the same traceback point. The internal diagnostic records both
`FALSE` inputs in the same order and drops the HSP before extension. With the
same artificial first HSP, the final ordered HSP lists match both query
contexts after the fence retry. The score calculation uses the existing
NCBI-derived standard protein matrix; no NCBI code is linked or invoked from
LOSAT runtime or build paths.

A second comparison-only input changes that HSP to query `0..120`, subject
`+1 200..320`, with both gapped starts zero. At the same NCBI real-path call,
the sliding-window function returns `TRUE` with `(34,234)` on both passes.
The first partial pass hits a fence; on the full retry NCBI writes these
acquired starts into the surviving HSP. Rust matches both function rows,
the final first HSP including its acquired starts, and the ordered final
HSP lists for both query contexts. The final NCBI output bytes equal the
unprobed baseline in this success case. This is an artificial input-state
check; natural candidate coverage is still a Stage C gate.

## Reproduce and verification

From the repository root, with a fresh output directory:

    python3 docs/evidence/tlosan_stage_c/run_start_failure_trace.py /tmp/tlosan-c-start-fail-new
    python3 docs/evidence/tlosan_stage_c/run_start_failure_trace.py /tmp/tlosan-c-start-pass-new positive
    (cd docs/evidence/tlosan_stage_c/start_failure_real_path_20260924 && sha256sum -c retained.sha256)
    (cd docs/evidence/tlosan_stage_c/start_success_real_path_20260924 && sha256sum -c retained.sha256)
    cargo test --manifest-path LOSAT/Cargo.toml --lib algorithm::tblastn::search_gapped::tests

The runner requires the pinned source commit and binary checksum, confirms
plain output against the saved baseline, and records exactly two start results
of the selected kind. The failure output differs; the success output is
unchanged. The Rust tests compare NCBI function-input/result rows, both
final query-indexed HSP lists, and the successful HSP start writeback.

All five retained files in each fixture reproduced byte for byte in fresh
`/tmp` runs; both SHA-256 manifests passed.
The affected `algorithm::tblastn::search_gapped::tests` suite passed 23/23,
including both real-path start fixtures, direct short/sliding/terminal-window
boundary tests, the earlier real-path HSPTest fixture, and the default
traceback regressions. `cargo fmt --check`,
`cargo clippy --lib --tests -- -D warnings`, the retained SHA-256 manifest,
runner syntax, and `git diff --check` passed. A full suite, release build,
benchmark, and release audit are deferred while Stage C remains open.

## Remaining Stage C gate

1. Produce positive **chunk-local** endpoint-purge deletion and a chunk-path
   `Blast_HSPListAppend` cap boundary, comparing exact NCBI input state,
   execution order, and removed HSPs. Existing direct API cap tests and
   unlimited-cap append snapshots remain bounded.
2. Expand long-subject candidate/HSP comparison to multiple query contexts,
   including invalid contexts and ranking, and establish positive
   interval-tree containment on a complete-path NCBI fixture. The current
   short multi-query, alternate-matrix, and traceback replay fixtures do not
   by themselves close the complete gate.
3. Sweep every Stage C fixture and resolve the first difference, if any,
   against the pinned NCBI function/state/order before accepting C.

Only after C is exact should D's function/state/order table cover statistics,
linking, composition adjustment, Kappa redo, ranking, and deletion. Only after
C and D match should E byte-compare local `-subject` outfmt 0/6/7. Code 32
requires the comparison-only NCBI C++ API oracle; `-db` statistics and
headers are not local-subject expected bytes.
