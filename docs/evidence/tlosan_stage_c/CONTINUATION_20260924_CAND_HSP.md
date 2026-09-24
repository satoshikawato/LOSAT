# TLOSAN Stage C continuation: long candidate order and real HSPTest deletion

Starting branch and commit: `feature/tlosan-tblastn-v0.2.0` at
`f9832b6ce5fc8cb8dbde16286542186525819786`, clean worktree. Read the
[preceding checkpoint](CONTINUATION_20260924_9077.md) and its linked records.
Pinned NCBI source: `598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4`.
Comparison-only `tblastn` SHA-256:
`e3956f1e107a30439d56c8f72fae4267a7d62ebcc16e4f99d4baf2bdf00402e0`.

**Stage C remains open. D and E have not begun. This record makes no output,
performance, completion, release, or certification claim.** The public TBLASTN
CLI still returns the explicit Stage C-E unimplemented error.

## Exact long-subject candidate streams

The [retained compact fixture](ordered_candidates_20260924/) holds every
ordered NCBI scan callback candidate for four existing generated subjects.
The unmasked and boundary-masked subjects are 15,000,962 nt, or 5,000,320
translated +1 residues; each has 12 WordFinder calls. The two middle-mask
subjects are 30,001,262 nt, have 15 surviving WordFinder calls, and exercise
an omitted positive middle chunk. The two-hit variant has matches on both
sides of the masked span. Inputs, exact options, commands, hashes, frame/chunk
states, and generation rules are in the manifests and runners. No long FASTA
is committed.

| Fixture | NCBI ordered candidates | Rust comparison |
| --- | ---: | --- |
| `long_chunk_20260924` | 259 | all 259 match |
| `masked_chunk_boundary_20260924` | 113 | all 113 match |
| `no_range_middle_20260924` | 189 | all 189 match |
| `no_range_two_hits_20260924` | 381 | all 381 match |

Pinned `blast_engine.c:221-325,478-500,804-841` supplies the strict
5,000,000-residue split, 100-residue overlap, omitted no-range calls, frame
order, and WordFinder timing. `aa_ungapped.c:478-505` calls `scansub` for each
range; the comparison-only `ncbi_candidate_trace.c` records every callback
pair before two-hit extension. Rust's `search_seed.rs` comparison reconstructs
the four exact inputs, checks the 12/15 call states against saved NCBI ranges,
and compares frame, chunk offset, query offset, subject offset, and global
candidate order without filtering. Each probe preserved the unprobed NCBI
final bytes. No new candidate difference was found in these four profiles.
This does not cover multiple query contexts on a long subject.

## Positive HSPTest deletion in a real traceback call

The [real-path fixture](hsp_test_real_path_20260924/) retains the raw probe
stream as deterministic gzip and its ordered event rows as TSV. It reuses the saved
three-query/6,377-nt subject. NCBI `blast_traceback.c:585-605` updates an HSP,
calculates alignment identity and length, calls `Blast_HSPTest`, and frees an
HSP when the result is true. `blast_hits.c:993-1001` compares identity against
`BlastHitSavingOptions.percent_identity`. The comparison-only
`ncbi_hsp_test_inject.c` sets that option to the valid C API value 100.0 only
during each actual `Blast_HSPTest` call, then restores it. It is preloaded
before the unmodified `ncbi_gapped_trace.c` probe. This is a characterized
intervention in NCBI input state, not a CLI/default-search parity oracle.

The real NCBI traceback made 31 HSPTest calls, deleted 17 HSPs, and changed its
output from 3,858 to 2,886 bytes. Rust's internal diagnostic received the
same preliminary lists and 100.0 identity option. Its 31 ordered call records
match NCBI in delete flag, alignment length, query start/end, and subject
start/end; its final surviving HSP lists match both query contexts after the
fence retry. The test observes events at the same point as NCBI, after
traceback and before subject-offset adjustment or interval-tree insertion.
The comparison-only NCBI code and libraries are absent from LOSAT runtime and
build paths.

## Reproduce and verification

From the repository root, using fresh output directories:

    python3 docs/evidence/tlosan_stage_c/run_long_chunk_trace.py /tmp/tlosan-c-long-new
    python3 docs/evidence/tlosan_stage_c/run_long_chunk_trace.py /tmp/tlosan-c-mask-new --mask-boundary
    python3 docs/evidence/tlosan_stage_c/run_no_range_middle_trace.py /tmp/tlosan-c-middle-new
    python3 docs/evidence/tlosan_stage_c/run_no_range_middle_trace.py /tmp/tlosan-c-two-new --early-hit
    python3 docs/evidence/tlosan_stage_c/run_hsp_test_deletion_trace.py /tmp/tlosan-c-hsp-test-new
    (cd docs/evidence/tlosan_stage_c/ordered_candidates_20260924 && sha256sum -c retained.sha256)
    (cd docs/evidence/tlosan_stage_c/hsp_test_real_path_20260924 && sha256sum -c retained.sha256)
    cargo test --manifest-path LOSAT/Cargo.toml --lib long_subject_every_candidate_matches_ncbi_in_chunk_order
    cargo test --manifest-path LOSAT/Cargo.toml --lib real_path_hsp_test_deletion_order_matches_ncbi

The four full NCBI candidate traces were regenerated in fresh `/tmp`
directories and matched their previously retained final-output bytes. The
candidate Rust test passed. The HSPTest runner recorded 31 calls and 17
positive deletions; the exact Rust test passed. The preexisting multi-query
traceback replay also passed with the optional observer disabled.
`cargo fmt --check`, `cargo clippy --lib --tests -- -D warnings`, both retained
SHA-256 manifests, all three runner syntax checks, `git diff --check`, and a
fresh debug build passed. The freshly built public CLI with `-db_gencode 32`
exited 1, wrote zero stdout bytes, and reported `TBLASTN local search is
unimplemented (Stages C-E)`. A whole-suite test, release build, benchmark,
and independent release audit were not run while Stage C remains open.

## Remaining Stage C gate

1. Produce a complete traceback path with positive start-offset acquisition
   failure and compare its actual NCBI input, deletion, and Rust order.
2. Produce positive **chunk-local** endpoint-purge deletion and a chunk-path
   `Blast_HSPListAppend` cap boundary; direct API controls and unlimited-cap
   snapshots do not close these cases.
3. Expand long-subject candidate and HSP comparison to multiple query
   contexts, including invalid contexts, ranking, and complete-path positive
   interval-tree containment. The existing short multi-query and replay
   fixtures remain bounded evidence.
4. Sweep all Stage C profiles and resolve any first difference with the pinned
   NCBI function, input state, and execution order before accepting C.

Only after C matches should D's function/input-state/order table address
statistics, linking, composition adjustment, Kappa redo, ranking, and deletion.
Only after C and D match should E byte-compare local `-subject` outfmt 0/6/7.
Code 32 requires the comparison-only NCBI C++ API oracle; `-db` statistics and
headers are not local-subject expected bytes.
