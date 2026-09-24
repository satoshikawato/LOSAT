# TLOSAN Stage C continuation: chunk deletions and long multi-query state

Branch `feature/tlosan-tblastn-v0.2.0`, resumed after `2261f8da`. Read the
[start-offset checkpoint](CONTINUATION_20260924_START_FAILURE.md) and its
predecessors. Pinned NCBI source is
`598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4`; comparison-only NCBI
`tblastn` SHA-256 is
`e3956f1e107a30439d56c8f72fae4267a7d62ebcc16e4f99d4baf2bdf00402e0`.

**Stage C remains open. D and E have not begun.** These are bounded native,
one-thread, local `-subject` checks, not complete TBLASTN, outfmt, speed,
release, or certification claims. The public TBLASTN CLI still has its
explicit Stage C-E unimplemented error.

## NCBI function, input state, and order

| Order and pinned source | Comparison fixture and exact observation | Rust check |
| --- | --- | --- |
| `aa_ungapped.c:478-505`, `blast_engine.c:804-841`: each translated frame/chunk calls `scansub` on the combined protein-query lookup | [Long multi-query](long_multi_query_20260924/) reconstructs a 15,000,962-nt subject (5,000,320 translated +1 residues) and uses three queries, with contexts 0/1 valid and 2 invalid. All 12 scan calls, 394 ordered candidate pairs, frame/chunk lengths, and per-call candidate counts were saved without an output filter. | `long_multi_query_every_candidate_and_context_matches_ncbi` compares every ordered pair and context; it checks all `PARAM`, `FRAME_CHUNK`, and `END` rows. |
| `aa_ungapped.c:562-590`: `BSearchContextInfo` selects `word_params->cutoffs + curr_context`; its x-drop and cutoff govern extension and HSP saving | The new comparison-only [WordFinder context probe](ncbi_wordfinder_context_cutoff_trace.c) records 36 rows across those 12 calls: valid contexts use `(16,16,28)` and `(16,16,25)` for initial x-drop, x-drop, cutoff; invalid context uses `(0,0,INT4_MAX)`. | The internal Rust WordFinder now indexes both x-drop and cutoff by query context. The long fixture passes the recorded arrays; an equality/just-above score boundary test covers separate contexts. |
| `aa_ungapped.c:562-590`: distinct context x-drop changes extension on the real WordFinder path | [Real-call x-drop intervention](word_xdrop_real_path_20260924/) sets context 1 x-drop from 16 to 1 during each of six actual calls, restoring 16 afterward. Its 32 ordered init HSP rows differ from the unmodified run, while final NCBI bytes remain the same. | `real_wordfinder_distinct_context_xdrop_matches_ncbi` compares every changed init HSP row in order using Rust x-drop array `[16,1,0]`. This input is artificial and comparison-only. |
| `blast_gapalign.c:3924-3927`, `blast_engine.c:478-586,840-850`: chunk-local gapped HSPs are purged, offset-adjusted, merged and appended in frame order | The same long input has context cutoffs 28/25/`INT4_MAX`, four ordered WordFinder init HSPs, four chunk-local gapped HSPs, two merged +1 HSPs, six append calls, and two query-indexed full-translation traceback retries. | `long_multi_query_chunk_hsps_and_append_match_ncbi` compares all four init HSPs before gapped scoring, all gapped HSP fields, final merged HSP gapped starts, every append snapshot, and ordered per-query post-retry HSPs. |
| `blast_engine.c:539-552,840-850`, `blast_hits.c:2455-2537,2809-2864`: endpoint purge precedes sort/merge and `Blast_HSPListAppend` | [Real-call intervention](chunk_purge_cap_real_path_20260924/) changes a lower-scoring HSP's starts at the first actual chunk purge, causing 3→2 removal; it passes cap 3 to each of the six actual frame append calls. The second append reduces 4 HSPs to 3; all six ordered snapshots are saved. | `real_chunk_endpoint_deletion_and_append_cap_match_ncbi` uses the same artificial function input and cap, compares the removed HSP and all append outputs. The unmodified NCBI command reproduces the earlier baseline bytes. |
| `blast_traceback.c:401-409,675-693`: interval-tree containment before traceback and again after endpoint purge/score sort | The natural BLOSUM45/word-2 fixture has 12 positive traceback-side containment events, including full-translation retry. | `alternate_matrix_traceback_positive_containment_matches_ncbi` now creates Rust HSPs from the natural query/subject through seed, gapped, chunk merge and append. It compares full append input HSPs including gapped starts, every containment predicate/input in call order, and final HSPs. |

The earlier start-offset and real `Blast_HSPTest` deletion probes remain in the
preceding records. Their interventions, like this purge/cap case, characterize
specific real NCBI call paths with artificial function-input states. No NCBI
source, executable or library is a LOSAT runtime or build dependency.

## First newly found input-state difference

The older long WordFinder tests passed cutoff 0 to Rust while NCBI's retained
`PARAM` rows used 28 for the 15,000,962-nt cases and 30 for the 30,001,262-nt
masked-middle cases. Strong HSPs masked that mismatch in the output. In the
new long three-query case, NCBI used cutoff 28 in context 0 and 25 in context
1. The pinned `aa_ungapped.c:562-590` selects both x-drop and cutoff by
`curr_context`. Rust previously used scalars; it now accepts arrays indexed at
the same extension point. The long single-query diagnostic calls now pass
their saved 28/30 values. The affected initial-HSP and gapped comparisons
remain exact after the correction. This is a function-input-state correction,
not a new NCBI output exception.

## Reproduce and current verification

From the repository root, with fresh output directories:

    python3 docs/evidence/tlosan_stage_c/run_chunk_purge_cap_trace.py /tmp/tlosan-c-purge-cap-new
    python3 docs/evidence/tlosan_stage_c/run_long_chunk_trace.py /tmp/tlosan-c-long-multi-new --multi-query
    python3 docs/evidence/tlosan_stage_c/run_word_xdrop_trace.py /tmp/tlosan-c-word-xdrop-new
    (cd docs/evidence/tlosan_stage_c/chunk_purge_cap_real_path_20260924 && sha256sum -c retained.sha256)
    (cd docs/evidence/tlosan_stage_c/long_multi_query_20260924 && sha256sum -c retained.sha256)
    (cd docs/evidence/tlosan_stage_c/word_xdrop_real_path_20260924 && sha256sum -c retained.sha256)
    cargo test --manifest-path LOSAT/Cargo.toml --lib algorithm::tblastn

The long runner generates the 15-MB subject in the fresh output directory;
only compact trace/output files and its exact subject checksum are retained.
The chunk and long compact streams reproduced byte for byte in fresh NCBI
runs before the WordFinder context-probe addition; the regenerated long
stream preserved all earlier candidate, gapped, chunk, output, and context
files byte for byte. All three retained SHA-256 manifests and fresh byte
replays passed. The independent auditor identified a missing four-INIT-row
check and a distinct-x-drop coverage limit. Both were addressed with exact
NCBI traces. The combined Stage C Rust suite passed 54/54, and `cargo clippy
--lib --tests -- -D warnings` passed. Production, test, and evidence diffs
were reviewed at this checkpoint. The independent auditor found no concrete
new fixture mismatch; Stage C acceptance remains open.

## Still required before Stage C acceptance

The 54-test fixture sweep compares saved function states, but the current
Rust diagnostics still call the full WordFinder pass before the separate
gapped/chunk-merge pass. Pinned `blast_engine.c:478-586,804-850` performs
WordFinder → GetGappedScore → chunk endpoint purge → offset adjustment →
merge for each chunk, then appends each frame. Integrate that execution order
in one internal Stage C path and compare its ordered event stream with NCBI
before accepting C. The positive chunk purge and append-cap paths above were
induced by comparison-only input changes; they establish those function
boundaries, not a natural high-count search. Source-backed coverage remains
bounded by the saved matrices, codes, subjects, query contexts, and one native
thread. Obtain independent read-only review of any final gate claim.
Do not start D's statistics/linking/composition/Kappa mapping or E's output
comparison until C is accepted. Code 32 remains a comparison-only NCBI C++
API oracle; database statistics and headers are not local `-subject` expected
output.
