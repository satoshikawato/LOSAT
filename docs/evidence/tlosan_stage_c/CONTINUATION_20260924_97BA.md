# TLOSAN Stage C continuation from 97ba376

Starting branch: feature/tlosan-tblastn-v0.2.0; commit:
97ba376499a6c5a28b72c6f45d80a00fb4f002fd; clean worktree.
Pinned NCBI source: 598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4.
Comparison tblastn binary SHA-256:
e3956f1e107a30439d56c8f72fae4267a7d62ebcc16e4f99d4baf2bdf00402e0.
This follows [the cf29e853 checkpoint](CONTINUATION_20260924_CF29.md).
The public TBLASTN CLI still exits with the explicit Stage C-E unimplemented error.

**Stage C remains open. D and E have not begun. This is not TBLASTN
completion, output parity, performance certification, or release evidence.**

## First difference, NCBI owner, and measured order

[The three-query fixture](multi_query_20260924/) uses the prior 6,377-nt
six-frame subject, a full 120-aa query, its residues 25..95 as a second
70-aa query, and a 120-aa no-hit X query. Options: local -subject, code 1,
BLOSUM62/word 3, threshold 13, window 40, gap 11/1, E-value 10000,
SEG off, composition 0, sum statistics off, one thread. The
[runner](run_multi_query_trace.py) verifies the source/binary pins,
records commands and hashes, and requires each of three comparison probes
to preserve the unprobed NCBI final output bytes.

The first Rust difference was the eighth ordered WordFinder HSP (index 7):
Rust had frame +3, raw 13, query/subject starts 50/1011; NCBI had frame
+3, raw 18, starts 31/1150. Both had 32 HSPs. The pre-fix Rust
row was observed during this diagnostic run but was not retained as a
separate raw trace; the saved NCBI trace reproduces its reference row.
NCBI
aa_ungapped.c:200-234 calls Blast_InitHitListSortByScore after each frame's
WordFinder call; blast_extend.c:273-313 defines its comparator. The Rust
diagnostic now sorts at this per-frame boundary. Its shared query lookup,
context starts, and global-to-local HSP conversion follow
blast_query_info.c:68-96, aa_ungapped.c:547-583, and
blast_gapalign.c:2371-2468.

| NCBI function and input state | NCBI trace | Rust internal comparison |
| --- | --- | --- |
| BlastQueryInfoNew, protein contexts | offsets 0/121/192; lengths 120/70/120; validity 1/1/0 | all contexts exact |
| BlastAaWordFinder, shared lookup and diagonal table | 2,832 ordered candidates; 32 ordered initial HSPs | all candidate and HSP rows exact |
| BLAST_GetGappedScore, context-local query and frame-local tree | 31 ordered HSPs with context, score, frame, offsets, gapped starts; positive preliminary containment in frame -3 | all 31 rows exact, including omission of the contained candidate |
| Endpoint purge and Blast_HSPListAppend | cumulative counts 3/5/11/19/25/31 | all six ordered seven-field append snapshots exact (context, score, frame, q/s starts and ends); gapped starts and edit scripts are not in this trace |
| Blast_TracebackFromHSPList, query-indexed lists | context 1: 13 input, fence retry, 11 final; context 0: 18 input, fence retry, 17 final | both final lists' score, frame, offsets, order exact |

The NCBI traceback call order here is context 1 then context 0. Rust
compares the already separated NCBI lists in that order. This does not
establish the public HSP-stream partition or later display order.
The independent read-only parity audit found that the preliminary save
condition still used score >= 0. Pinned blast_gapalign.c:3924-3927,4058
uses hit_params->cutoffs[context].cutoff_score. The internal Rust
comparison now takes a measured cutoff for each query context. The
[the independent per-context probe](context_cutoffs_20260924/)
measures 0/0/INT4_MAX for the three-query BLOSUM62 and BLOSUM45
profiles: the no-hit X context is invalid. All six frame calls agree,
and the Rust tests read those measured values. The saved code-32 C++
API fixture supplies 1; the long-subject trace reports 28. The long-subject search
path is still unported, so the cutoff fix does not validate that path.

[The direct comparison-only NCBI C API cap oracle](append_cap_20260924/)
calls Blast_HSPListAppend with cap 3 on two three-HSP lists. It retains
scores 100, 90, then the old-list 80 when old/new 80 tie on every
ScoreCompareHSPs field. The Rust preliminary append diagnostic agrees on
context, score, frame, offsets, and order. The C harness is never in
LOSAT's runtime or build.

## Alternate matrix, word size, and traceback deletion

[The BLOSUM45/word-2 fixture](alternate_matrix_word2_20260924/)
uses the same three query contexts against the saved 362-nt +1 subject,
with threshold 16, window 60, gap 14/2, and E-value 10000. Its pinned
NCBI probes preserve the unprobed output bytes. NCBI records 558
ordered candidate rows, 27 WordFinder initial HSPs, and 22 preliminary
gapped HSPs. Rust's lookup, two-hit WordFinder, and preliminary gapped
entry points still explicitly cover BLOSUM62/word 3 only. They have
not been compared to this alternate profile.

The selected-matrix Rust traceback diagnostic replays NCBI's 22
preliminary HSPs, including their gapped starts, rather than claiming
Rust produced them. NCBI's traceback trace has positive containment;
Rust's final per-query ordered score/frame/coordinate rows match the
3 and 7 NCBI final rows. This verifies a positive traceback containment
case for those inputs, not the general deletion boundary.

[A direct pinned NCBI C API oracle](hsp_test_20260924/) exercises
Blast_HSPTest at identity-threshold equality, just above equality,
and minimum-length boundaries. Rust's filter matches all six direct
API results and now runs after its traceback alignment, before adding
a retained HSP to the interval tree. All retained search fixtures use
the default 0%/0-length filter and their traced HSP_TEST calls return
zero, so a positive deletion inside the complete NCBI traceback path
is still missing.

## Long translated subject

[The compact long-chunk evidence](long_chunk_20260924/) retains query,
NCBI output and traces, command, generated-subject SHA-256, and recipe.
[Its runner](run_long_chunk_trace.py) creates 15,000,962 nt from ATG x
4,999,950 codons, the saved 362-nt +1 insert, and ATG x 250 codons.
The first frame has 5,000,320 translated residues. The generated 15-MB
subject remains outside Git; the manifest pins its checksum. Both probes
left NCBI output bytes unchanged.

NCBI blast_engine.c:221-325,430-586 makes 5,000,000 and 420-residue
WordFinder chunks in every frame (12 calls). The +1 first chunk has
raw-279 HSP at subject 4,999,950..5,000,000. The overlapping second
has raw-656 HSP at chunk-local 50..170. NCBI adjusts the latter to
4,999,950..5,000,070, then calls Blast_HSPListsMerge with split
4,999,900 and overlap 100. Pinned blast_hits.c:1450-1537,2857-3035
merges the two into one raw-660 preliminary HSP at
4,999,950..5,000,070. Rust's isolated merge diagnostic matches both
captured NCBI input/output transitions, including score and gapped starts.
Later full-subject NCBI traceback returns raw 656.

Rust preliminary search still scans whole translated frames. It does
not yet produce both chunk-local WordFinder/GetGappedScore inputs.
**The first remaining long-subject difference is the split before
WordFinder**, while this pairwise overlap-merge transition matches.

## Remaining Stage C gates, in order

1. Port and compare the general translated-subject chunk lifecycle at
   blast_engine.c:221-325,430-586, including every frame, mask range,
   chunk offset, overlap list, endpoint purge, and cap. The pairwise
   merge test does not prove the entire path.
2. Port and compare the saved BLOSUM45/word-2 candidate, WordFinder,
   and preliminary gapped inputs and rows. Traceback replay alone does
   not close this profile.
3. Exercise start-offset acquisition failure and positive Blast_HSPTest
   deletion within a complete NCBI traceback call, with exact inputs
   and deletion order. Positive *preliminary* and *traceback*
   containment are captured, and the direct HSPTest API agrees, but
   those other complete-path deletions have not been observed.
4. Compare general preliminary endpoint-purge deletion and HSP-stream
   query-context partition/rank over varied inputs. Current tests prove
   only their observed fields and lists.

Only after Stage C matches fully, create the Stage D
function/input-state/execution-order table and compare effective lengths,
statistics, linking, composition mode 2, Kappa redo, rank, and deletion.
Only after C and D agree, compare local -subject outfmt 0/6/7 bytes in E.
Code 32 remains a comparison-only NCBI C++ API oracle case; database
statistics/headers cannot serve as the local-subject oracle.

## Reproduce and verification

Run each fixture into a fresh path from repository root:

    python3 docs/evidence/tlosan_stage_c/run_multi_query_trace.py /tmp/tlosan-c-multi-query-new
    python3 docs/evidence/tlosan_stage_c/run_long_chunk_trace.py /tmp/tlosan-c-long-chunk-new
    python3 docs/evidence/tlosan_stage_c/run_multi_query_trace.py /tmp/tlosan-c-alternate-new --blosum45-word2
    python3 docs/evidence/tlosan_stage_c/run_context_cutoff_trace.py /tmp/tlosan-c-context-cutoffs-new
    gcc -std=c11 -O2 -Wall -Wextra -o /tmp/tlosan-c-append-cap-oracle docs/evidence/tlosan_stage_c/ncbi_append_cap_oracle.c -ldl
    /tmp/tlosan-c-append-cap-oracle /home/kawato/micromamba/lib/ncbi-blast+/libxblast.so
    (cd docs/evidence/tlosan_stage_c/multi_query_20260924 && sha256sum -c outputs.sha256)
    (cd docs/evidence/tlosan_stage_c/long_chunk_20260924 && sha256sum -c retained.sha256)
    (cd docs/evidence/tlosan_stage_c/append_cap_20260924 && sha256sum -c retained.sha256)
    (cd docs/evidence/tlosan_stage_c/alternate_matrix_word2_20260924 && sha256sum -c outputs.sha256)
    gcc -std=c11 -O2 -Wall -Wextra -o /tmp/tlosan-c-hsp-test-oracle docs/evidence/tlosan_stage_c/ncbi_hsp_test_oracle.c -ldl
    /tmp/tlosan-c-hsp-test-oracle /home/kawato/micromamba/lib/ncbi-blast+/libxblast.so
    (cd docs/evidence/tlosan_stage_c/hsp_test_20260924 && sha256sum -c retained.sha256)
    (cd docs/evidence/tlosan_stage_c/context_cutoffs_20260924 && sha256sum -c retained.sha256)
    cargo fmt --manifest-path LOSAT/Cargo.toml --check
    cargo test --manifest-path LOSAT/Cargo.toml --lib algorithm::tblastn::search_seed::tests::multi_query_candidate_contexts_and_order_match_ncbi -- --exact
    cargo test --manifest-path LOSAT/Cargo.toml --lib algorithm::tblastn::search_init::tests::multi_query_wordfinder_hsps_match_ncbi_in_order -- --exact
    cargo test --manifest-path LOSAT/Cargo.toml --lib algorithm::tblastn::search_gapped::tests

Program/task: TBLASTN local -subject code 1; cap test: direct comparison-only
NCBI C API. LOSAT owners: search_seed.rs, search_init.rs, search_gapped.rs.
NCBI functions/lines: above. First observed difference: initial HSP order
at index 7; first current reproducible unported difference: subject
chunk split before WordFinder on the long-subject fixture. The saved
NCBI FRAME_CHUNK trace records two calls per frame while the Rust
search_seed/search_init source still iterates one whole translated frame. Accepted exception: none on these code-1 runs. Broader
output-byte gate, native/Wasm, benchmark, Stage D/E, release audit: not
applicable while C remains open.

Focused checks completed for this checkpoint:
- Retained multi-query, compact long-chunk, append-cap, alternate-profile,
  HSPTest direct API, and context-cutoff SHA-256 checks: all passed.
- cargo fmt --check passed. Source/document diff whitespace checks passed.
  The six retained raw NCBI stderr files include one trailing space in
  NCBI's invalid-query warning each; plain git diff --check reports
  those unmodified oracle bytes.
- TBLASTN seed 10/10, WordFinder 5/5, gapped/traceback/append/merge/
  HSPTest 15/15 focused Rust tests: passed.
- cargo build (debug): passed.
- Newly built public CLI with the multi-query fixture exited 1 with
  "Error: TBLASTN local search is unimplemented (Stages C-E)".
- No formal benchmark was run because Stage C remains open.
