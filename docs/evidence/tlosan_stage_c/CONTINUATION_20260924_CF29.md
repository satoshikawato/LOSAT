# TLOSAN Stage C continuation from cf29e853 — partial translation and merged HSP lists

Starting LOSAT commit: `cf29e853457e8226c0144a841f31b7dc75f4f308`, branch
`feature/tlosan-tblastn-v0.2.0`, clean worktree. Pinned NCBI source:
`598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4`; comparison-only TBLASTN
binary SHA-256:
`e3956f1e107a30439d56c8f72fae4267a7d62ebcc16e4f99d4baf2bdf00402e0`.
This follows [the earlier 2026-09-24 record](CONTINUATION_20260924.md).
The public TBLASTN CLI still returns `TBLASTN local search is unimplemented
(Stages C-E)`. This record is **not Stage C acceptance**, and no Stage D, Stage E,
completion, performance, or certification claim is made.

## Current comparison and first differences

The fresh fixtures use one protein query, one local nucleotide `-subject`,
BLOSUM62/word 3, threshold 13, window 40, gaps 11/1, E-value 10000,
SEG off, composition 0, sum statistics off, code 1, one thread. The
[runner](run_multi_hsp_trace.py) derives both inputs from the pinned six-frame
fixture, checks the NCBI source and binary hashes, runs NCBI with and without a
comparison-only [probe](ncbi_gapped_trace.c), and requires equal final bytes.
Each evidence directory contains input/output SHA-256, exact commands, raw
trace, and structured gapped, target-translation, traceback, and endpoint rows.

| Fixture | NCBI first boundary | Rust comparison after this change |
| --- | --- | --- |
| [Two copies, 2,364 nt](multi_hsp_20260924/) | Three gapped HSPs. First partial `+3` target window AA `113..298` hits a fence; NCBI rolls the entire three-HSP list back and retries on the full translated subject. | All three ordered gapped score/offset rows, translation ranges, and ordered final traceback score/offset rows exact. |
| [No-fence partial window, 2,364 nt](multi_hsp_no_fence_20260924/) | Same two-copy input, with both X-drop options set to 5 bits (raw 12/12). NCBI uses three partial translation windows (`+3` 113..298 and 467..652; `−1` 60..133), completes traceback with `fence_hit=0`, and returns three HSPs. | Window starts/stops and all three ordered score/offset rows exact; Rust confirms no full-subject retry. |
| [Six frames merged, 6,377 nt](six_frame_merged_20260924/) | Six expected raw-656 HSPs in one subject, plus lower HSPs. NCBI `GetGappedScore` emits 18 HSPs; traceback input is score ordered. A `−1` HSP grows from raw 38 to raw 236, query `0..82`, subject `842..990`; the endpoint purge deletes it because it shares query/subject start with a raw-656 `−1` HSP. The list goes 18→17. | All 18 per-frame ordered gapped rows, raw-236 pre-purge HSP, and all 17 ordered post-purge traceback rows exact. |

The first Rust divergence on the merged fixture was an extra raw-32 `−3`
HSP and frame-local order differences before traceback. The pinned NCBI
`aa_ungapped.c:234`, `blast_extend.c:273-313`, and
`blast_gapalign.c:3827-3834,3908-3919,4057-4091` show that WordFinder sorts
initial HSPs by score before `GetGappedScore`, which uses a new interval tree
for each frame and checks containment before gapped extension. Rust now uses
those same boundaries and reuses the existing NCBI-derived interval tree.

The pinned NCBI `blast_engine.c:522-552,804-850` then purges common endpoints,
sorts each frame list by score, and appends into the six-frame list. During
traceback, `blast_traceback.c:303-312,408-440,503-536,635-693,709-721,1644-1684`
gets a partial target translation per HSP, rolls back the original list on
fence, retries with full subject, purges endpoints, and sorts again.
`blast_hits.c:1154-1228,2268-2379,2455-2537` defines translation windows and
the two endpoint comparators. The Rust internal score/offset diagnostic now
reproduces the observed window, successful nonzero-base partial traceback, and list retry with translation-window-sized
storage (the late subject offset does not allocate an unused prefix). It applies the TBLASTN
`purge=TRUE` common-endpoint ordering after traceback, and follows the pinned `blast_traceback.c:352-360,401-405,583-605,675-693` interval-tree ordering before and after that purge. The source snippets
and line numbers are directly above the Rust ports.

For this default non-greedy TBLASTN traceback, NCBI
`blast_hits.c:2455-2537` forces endpoint `purge=TRUE` for programs other than
BLASTN. `blast_traceback.c:635-666` then receives `extra_start == hspcnt`, so
the subsequent `Blast_HSPReevaluateWithAmbiguitiesGapped` loop has no members.
The comparison probe records zero calls on all three new fixtures. This explains
these runs; it does **not** establish other modes or the separate Kappa redo
path. The new `−1` HSP deletion is at the endpoint purge, after its traceback
score change, not at the ambiguity reevaluation loop.

## Verification

```bash
python3 docs/evidence/tlosan_stage_c/run_multi_hsp_trace.py /tmp/tlosan-c-two-new
python3 docs/evidence/tlosan_stage_c/run_multi_hsp_trace.py /tmp/tlosan-c-no-fence-new --no-fence
python3 docs/evidence/tlosan_stage_c/run_multi_hsp_trace.py /tmp/tlosan-c-six-new --six-frame
(cd docs/evidence/tlosan_stage_c/multi_hsp_20260924 && sha256sum -c outputs.sha256)
(cd docs/evidence/tlosan_stage_c/multi_hsp_no_fence_20260924 && sha256sum -c outputs.sha256)
(cd docs/evidence/tlosan_stage_c/six_frame_merged_20260924 && sha256sum -c outputs.sha256)
cargo fmt --manifest-path LOSAT/Cargo.toml --check
cargo test --manifest-path LOSAT/Cargo.toml --lib algorithm::tblastn::search_gapped::tests
```

The retained checksums passed. The focused eight Rust oracle tests passed.
`LOSAT/target/debug/LOSAT tblastn -query <fixture/query.faa> -subject
<fixture/subjects.fna> -outfmt 6` exited 1 with
`Error: TBLASTN local search is unimplemented (Stages C-E)`.
These tests exercise the stated BLOSUM62/word-3 score/offset boundary and do
not certify the full TBLASTN pipeline. No formal benchmark was run because
Stage C/D/E parity gates remain open.

## Remaining gates, in order

1. **Stage C remains open.** Complete a direct NCBI comparison of every
   candidate and HSP on multiple protein queries, alternate matrix/word-size
   paths, and subjects above NCBI's 5,000,000 translated-residue chunk limit
   (`blast_gapalign.h:54`, `blast_engine.c:221-325`). Complete the general
   `Blast_HSPListAppend` and its cap, preliminary endpoint purge, positive
   traceback containment and `Blast_HSPTest` deletion cases, and query-context
   state across those profiles. NCBI `blast_traceback.c:436-445,585-605` frees
   HSPs on start-offset failure or `Blast_HSPTest`; the current Rust internal
   boundary does not yet prove these deletion paths. The current
   Rust entry is internal and limited to one query, BLOSUM62/word 3, and the
   measured traceback fields; it is not the public search engine.
2. **Stage D has not begun.** Only after all Stage C cases match, create the
   function/input-state/execution-order table and implement effective lengths,
   sum statistics, default intron/linking, composition mode 2, Kappa redo,
   E-value/bit score, filtering, rank, ties, and deletion order. Code 32 uses
   the comparison-only NCBI C++ API oracle; `-db` statistics and headers are
   not a local `-subject` oracle.
3. **Stage E is gated by C and D.** Only then byte-compare outfmt 0/6/7 against
   the local `-subject` oracle and retain the explicit public CLI
   unimplemented error until all three formats are complete.

Program/task: TBLASTN local `-subject`. Target: native Linux, one thread.
LOSAT owner: `LOSAT/src/algorithm/tblastn/search_gapped.rs`; existing shared
interval tree: `LOSAT/src/algorithm/blastn/interval_tree.rs`. Accepted
exception: none on these code-1 runs. Native/Wasm, release certification, and
performance: not applicable because later gates are open.
