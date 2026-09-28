# TLOSAN Stage C continuation after 2e607cca

Branch: `feature/tlosan-tblastn-v0.2.0`; starting commit:
`2e607cca123dd678107022b9241b4f08f158ba53` (the checkpoint after
`97ba376499a6c5a28b72c6f45d80a00fb4f002fd`). Pinned NCBI C/C++
source: `598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4`; comparison-only
`tblastn` binary SHA-256:
`e3956f1e107a30439d56c8f72fae4267a7d62ebcc16e4f99d4baf2bdf00402e0`.
Read [the prior continuation](CONTINUATION_20260924_97BA.md) for the three
starting fixtures, including the original long-subject difference.

**Stage C remains open. Stages D and E have not begun. This is neither
TBLASTN completion nor a parity, performance, or release certification.**
The public TBLASTN CLI still reports the explicit Stage C-E unimplemented
error.

## Reproduced first difference and NCBI owner

The prior long-subject trace found the first difference at the
5,000,000-residue split before WordFinder. The Rust diagnostic now follows
`blast_engine.c:221-325,430-586,747-850`: six translated frames, strict
`offset + 5,000,000 < right` split, 100-residue overlap, per-chunk scan and
WordFinder, per-chunk gapped call and endpoint purge, offset adjustment,
`Blast_HSPListsMerge`, then frame-level `Blast_HSPListAppend`. The saved
[unmasked 15-MB fixture](long_chunk_20260924/) has two chunk-local HSPs
(raw 279 and 656), one merged raw-660 HSP, and six append snapshots. Rust
matches these recorded rows. This is an internal diagnostic path; the
public CLI has not been connected to it.

A new direct [NCBI input-state probe](ncbi_chunk_ranges_trace.c) then
exposed the next difference in the [masked chunk-boundary fixture](masked_chunk_boundary_20260924/).
At WordFinder call 0, NCBI supplied three `seq_ranges`:
`(0,4999900)`, `(4999999,4999900)`, `(4999999,5000000)`.
Rust initially supplied only the first and third. The saved probe records
`SET_RANGES=3`, all initialized setup rows after `BlastSeqBlkSetSeqRanges`,
and every later `RANGE_CALL`/`RANGE` row. The first differing NCBI/Rust
call state is reproducible by the fixture runner and the focused Rust
assertion; the initial failed assertion is described here, while the
retained oracle is the NCBI trace.

The extra inverted range follows pinned source rather than a heuristic:
`blast_fasta_input.cpp:489-502` requests masks on both strands;
`blast_aux_priv.cpp:436-454` pushes plus and minus copies of each lowercase
interval; `blast_setup_cxx.cpp:689-707,813-826` packs both;
`blast_util.c:211-215` fills the first/last complement boundaries;
`blast_engine.c:815-829,259-310` converts and slices those ranges for each
frame/chunk. Rust now preserves both copies and their inverted complement.
Its masked-boundary `seq_ranges` match all 12 WordFinder calls exactly,
including zero/negative-width entries. The earlier final-HSP comparison
had missed this input-state difference.

## Function, input state, and execution-order comparisons

| Pinned NCBI function and state | Saved NCBI oracle | Rust diagnostic result |
| --- | --- | --- |
| `s_GetNextSubjectChunk`, `blast_engine.c:221-325`: frame length, soft ranges, chunk offset | Unmasked 5,000,000 + 420; masked boundary 12 calls with exact ranges; fully masked middle omitted from positive frames | All recorded call lengths and ranges exact for these fixtures |
| `BlastAaWordFinder`, `aa_ungapped.c:200-234,492-614`: shared query lookup, diagonal state, chunk-local translated sequence | Unmasked boundary 2 initial HSPs; masked boundary one retained HSP; 30-MB middle-mask 15 calls; two-sided fixture first and later positive chunks each one raw-656 HSP | Recorded initial HSP fields and order exact; long-subject full candidate stream is still untraced |
| `BLAST_GetGappedScore` and endpoint purge, `blast_engine.c:478-552`: chunk-local HSP list and tree | Unmasked boundary raw 279/656; later chunk of middle-mask raw 656; two-sided fixture raw 656 in both called positive chunks | Recorded local score, frame, offsets, and gapped starts exact; these cases do not trigger positive endpoint deletion |
| `Blast_HSPListsMerge`, `blast_engine.c:572-586`, `blast_hits.c:2857-3035`: adjusted chunk offsets and overlap | Unmasked overlap pair becomes raw 660; two-sided nonoverlapping equal-score HSPs remain in first/later order | Recorded merge input/output fields exact in the tested cases |
| `Blast_HSPListAppend`, `blast_engine.c:842-850`, `blast_hits.c:213-232`: per-frame list and cap | Six 7-field snapshots; two-sided fixture has two HSPs in every snapshot | All six recorded snapshot rows exact at observed `INT4_MAX` cap; separate direct NCBI C API cap-3 oracle remains in the prior checkpoint |

The [fully masked middle fixture](no_range_middle_20260924/) puts the
lowercase interval at nucleotide `[14997000,30000000)`. For positive
frames, NCBI calls WordFinder on the first and last chunks and returns
`SUBJECT_SPLIT_NO_RANGE` on the middle chunk. Rust matches all 15
`RANGE_CALL` states and the later raw-656 HSP at local subject 250..370.
The [two-sided fixture](no_range_two_hits_20260924/) adds an unmasked
copy at translated offset 4,998,450, before that masked span, while
retaining the later copy at 10,000,050. NCBI returns raw-656 HSPs in
both chunks, merges them in first/later order, and appends two rows in
all six frame snapshots. The Rust test compares initial HSPs, local
gapped HSPs, and all append rows, including equal-score order. The
30-MB subject sequences are generated only in `/tmp`; compact retained
traces, manifests, and SHA-256 files are in these directories.

Every comparison probe is external to LOSAT's runtime/build and preserves
the unprobed NCBI final output bytes for its fixture. `SET_RANGE` rows are
captured after NCBI initializes boundaries, so the retained values do not
include uninitialized buffer contents.

## Remaining Stage C work

1. Compare every candidate and WordFinder/gapped row for the saved
   [BLOSUM45/word-2 fixture](alternate_matrix_word2_20260924/). The Rust
   lookup, extension, and gapped entry points remain BLOSUM62/word 3 only;
   traceback replay of NCBI preliminary HSPs does not close this gap.
2. Capture a complete-path NCBI fixture where start-offset acquisition
   fails and one where `Blast_HSPTest` deletes an HSP. The direct API
   HSPTest boundary oracle and positive containment fixtures do not
   establish those full-path timing/deletion cases.
3. Find and compare a positive chunk-local endpoint-purge deletion, a
   chunk-path `Blast_HSPListAppend` cap boundary, full query-context
   partition/rank, and long-subject ordered candidate streams. The
   current fixtures cover the exact observed rows, not every possible
   input state.
4. Sweep all Stage C fixture differences together before acceptance;
   source-specific changes here establish only the bounded profiles
   above. Do not advance to D until every Stage C function/input-state/
   execution-order comparison is exact.

Only after Stage C is exact, construct Stage D's function, input-state,
and order table for effective lengths/statistics, linking, composition
adjustment, Kappa redo, rank, and deletion. Only after C and D agree,
compare local `-subject` outfmt 0/6/7 bytes in E. Code 32 uses the
comparison-only NCBI C++ API oracle; `-db` statistics and headers are not
local-subject expected output.

## Reproduce

From repository root, use fresh output directories:

    python3 docs/evidence/tlosan_stage_c/run_long_chunk_trace.py /tmp/tlosan-c-long-new
    python3 docs/evidence/tlosan_stage_c/run_long_chunk_trace.py /tmp/tlosan-c-mask-boundary-new --mask-boundary
    python3 docs/evidence/tlosan_stage_c/run_no_range_middle_trace.py /tmp/tlosan-c-no-range-new
    python3 docs/evidence/tlosan_stage_c/run_no_range_middle_trace.py /tmp/tlosan-c-two-hit-new --early-hit

The scripts require the pinned source commit and binary hash, build
comparison-only probes in temporary directories, and reject any probe
that changes unprobed output bytes. The retained fixture manifest gives
exact inputs and options. The source sequences are generated from saved
small inputs so no 15- or 30-MB FASTA is committed.

## Local verification

- `cargo test --lib algorithm::tblastn::search_seed::tests`:
  11 passed, 0 failed.
- `cargo test --lib algorithm::tblastn::search_init::tests`:
  8 passed, 0 failed, including exact masked-boundary and no-range
  `seq_ranges` rows against the direct NCBI probe.
- `cargo test --lib algorithm::tblastn::search_gapped::tests`:
  18 passed, 0 failed, including two nonoverlapping equal-score HSPs
  across the skipped middle chunk and six exact append snapshots.
- `cargo fmt --check`, `cargo build --release`, and
  `cargo clippy --lib -- -D warnings`: passed. Both modified fixture
  runners passed `python3 -m py_compile`.
- Rebuilt public `LOSAT tblastn -query ... -subject ... -db_gencode 32`:
  exit 1, zero stdout bytes, explicit
  `Error: TBLASTN local search is unimplemented (Stages C-E)`.
- The three new compact fixture `retained.sha256` files passed all
  entries; `git diff --cached --check` passed. No generated 15-/30-MB
  subject or NCBI library/binary is staged.

## Independent read-only audit

The `ncbi_parity_auditor` rechecked the pinned source and binary hash,
the duplicate-mask source chain, WordFinder input ranges, two-sided
chunk ordering, six append snapshots, and all three retained checksum
files. It found no concrete new NCBI/Rust mismatch in these bounded
one-mask fixtures. It confirmed that the per-interval plus/minus copy
order follows `blast_aux_priv.cpp:436-454` and
`seqlocinfo.cpp:126-139`. The exact evidence limits remain the open
Stage C items above: long ordered candidates, positive endpoint deletion,
chunk-path cap, and complete-path deletion cases are unverified.
Multiple separated lowercase intervals are source-supported but have no
independent fixture comparison in this checkpoint. The auditor also
identified the prior continuation's stale whole-frame statement; it is
now marked historical and superseded there.
