# LOSATP hot path performance priorities

Date: 2026-05-11

Scope: LOSATP means the Rust `blastp` implementation under
`LOSAT/src/algorithm/blastp/`.  The target is NCBI BLASTP speed parity while
preserving bit-perfect output parity.  NCBI BLAST remains reference-only; none
of the items below require or permit invoking NCBI from LOSAT runtime code.

## Reference map

The Rust pipeline currently mirrors these NCBI stages:

| Stage | LOSAT code | NCBI reference |
|---|---|---|
| Subject scan and offset pair emission | `LOSAT/src/algorithm/blastp/blast_engine.rs:1506`, `LOSAT/src/algorithm/tblastx/blast_aascan.rs:184` | `c++/src/algo/blast/core/blast_aascan.c:48`, `:234`; `aa_ungapped.c:440`, `:713` |
| One-hit/two-hit ungapped extension | `LOSAT/src/algorithm/blastp/blast_engine.rs:3688`, `LOSAT/src/algorithm/blastp/extension.rs:151` | `c++/src/algo/blast/core/aa_ungapped.c:831`, `:886`, `:1020`, `:1089` |
| Initial HSP score sort, chaining, restricted alignment setup | `LOSAT/src/algorithm/blastp/blast_engine.rs:3972`, `:3981`, `:4015` | `c++/src/algo/blast/core/blast_gapalign.c:3517`, `:3681`, `:3736` |
| Preliminary score-only gapped extension | `LOSAT/src/algorithm/blastp/blast_engine.rs:4065`, `LOSAT/src/algorithm/blastp/gapalign.rs:1055` | `c++/src/algo/blast/core/blast_gapalign.c:736`, `:4209` |
| Composition-based redo and traceback | `LOSAT/src/algorithm/blastp/blast_engine.rs:4385`, `LOSAT/src/algorithm/blastp/kappa.rs:1269` | `c++/src/algo/blast/core/blast_kappa.c:3009`, `:3450`, `:3658`; `composition_adjustment/redo_alignment.c:1101` |
| HSP stream ordering and final collection | `LOSAT/src/algorithm/blastp/hsp.rs:991`, `LOSAT/src/algorithm/blastp/blast_engine.rs:4653` | `c++/src/algo/blast/core/blast_hspstream.c:271`, `:339`; `blast_kappa.c:2500`, `:3817` |

## Current timing snapshot

Diagnostic command, run from repository root:

```bash
LOSAT_TIMING=1 LOSAT/target/release/LOSAT blastp \
  -query LOSAT/tests/fasta/AP027131.faa \
  -subject LOSAT/tests/fasta/AP027133.faa \
  -outfmt 6 -evalue 1e-5 -max_target_seqs 10 \
  -num_threads 1 \
  -out /tmp/losatp_profile_ap027131_ap027133_n1.out
```

| Mode | Threads | Total | Main observed costs |
|---|---:|---:|---|
| default `blastp` | 1 | 5.033 s | scan/ungapped 1.733 s, preliminary gapped 1.997 s, kappa redo 1.110 s |
| default `blastp` | 8 | 1.914 s | parallel scan/gapped aggregate 4.816 s, serial kappa redo 1.068 s |
| `blastp-fast` | 8 | 1.261 s | compressed lookup build 0.954 s, kappa redo 0.264 s |

The per-stage scan/gapped numbers under `-num_threads 8` are aggregate atomic
durations across workers, so wall time is the `blastp_total` line.  The key
finding is still clear: default BLASTP spends most work in scan/ungapped plus
preliminary gapped extension, while multi-threaded default BLASTP is then
limited by serial composition-based redo.  `blastp-fast` is dominated by
compressed lookup construction.

## Priority order

### P0. Parallelize composition-based redo

Expected effect: highest for `-num_threads > 1`; in the measured 8-thread
default run, serial kappa redo was about 56% of wall time.

NCBI does this work with thread-local query info, score blocks, gap-align
scratch, composition workspaces, and redone-match heaps, then consolidates
thread-local results.  See `blast_kappa.c:3009-3335`, `:3450-3678`, and
`:3817-3850`; heap insertion/pop order is in `compo_heap.c:330-465`.

LOSAT currently enters kappa after preliminary search and iterates query/local
matches serially in `blast_engine.rs:4385-4626`, using one
`BlastCompositionWorkspace`, one `GapAlignScratch`, and one `redone_matches`
vector.  The direct porting target is:

- partition preliminary `BlastpHitList` work into deterministic per-thread
  batches after `preliminary_hit_lists` have been sorted by NCBI e-value order;
- allocate thread-local `BlastCompositionWorkspace`, `GapAlignScratch`,
  `kappa_preliminary_hits`, and `Vec<BlastCompoHeap>`;
- run `postprocess_preliminary_hits` in parallel without sharing heap state;
- consolidate thread-local hit lists with NCBI heap/e-value comparators, matching
  `SThreadLocalDataArrayConsolidateResults` behavior.

Parity risk: high.  The redone heap tie-breakers, early termination, and final
read order must match NCBI exactly.  Add counters before changing behavior:
local matches redone, early-terminated matches, per-query heap insertions,
discarded heap records, and final per-query subject order.

### P0. Specialize preliminary gapped DP hot loops

Expected effect: high for default `blastp`; preliminary score-only gapped
extension was the largest single-thread cost in the measured default run.

NCBI `Blast_SemiGappedAlign` keeps a rolling `BlastGapDP*` window and uses
matrix-row pointers in the inner loop (`blast_gapalign.c:736-957`).  LOSAT
matches the algorithm in `gapalign.rs:1055-1228`, but the Rust loop still pays
for repeated slice indexing, branchy helper calls, and `Vec` growth checks on
the hot path.

Recommended implementation sequence:

- add DP-cell counters for exact and restricted preliminary alignments;
- split a dedicated BLOSUM62/non-reverse score-only function for the dominant
  default path, with raw pointers for `dp_mem`, query, subject, and the matrix
  row after proving all ranges from the NCBI scan contract;
- keep the generic safe path for non-default matrices and tests until the
  specialized path has parity coverage;
- apply the same treatment to restricted gapped alignment only after counters
  show it is frequently selected.

Parity risk: medium-high.  The x-drop loop boundaries, `b_size`, `first_b_index`,
and `best_score` update order must remain byte-for-byte equivalent to
`blast_gapalign.c:736-957`.

### P0. Specialize scan plus ungapped extension for the default two-hit path

Expected effect: high for default `blastp`; scan/ungapped was 1.733 s of the
5.033 s single-thread run.

NCBI scans subject words into `BlastOffsetPair` batches, then immediately walks
those hits through the two-hit diagonal table and `s_BlastAaExtendTwoHit`
(`blast_aascan.c:48-131`; `aa_ungapped.c:440-588`, `:1089-1177`).  LOSAT mirrors
this in `blastp_for_each_offset_pair` and a closure in `blast_engine.rs:1506`
and `:3688-3958`.

Best first change is not to alter NCBI timing/order, but to create a
monomorphized default-path worker:

- `BLOSUM62`, word size 3, two-hit, standard lookup;
- inline scan result walking, diagonal update, context cutoff lookup, and
  `extend_two_hit_blosum62`;
- raw-pointer residue scoring in `extend_left`/`extend_right`, using the
  existing NCBI comments in `extension.rs` as the safety boundary.

Keep the existing generic path for one-hit, compressed lookup, non-BLOSUM62,
and tests.  This should remove closure dispatch, repeated enum matches, and
some slice bounds checks without changing the NCBI batch order.

Parity risk: medium.  The diagonal flag/reset order in `aa_ungapped.c:509-588`
is fragile and must be preserved exactly.

### P1. Optimize `blastp-fast` compressed lookup construction

Expected effect: very high for `blastp-fast`; measured 8-thread
`blastp-fast` spent 0.954 s of 1.261 s in lookup construction.

NCBI compressed lookup generation is in `blast_aalookup.c:1137-1242` and
scan-time layout is consumed by `blast_aascan.c:234-384`.  LOSAT constructs the
compressed table in `lookup/compressed.rs:1178`, recursively generates
neighbors in `:430-540`, and finalizes the PV in `:860`.

Likely wins, in order:

- flatten `CompressedNeighborInfo.matrix_sorted` and `matrix_sorted_char` from
  `Vec<Vec<_>>` into fixed `[BLASTAA_SIZE][15]`-style arrays;
- remove the `CompressedBackbonePayload` enum clone in `add_word_hit`
  (`lookup/compressed.rs:272`) by storing an NCBI-shaped inline/overflow struct
  or using a small tagged representation;
- reserve overflow cells from a pre-count or growth heuristic before recursive
  neighbor insertion;
- precompute compressed query-word windows once in
  `prepare_blosum62_lookup_query_for_word_size` and pass word indices into the
  recursive generator.

Parity risk: medium.  Neighbor append order is output-relevant because scan
emission order is subject-major then query-offset order.

### P1. Reduce per-redo range copying and SEG work

Expected effect: medium, and higher on long subjects with many CBS redo
alignments.

NCBI obtains a protein range in `blast_kappa.c:1573-1646` and then moves the
range pointer so the returned data has a fresh leading sentinel.  LOSAT
`build_subject_range_data` currently starts from the full subject
`BlastCompoSequenceData::from_ncbistdaa(matching_seq.data)` and then calls
`fit_protein_range_in_place` (`kappa.rs:701-785`).

Optimization target:

- add counters for bytes copied in query range and subject range construction;
- build the subject range directly with one leading sentinel and one trailing
  sentinel;
- apply SEG only to the requested subject range, but keep NCBI's
  `subject_maybe_biased` and near-identical timing from `blast_kappa.c:1626-1646`.

Parity risk: medium-high.  SEG interval coordinates and near-identical bypass
must stay aligned with NCBI.

### P1. Replace collect-sort preliminary merge with NCBI-style streaming merge

Expected effect: medium for large databases and high subject counts; low in the
measured 543-subject case.

NCBI writes each subject HSP list into an HSP stream under a lock
(`blast_engine.c:1409-1554`, `blast_hspstream.c:339-385`), then closes/reads
the stream in NCBI order (`blast_hspstream.c:271-327`).  LOSAT parallel search
currently collects all `BlastpSubjectResult`s, sorts by `subject_index`, and
then merges (`blast_engine.rs:4301-4334`).

This is parity-safe only if the resulting hit-list update order is identical.
The implementation should use per-worker buffers plus a deterministic reducer
rather than pushing directly into shared hit lists.

Parity risk: medium.  Main failure mode is subject-order drift on equal
e-values/scores.

### P2. Optimize common endpoint purge and interval-tree rebuilds

Expected effect: workload-dependent; likely medium only on high-hit datasets.

LOSAT uses NCBI-compatible sort/purge passes for preliminary endpoints
(`blast_engine.rs:2829-2886`) and rebuilds the preliminary interval tree during
restricted-alignment exact retry (`blast_engine.rs:4161-4174`).  NCBI performs
the same retry control flow in `blast_gapalign.c:3955-4007`.

Before changing this, add counters for:

- exact retry count;
- number of HSPs removed by query on retry;
- number of interval-tree rebuilds and inserted nodes;
- number of endpoint duplicates purged.

Only optimize if those counters are material.  Possible safe shape is a
per-query preliminary list index that lets the exact retry remove/rebuild the
current query without scanning unrelated HSPs.

Parity risk: high, because retry order and interval containment directly affect
which HSPs survive.

### P2. Trim final allocation and formatting overhead

Expected effect: low for default outfmt 6 in the measured run, but useful for
large output and pairwise output.

Current finalization clones/converts HSPs through `collect_hits_from_hit_lists`
(`hsp.rs:1488-1496`) and wraps every tabular hit into a `PairwiseHit`
(`blast_engine.rs:4653-4700`).  NCBI stream reads HSP lists directly and only
builds alignment text when formatting needs it.

For outfmt 6 with default fields, write directly from `BlastpHitList`/`BlastpHsp`
without constructing alignment strings or cloning gap scripts.  Keep the
existing rendered path for outfmt 0, `qseq`, `sseq`, and `btop`.

Parity risk: low-medium.  Formatting order and numeric formatting must remain
unchanged.

### P2. Build/profile configuration work

Expected effect: variable; useful after algorithmic hot loops are fixed.

Release already uses `lto = true`, `codegen-units = 1`, `panic = "abort"`, and
`opt-level = 3` in `LOSAT/Cargo.toml`.  After P0/P1 changes, benchmark:

- default release;
- `RUSTFLAGS="-C target-cpu=native"` release;
- PGO, if the build pipeline can support it without changing runtime behavior.

Parity risk: low, but floating-point output must be rechecked because CBS and
E-value code is sensitive.

## Measurement checklist for the next sweep

Add diagnostics before major rewrites:

- subject residues scanned and offset pairs emitted;
- one-hit/two-hit extensions attempted and saved;
- `init_hsps` count before sort, after chaining, and after cutoff;
- preliminary exact/restricted DP calls and DP cells visited;
- restricted exact-retry count and interval-tree rebuild count;
- kappa local matches considered, early-terminated, redone, and heap-inserted;
- subject/query range bytes copied in CBS redo;
- per-thread subject count, local-match count, and final heap sizes.

Suggested benchmark matrix:

| Case | Why |
|---|---|
| default `blastp`, CBS 2, `-num_threads 1` | single-thread hot-loop baseline |
| default `blastp`, CBS 2, `-num_threads N` | exposes serial kappa and merge overhead |
| default `blastp`, `-comp_based_stats 0` | isolates preliminary search from CBS redo |
| `blastp-fast`, `-num_threads 1/N` | isolates compressed lookup build and scan |
| outfmt 0 and outfmt 6 | separates final rendering from search |

Acceptance remains bit-perfect parity with NCBI output.  Timing improvements
should be accepted only after parity comparisons pass for hit counts, raw
scores, bit scores, E-values, coordinates, identities, positives, and output
ordering.

## Next session handoff

Skip follow-up work for the completed compressed lookup P1 item.  Do not spend
the next session on extra diff review, parity comparison setup, or repeated
benchmarking for that item unless explicitly requested.

Proceed directly to the remaining P1 work, starting with "Reduce per-redo range
copying and SEG work" in `LOSAT/src/algorithm/blastp/kappa.rs`.  Begin by
adding the requested counters for query/subject range bytes copied and SEG range
work, then use those counters to guide the direct subject-range construction
optimization against the NCBI reference in `blast_kappa.c:1573-1646`.
