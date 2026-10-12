# LOSAT debug and diagnostics environment variables

Moved from `AGENTS.md` on 2026-10-08 so that it is read only when needed.

- `LOSAT_TRACE_HSP="qstart,qend,sstart,send"` trace a specific TBLASTX HSP.
- `LOSAT_TRACE_HSP_MASKS=1` print mask coverage for the traced TBLASTX HSP.
- `LOSAT_TRACE_CHAIN_HSP="qstart,qend,sstart,send"` trace TBLASTX chain
  selection for a specific HSP.
- `LOSAT_TRACE_LINK_SELECTIONS=1` print TBLASTX link-selection details.
- `LOSAT_LINKING_LEGACY=1` TBLASTX sum-statistics linking by the literal port
  `link_hsp_group_ncbi` instead of the default kernel (the incremental kernel of
  `tblastx/sum_stats_linking/linking_incr.rs`; same result). The linking traces
  (`LOSAT_TRACE_HSP`, `LOSAT_TRACE_CHAIN_HSP`, `LOSAT_TRACE_LINK_SELECTIONS`,
  `LOSAT_DEBUG_CHAINING`) also run the literal port, which prints them.
- Checks of the default TBLASTX linking kernel: `LOSAT_LINK_FAST_VERIFY=1` (every
  pass against a pass computed from scratch by NCBI's scans; slow on large
  groups), `LOSAT_LINK_FAST_SHADOW=1` (every group also linked by the literal
  port, all fields compared), `LOSAT_LINK_WIDE_KEYS=1` (u128 keys for every
  group, to measure their cost), `LOSAT_LINK_STATS=1` (per-group counters,
  `[LINK_INCR_STATS]` lines).
- `LOSAT_DUMP_TBLASTX_STAGE=<dir>` append TBLASTX stage snapshots as TSV files.
- `LOSAT_TRACE_BLASTN_HSP="qstart,qend,sstart,send"` trace a specific BLASTN
  HSP.
- `LOSAT_TRACE_BLASTN_SEED="q,s"` trace a specific BLASTN seed.
- `LOSAT_TRACE_BLASTN_CONTEXT=<context_idx>` restrict BLASTN tracing by context.
- `LOSAT_TRACE_BLASTN_SUBJECT=<subject_id_or_index>` restrict BLASTN tracing by
  subject.
- `LOSAT_TRACE_BLASTN_STAGE=<seed|ungapped|prelim|traceback|purge|hitlist|all>`
  restrict BLASTN tracing by stage.
- `LOSAT_DEBUG_CUTOFFS=1` cutoff calculations (tblastx + blastn).
- `LOSAT_DEBUG_CUTOFFS_ALL=1` verbose TBLASTX cutoff diagnostics.
- `LOSAT_DEBUG_CHAINING=1` chaining debug (legacy; tblastx).
- `LOSAT_DEBUG_EXTENSION=1` tblastx extension debug.
- `LOSAT_DEBUG_HSP_SAVING=1` TBLASTX HSP-save diagnostics.
- `LOSAT_DEBUG_OUTPUT_FILTER=1` TBLASTX output filter diagnostics.
- `LOSAT_DEBUG_BLASTN=1` blastn hit loss diagnostics.
- `LOSAT_DEBUG_COORDS=1` blastn coordinate transforms.
- `LOSAT_DEBUG_COORDS_START=<int>` narrow selected BLASTN coordinate diagnostics.
- `LOSAT_DEBUG_SCAN_SOFF=<int>` tblastx scan debug center subject offset.
- `LOSAT_DEBUG_SCAN_WINDOW=<int>` tblastx scan debug window size.
- `LOSAT_TIMING=1` timing breakdown.
- `LOSAT_DIAGNOSTICS=1` general diagnostics counters.
- `LOSAT_STARTUP_TRACE=1` startup trace.
- `LOSAT_WASI_THREADS_DEBUG=1` threaded-WASI scheduling diagnostics.
- `LOSAT_WASI_PRESPAWN=auto|N` (Node threaded runner, opt-in) boot N thread
  workers (`auto`: `-num_threads` - 1) before `_start`; same output, extra
  `standby_*` debug events. See `LOSAT/tests/README.md`, "Threading contract gates".
- `LOSAT_WASI_GUARD_CACHE=<dir>` (Node threaded runner, opt-in) reuse the
  shared-memory guard's output from `<dir>`, checked by SHA-256 and rebuilt
  when it does not match.
- `LOSAT_TBLASTX_PARALLEL_CHUNKS=1` force TBLASTX subject-chunk parallel path
  for diagnostics.
- `LOSAT_TBLASTX_SERIAL_SCAN_CHUNKS=1` diagnostic-only sequential TBLASTX
  scan-interior chunking; this does not enable parallel scan work.

## Performance switches (on by default)

Ported from the round 3-5 performance work (2026-10, `perf/strict-set` and the packages after it).
Each switch of the strict set gives the reference result by construction (the same integer
results, or the same floating-point operations in the same order, or only a change of order,
placement or reuse), so none of them changes the output. Until main 62102f7a every switch was
off unless its variable was set. Since `perf/strict-default` (2026-10-12) the 43 switches of the
strict set are on by default, and their variables mean (`LOSAT/src/utils/xswitch.rs`):

| `LOSAT_X_<NAME>` | `LOSAT_X_DEFAULT` | the switch is |
| --- | --- | --- |
| `0` | any | off |
| any other non-empty value | any | on (a value switch takes the number if it is in range, otherwise its built-in value) |
| not set, or empty | `0` | off, the default before the change |
| not set, or empty | not set, empty, or anything but `0` | on |

- `LOSAT_X_DEFAULT=0` runs the previous (reference) paths everywhere. With some variables set
  as well, it gives the run those variables gave before the default changed.
- `LOSAT_X_<NAME>=0` turns one switch off, to find which path a problem comes from. Single
  switches off are a diagnostic aid: `DPAVX2` and `DPW128` act only with `DPFAST`; `TBNSSIDE`,
  `TBNQSIDE` and `TBNBATCH` only with `TBNPAR` (`TBNEVENTS` only without it); `TBNBATCH` needs
  `TBNSSIDE`; `BXBATCH` needs `BXPOOL`; `AHEAD` needs `GREEDYSPEC` for megablast; `MBBATCH` runs
  where `MBDELAY` does not.
- White space around a value is ignored (`set LOSAT_X_DEFAULT=0 && ...` in cmd.exe stores `0 `).
  Only `0` means off; `00`, `false` or `off` are "any other value".
- The variables are read once per process (once per Wasm instance), except `GREEDYSPEC`,
  `MBDELAY`, `MBBATCH` and `TBNBATCH`, which are read once per search or lookup build. The Web
  reactors get an empty WASI environment, so they always run the default.
- The `...SHADOW` / `...VERIFY` switches run the new and the reference computation side by side
  and stop with an assertion on the first difference. They, the counters, the timing aids and the
  tuning variables are not part of the set: they are off unless set to any value, and a shadow
  switch turns its switch on whatever the two variables above say (except
  `LOSAT_X_AHEADSHADOW`, which compares only while `AHEAD` is on). A shadow compares the kernel
  the default selects (for `LOSAT_X_DPSHADOW` the 256-bit, `DPW128` or base kernel that takes the
  call).
- The names are provisional (`x_` prefixes, `EXPERIMENT` comments). The 43 single variables are
  to be removed once the default has run on every CI environment and through a few releases;
  `LOSAT_X_DEFAULT` stays. In the descriptions below `LOSAT_X_<NAME>=1` shows the switch, not a
  setting that is needed, and "strict-set candidate" marks the round-5 switches, which are in
  the set.

The strict set (the variables that are on; `STRICT_SET` in `xswitch.rs`), written as the
environment that turned it on before the change. The three numbers are the built-in values:

```bash
export LOSAT_X_DPFAST=1 LOSAT_X_COMPFAST=1 LOSAT_X_PURGEFAST=1 LOSAT_X_ITREEFAST=1 LOSAT_X_CTXFAST=1
export LOSAT_X_GREEDYFAST=1 LOSAT_X_GREEDYSPEC=1 LOSAT_X_AHEAD=16 LOSAT_X_ENVCACHE=1
export LOSAT_X_DUSTFAST=1 LOSAT_X_DUSTRING=1 LOSAT_X_MBDELAY=1 LOSAT_X_MBBATCH=64
export LOSAT_X_SEGFAST=1 LOSAT_X_SEGMEMO=1 LOSAT_X_SEGSHARE=1 LOSAT_X_KARLINFAST=1 LOSAT_X_IDEALMEMO=1
export LOSAT_X_ERFMEMO=1 LOSAT_X_DEKKER=1 LOSAT_X_LUTARENA=1 LOSAT_X_LUTSPLIT=1
export LOSAT_X_TBNPAR=1 LOSAT_X_TBNEVENTS=1 LOSAT_X_TBNQSIDE=1 LOSAT_X_TBNSSIDE=1 LOSAT_X_TBNBATCH=1
export LOSAT_X_BXPAR=1 LOSAT_X_BXLAZYCTX=1 LOSAT_X_BXLEAN=1 LOSAT_X_BXPOOL=1 LOSAT_X_BXCHUNK=64 LOSAT_X_BXBATCH=1
export LOSAT_X_CODONFAST=1
export LOSAT_X_NEWTONEXACT=1
export LOSAT_X_DPAVX2=1 LOSAT_X_DPW128=1
export LOSAT_X_BXLUT=1 LOSAT_X_BXSCAN=1 LOSAT_X_BXSEGMEMO=1
export LOSAT_X_PAIRPAR=1 LOSAT_X_NEWTONLANES=1 LOSAT_X_GREEDYPRELIM=1
```

Value switches: `LOSAT_X_AHEAD` (2 to 65,536, built-in 16), `LOSAT_X_MBBATCH` (2 to 4,096,
built-in 64), `LOSAT_X_BXCHUNK` (1 or more, built-in 64).

`LOSAT_X_SEEDBUCKET` was in the block of the round 3-4 set and is not on by default (it stays
"set to any value"). When the set became the default it numbered the lookup hits of one chunk
scan in a 32-bit counter, exact only up to 2^32 hits per scan, which a query and a subject of a
few megabases pass. It now numbers the hits per flush, which is exact for any number of hits. It
stays opt-in until it is measured on the searches it would act on with the set on (the TBLASTX
searches that `LOSAT_X_PAIRPAR` leaves to the scan loop).

`LOSAT_X_PAIRLUT` (added after the default changed) is not in the set either: it is off unless
set to any value, and acts while `LOSAT_X_PAIRPAR` is on.

What runs per target:

- x86_64: the vector X-drop DP (`DPFAST`: SSE4.1 or AVX 128-bit; `DPAVX2`: 256-bit with AVX2;
  `DPW128`), the AVX2 builds of the greedy rows and of the Newton kernels where the CPU has AVX2,
  `MBDELAY` with the prefetch instruction.
- wasm32 (built with simd128): the simd128 X-drop DP (`DPW128` replaces the first kernel),
  `MBBATCH` instead of `MBDELAY`, `DEKKER`; everything else is the same Rust as on x86_64
  without AVX2.
- aarch64 and any other target: no vector X-drop DP is compiled, so `DPFAST`, `DPAVX2` and
  `DPW128` do nothing (the scalar port runs); `MBBATCH` instead of `MBDELAY`; `DEKKER` does
  nothing (it is compiled for wasm32 only); everything else as on wasm32.
- One thread (`-num_threads 1`, a build without the `parallel` feature, wasm32 without
  `wasm-threads`): the scheduling switches (`AHEAD`, `GREEDYSPEC`, `PAIRPAR` for
  blastn/megablast, `SEGSHARE`, `TBNBATCH`, `BXPAR`, `BXPOOL`, `BXBATCH`, `BXCHUNK`) leave the
  reference code running. The unit paths have no thread-count gate: the tblastn `TBNPAR` units
  (with `TBNSSIDE` / `TBNQSIDE`; only its parallel redo needs a pool) and the TBLASTX `PAIRPAR`
  units run in order on the calling thread.

Same values, faster computation:

- `LOSAT_X_DPFAST=1` X-drop gapped DP (score-only and traceback) in 128-bit vectors (16 x i8 or
  8 x i16 lanes; AVX2/SSE on x86_64, simd128 on Wasm), blastn task blastn and the protein
  programs; calls that do not fit 8/16 bits use the original loop (`utils/xdrop_simd.rs`).
  Shadow: `LOSAT_X_DPSHADOW=1`. Timing aids: `LOSAT_X_DPNO8=1` (no 8-bit kernel),
  `LOSAT_X_DPNOAVX=1` (no AVX2).
- `LOSAT_X_DPAVX2=1` (with `LOSAT_X_DPFAST`; strict-set candidate) the same X-drop DP in 256-bit
  AVX2 vectors (32 x i8 or 16 x i16 lanes, x86_64 with AVX2 only): the same values, lane width
  and band steps as the 128-bit kernels, so the same result and edit script; a call the 256-bit
  kernel does not take runs on the 128-bit kernel (`utils/x_xdrop_avx2.rs`). Covered by
  `LOSAT_X_DPSHADOW=1`; `LOSAT_X_DPNOAVX=1` turns it off as well.
- `LOSAT_X_DPW128=1` (with `LOSAT_X_DPFAST`; strict-set candidate) the 128-bit X-drop DP kernels
  with the per-row savings of `LOSAT_X_DPAVX2` (first and last survivor kept while the row runs,
  row gap as a prefix maximum with a two-operation carry, the first lane's missing diagonal as a
  lane mask, column residues copied 16 at a time, row tables of `Rows28` / `Flat` matrices built
  with vector instructions, byte shifts in the form V8 compiles to `palignr`): the same values,
  lane width, blocks and result as the 128-bit kernels; a call it does not take (8-bit lanes with
  `16 * gap_extend - (gap_open + gap_extend) > 128`) runs on the 128-bit kernel
  (`utils/x_xdrop_w128.rs`). On x86_64 it runs where the 128-bit kernel would (after
  `LOSAT_X_DPAVX2` when both are set on an AVX2 machine, which keeps priority), on the SSE4.1 or
  AVX-128 build; on wasm32 it replaces the simd128 128-bit kernel (the only DP kernel there).
  Covered by `LOSAT_X_DPSHADOW=1`.
- `LOSAT_X_GREEDYFAST=1` megablast greedy alignment, one distance in three passes of 256 cells,
  match extension 8 bases at a time. Shadow: `LOSAT_X_GREEDYSHADOW=1`.
- `LOSAT_X_GREEDYPRELIM=1` (strict-set candidate, S-G, round 5) the preliminary (score-only)
  megablast greedy alignment reads a one-byte-per-base copy of the packed ncbi2na subject (the
  NCBI-resolved codes) with `compressed_subject = false` instead of the packed bytes: the same
  comparisons, bounds and results (`algorithm/blastn/alignment/x_greedy_prelim.rs`). Only with
  `LOSAT_X_GREEDYFAST` does the preliminary stage then use the fast rows; without it the ported
  loop runs on the copy. Memory: one byte per subject base while the subject is searched (made on
  the first chunk with ungapped hits, shared by its chunks and threads). Shadow:
  `LOSAT_X_GREEDYPRELIMSHADOW=1` (implies the switch) also runs the packed call and stops on the
  first difference; at exit `[X_SHADOW] GREEDYPRELIM calls_compared=N`. With `LOSAT_X_STATS=1`:
  `[X_STATS] GREEDYPRELIM calls= subjects_unpacked= bytes_unpacked=`.
- `LOSAT_X_COMPFAST=1` composition adjustment (Newton) loops reordered with the same element
  operations in the same order. Shadow: `LOSAT_X_COMPSHADOW=1`.
- `LOSAT_X_NEWTONEXACT=1` the Newton iteration of the composition adjustment with fixed-size
  storage, a structured right-looking column Cholesky and AVX2, same element operations
  (`core/composition_adjustment/x_newton_exact.rs`). Shadow: `LOSAT_X_NEWTONEXACTSHADOW=1`
  (convergence state and the 400 values compared bitwise). `LOSAT_X_NEWTONEXACT_LEVEL=0|2`
  (scalar | try AVX-512), `LOSAT_X_NEWTONEXACT_PROF=1` (cycle counters per phase).
- `LOSAT_X_NEWTONLANES=1` (blastp and tblastn) Newton problems of the composition adjustment
  solved ahead, four problems side by side (one per AVX2 lane, each lane with the operations of
  `x_newton_exact`; `core/composition_adjustment/x_newton_lanes.rs`): blastp the first adjustment of
  the next 16 matches of a query (serial and query-parallel redo; the single-query match-parallel
  redo is unchanged), tblastn the adjustments of the coming windows of a match (at least 32 problems;
  the window ranges, translation + SEG, computed ahead are handed to the redo loop). A call uses a
  stored result only when its input (q, row and column probabilities, relative entropy) has the same
  bits; otherwise the existing path runs. Shadow: `LOSAT_X_NEWTONLANESSHADOW=1` (every stored result
  that is used is compared with the reference port, status and the 400 values bitwise; every kept
  tblastn range with a fresh `get_range`). With `LOSAT_X_STATS=1` it prints
  `[X_STATS] NEWTONLANES prefetched= hits= misses= unused= ...` (no `xstats` feature needed).
  Strict-set candidate (S-F, round 5).
- `LOSAT_X_SEGFAST=1`, `LOSAT_X_SEGMEMO=1` SEG window as a residue histogram updated in O(1),
  entropy memoised per state vector. Shadow (`SEGFAST` path): `LOSAT_X_SEGSHADOW=1`.
- `LOSAT_X_DUSTFAST=1`, `LOSAT_X_DUSTRING=1` DUST perfect-interval list built in one merge,
  window kept in a ring buffer. Shadow: `LOSAT_X_DUSTSHADOW=1`.
- `LOSAT_X_KARLINFAST=1` Karlin K convolution in eight lanes (same per-lane order).
- `LOSAT_X_ERFMEMO=1` Spouge ErfC memoised per thread by the argument's bits.
- `LOSAT_X_DEKKER=1` (Wasm only) ErfC product error by Dekker's product instead of `fma`.
- `LOSAT_X_CODONFAST=1` an unambiguous codon translated by one table read.
- `LOSAT_X_PURGEFAST=1` blastn common-endpoint purge without shifting the array per removal.
- `LOSAT_X_ITREEFAST=1` blastn interval tree with side indexes (grid, endpoint bitmaps) that
  answer "none contained" / "no such endpoint" without walking the tree. Shadow:
  `LOSAT_X_ITREESHADOW=1`.
- `LOSAT_X_CTXFAST=1` blastn query-offset-to-context lookup by binary search instead of a per-base
  table. Shadow: `LOSAT_X_CTXSHADOW=1`.
- `LOSAT_X_LUTARENA=1` protein lookup-table chains carved from one arena per table.
- `LOSAT_X_MBDELAY=1`, `LOSAT_X_MBBATCH=64` megablast lookup construction touches (prefetches)
  hash cells ahead in batches of 64; same update order.
- `LOSAT_X_ENVCACHE=1` blastn debug variables read once instead of per extension.
- `LOSAT_X_BXLAZYCTX=1` blastx seed context lookup moved after the diagonal test.
- `LOSAT_X_BXLEAN=1`, `LOSAT_X_TBNEVENTS=1` do not build records nobody reads (candidate copies
  without an observer, test events); no tree or DP space for a (chunk, subject) without HSPs.
- `LOSAT_X_BXLUT=1` blastx lookup table filled in place from per-cell counts (same append order as
  `s_AddNeighboringWords` and `BlastAaLookupFinalize`), neighbour-word lists memoised per thread
  (`tblastx/lookup/x_lut_direct.rs`). Shadow: `LOSAT_X_BXLUTSHADOW=1` (every table compared field by
  field). Strict-set candidate (S-C, round 5).
- `LOSAT_X_BXSCAN=1` blastx seed loop reads the scanner's per-thread offset-pair buffer in place instead
  of a fresh buffer and a copy per subject (`blastx/x_seed_scan.rs`); not in the diagnostic stages, and
  only with the reference diagonal order (`LOSAT_X_SEEDBUCKET` mode 0). Shadow: `LOSAT_X_BXSCANSHADOW=1`
  (hit lists and diagonal table). Strict-set candidate (S-C, round 5).
- `LOSAT_X_SEEDBUCKET=1` TBLASTX two-hit stage: hits buffered per diagonal range and processed
  range by range, HSPs restored to scan order (`tblastx/x_seed_bucket.rs`); only for diagonal tables
  of at least `LOSAT_X_SEEDBUCKET_MIN_CELLS` cells (default 2^21). Shadow:
  `LOSAT_X_SEEDBUCKETSHADOW=1` (HSP lists and the whole diagonal table; its summary line gives
  `max_scan_hits`, the most lookup hits one scan went through). Tuning:
  `LOSAT_X_SEEDBUCKET_CELLS` (cells per bucket), `LOSAT_X_SEEDBUCKET_BUDGET` (hits per flush;
  values above 2^31 are read as 2^31, and memory, 12 bytes per hit reserved up front, limits a
  real setting long before that). A shadow run compares something only where the scan loop runs:
  with the default set use `LOSAT_X_PAIRPAR=0` and check that `chunks_compared` is not 0.
  BLASTX needs `LOSAT_X_BXSEEDBUCKET=1` as well (no gain measured).

Rebuilding less:

- `LOSAT_X_SEGSHARE=1` blastp with 2+ threads: the subject SEG masks shared by the whole search,
  as in the one-thread path. Shadow: `LOSAT_X_SEGSHARESHADOW=1`.
- `LOSAT_X_TBNQSIDE=1` tblastn: query SEG and lookup table once per query set.
- `LOSAT_X_TBNSSIDE=1` tblastn: subject six-frame translations once per search (up to 256 MiB).
- `LOSAT_X_LUTSPLIT=1` tblastn: build only the half (contexts or table) a caller uses.
- `LOSAT_X_IDEALMEMO=1` BLOSUM62 ideal Karlin block once per process.
- `LOSAT_X_BXSEGMEMO=1` blastx redo stage: the subject SEG once per subject instead of once per window
  (per-thread memo of the last subject, byte-equal input; `blastx/x_seg_memo.rs`). Shadow:
  `LOSAT_X_BXSEGMEMOSHADOW=1` (the masker also runs on every call). Strict-set candidate (S-C, round 5).

Scheduling (each value computed by one thread with the original function; order of output kept):

- `LOSAT_X_AHEAD=16` blastn/megablast ordered gapped-extension and traceback loops: helpers compute
  the next elements' alignments ahead, the tree is read and written in the original order by one
  thread (`utils/xahead.rs`). Shadow: `LOSAT_X_AHEADSHADOW=1`. `LOSAT_X_GREEDYSPEC=1` adds the
  megablast greedy traceback; `LOSAT_X_SPECBATCH` tunes the batch.
- `LOSAT_X_PAIRPAR=1` (strict-set candidate; S-E, round 5) blastn/megablast, one query-subject
  pair on a pool of 2+ threads. (1) The megablast lookup table is filled by the pool threads: every
  thread reads the whole query in order and applies the update only to the words whose cell lies in
  its blocks of `2^max(11, pv_array_bts)` cells, so each cell gets its words in the original order
  (`blastn/x_mb_lookup_par.rs`; not with database word counts or `BLEMIR_DEBUG`). (2) The one-hit
  seed stage of a subject chunk with the diagonal hash (query block over 8000 letters,
  `window_size` 0, `scan_range` 0, contiguous megablast table, unmasked subject, one subject, no
  small-query word extension, no debug or trace): the subject scan grid is cut into strips that
  the pool threads scan (lookup chains and word extension), and the owner thread folds the strips
  in scan order through a copy of the diagonal-hash code (`blastn/blast_engine/x_pair_seed.rs`).
  Anything else runs the original code. Shadow: `LOSAT_X_PAIRPARSHADOW=1` (lookup tables and, per
  chunk, the hit list and the whole diagonal hash compared with the original code; one
  `[X_PAIRPARSHADOW]` summary line on stderr at the end). With `LOSAT_TIMING` a
  `[TIMING] x_pairpar_seed` line gives the strips, seeds and the owner's fold, scan and wait time.
- `LOSAT_X_TBNPAR=1` tblastn: preliminary stage per (frame, chunk) and redo per match in
  parallel, one pool for the search. Shadow: `LOSAT_X_TBNPARSHADOW=1`.
- `LOSAT_X_TBNBATCH=1` (needs `TBNSSIDE`) tblastn query batches searched side by side, appended
  in input order.
- `LOSAT_X_PAIRPAR=1` TBLASTX seed stage (scan, two-hit test, ungapped extension) of each subject
  as units of (subject frame, chunk, query context), 36 per one-query pair, on the search pool (in
  order on the calling thread with one thread); each unit scans with a lookup table holding only its
  context's chains (order kept) on a fresh unmasked diagonal table, and the unit lists go through the
  reference's sort and merge in the reference order (`tblastx/x_pair_par.rs`, exactness argument in
  its header). Only when every query context, or every subject chunk, is at most `slack` = diagonal
  table size minus query length long, and `slack` plus twice the longest subject chunk stays under
  2^29 (it does unless `-window_size` is near 2^29); otherwise, and with the trace, debug,
  diagnostics and chunk switches, the reference loop runs (`LOSAT_X_STATS=1` counts `fallback_guard` / `fallback_gate`;
  `LOSAT_TIMING` prints one `x_pairpar` line). Strict-set candidate (S-E, round 5).
- `LOSAT_X_PAIRPARSHADOW=1` shadow of `LOSAT_X_PAIRPAR`: the unit path and the reference loop for
  every subject, the subject's combined ungapped HSP list compared field by field (the first
  difference stops the run with subject, frame, index and both values); one `[X_SHADOW] PAIRPAR`
  summary line at exit (subjects, units, HSPs compared, fallbacks).
- `LOSAT_X_PAIRLUT=1` (not in the strict set: off unless set; acts while `LOSAT_X_PAIRPAR` is on,
  not with its shadow) TBLASTX:
  the lookup tables of the unit path, one per context group, built directly from the query instead
  of being filtered out of the batch's full table, which is then not built: for each group the
  exact-match index of the group's lookup locations, then the count, layout and fill passes of
  `LOSAT_X_BXLUT` (`tblastx/lookup/x_pair_lut.rs`, exactness argument in its header), the groups
  in order on the calling thread (tuning only: `LOSAT_X_PAIRLUT_POOL=1` builds them on the search
  pool, a few ms faster per batch and more memory). The reference code holds a stub (the full
  table's header fields and `longest_chain`, no entry). Whether the unit path runs is decided before the lookup
  build by the function `prepare` decides with; the run stops if `prepare` disagrees or if the
  reference frame loop would scan the stub. The reference builder and `build_groups` run as before
  when the unit path does not run (gate, guard), in diagnostics mode, and when a lookup location
  does not lie inside one group. `LOSAT_X_STATS=1` prints
  `[X_STATS] PAIRLUT batches_direct= groups= fallback=`; `LOSAT_TIMING` prints one `x_pairlut` line
  (the direct build; it is part of `build_lookup`, and `sublookup` of the `x_pairpar` line is 0).
- `LOSAT_X_PAIRLUTSHADOW=1` shadow of `LOSAT_X_PAIRLUT` (while `LOSAT_X_PAIRPAR` is on, or with its
  shadow): the
  full table and `build_groups` as before, and the direct groups as well; every group table (cells,
  overflow array, presence vector, header, `longest_chain`), `q_min` / `q_max`, the stub's header
  and `longest_chain`, and the contexts are compared, and the first difference stops the run; the
  search uses the reference tables. One `[X_SHADOW] PAIRLUT batches= groups_compared= fallback=`
  line at exit.
- `LOSAT_X_BXPAR=1` blastx query-chunk, redo and per-context SEG parallelism;
  `LOSAT_X_BXPOOL=1` one pool for the search; `LOSAT_X_BXBATCH=1` (needs `BXPOOL`) small query
  batches searched together, output in input order; `LOSAT_X_BXCHUNK=64` subjects per worker;
  `LOSAT_X_BXWAVE` tuning only.

Counters and tests: `LOSAT_X_STATS=1` (counters on stderr; DP and Newton counts need the
`xstats` cargo feature), `LOSAT_X_ADJMEMO_STATS=1` (distinct composition-adjustment inputs),
`LOSAT_FUZZ_CASES=<n>` (random groups of the linking fuzz test; random problems of the X-drop DP
kernel tests).

X-drop DP capture and replay (diagnostics, not part of the strict set):

- `LOSAT_X_DPCAPTURE=<dir>` records every Nth call of each thread that reaches the vector X-drop
  DP (`LOSAT_X_DPFAST`) as a self-contained problem (residues, matrix, gap costs, x_drop, result
  and edit script) in `<dir>/xdp-<pid>-t<k>.xdp`, one file per thread, appended
  (`utils/x_xdrop_capture.rs`, format in its header). `LOSAT_X_DPCAPTURE_EVERY=<n>` sets N
  (default 25). The search output does not change.
- `LOSAT_X_DPREPLAY=<dir>` input of the ignored test `replay_captured_problems`
  (`cargo test --release replay_captured_problems -- --ignored --nocapture`): every captured
  problem through the scalar loops, the 128-bit, the `LOSAT_X_DPW128` and the 256-bit kernels and
  the dispatch, compared with the recorded result, then the three kernels timed (median ns per cell
  and per call, traceback and score-only, 8- and 16-bit lanes). `LOSAT_X_DPREPLAY_REPS` (default 7),
  `LOSAT_X_DPREPLAY_OUT=<file>` (append the table to a file). The ignored test
  `w128_replay_captured_problems` does the same for the 128-bit and the `LOSAT_X_DPW128` kernels
  only and also runs in the wasm32-wasip1 test build (the libtest `.wasm` under a WASI host such
  as Node's `wasi` module with a root preopen), which times the simd128 kernels.
