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
- `LOSAT_TIMING=1` timing breakdown. TBLASTX prints, per query batch, after its other lines:
  `[TIMING] x_search: setup .. subjects .. collect ..` (three wall times that add up to
  `search_total`: before the first subject, the subject loop, after it) and
  `[TIMING] x_subject: enter .. translate .. prepare .. units .. units_release .. merge .. frame_loop .. seed_release .. link1 .. reap .. reevaluate .. link2 .. convert .. emit .. end .. sum .. (subjects=N)`
  (the stages of the per-subject work in the order they run; one clock per subject is read at
  every stage boundary, so the stages add up to `sum`, and `end` includes the release of the
  subject's translation; with several subjects the times are sums over the subjects and their
  threads). `link1` is the first linking call, on every HSP the seed stage saved; `link2` is the
  call `sum_stats_linking` times (`tblastx/x_search_timing.rs`). Two more lines follow, one per
  linking call:
  `[TIMING] x_link1: reset .. sort .. groups .. kernels .. assemble .. final_sorts .. replay .. sum .. (calls=N hsps=N groups=N largest_group=N kernel_sum ..s kernel_max ..s)`
  and `[TIMING] x_link2:` in the same form (the steps of the default linking call,
  `x_hsp_inplace::link_in_place`, in the order they run; one clock per call, so the steps add up
  to `sum`. `kernels` is the wall time of the per-group linking, `kernel_sum` the per-group
  times added over the threads that ran them, `kernel_max` the longest group, which `kernels`
  cannot go below; `largest_group` is the number of HSPs in the largest group. A linking call
  that takes the reference body (`LOSAT_X_HSPINPLACE=0`, `LOSAT_X_DEFAULT=0`,
  `LOSAT_DIAGNOSTICS`, `LOSAT_DEBUG_CUTOFFS`) is not timed there: `calls=0`;
  `tblastx/x_link_timing.rs`). The clock is not read without the variable.
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
off unless its variable was set. Since `perf/strict-default` (2026-10-12) the switches of the
strict set are on by default (43 then, 48 since `LOSAT_X_PAIRLUT`, `LOSAT_X_HSPINPLACE`,
`LOSAT_X_SEEDBUCKET`, `LOSAT_X_TBNOIDCHUNK` and `LOSAT_X_SEEDRELEASE` joined), and their
variables mean (`LOSAT/src/utils/xswitch.rs`):

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
  where `MBDELAY` does not; `PAIRLUT` acts only while `PAIRPAR` is on; `SEEDBUCKET` acts only
  in the scan loop, which TBLASTX runs where `PAIRPAR` does not take the search.
- White space around a value is ignored (`set LOSAT_X_DEFAULT=0 && ...` in cmd.exe stores `0 `).
  Only `0` means off; `00`, `false` or `off` are "any other value".
- The variables are read once per process (once per Wasm instance), except `GREEDYSPEC`,
  `MBDELAY`, `MBBATCH` and `TBNBATCH`, which are read once per search or lookup build. The Web
  reactors get an empty WASI environment, so they always run the default.
- The `...SHADOW` / `...VERIFY` switches run the new and the reference computation side by side
  and stop with an assertion on the first difference. They, the counters, the timing aids and the
  tuning variables are not part of the set: they are off unless set to any value, and a shadow
  switch turns its switch on whatever the two variables above say (except
  `LOSAT_X_AHEADSHADOW`, which compares only while `AHEAD` is on, and
  `LOSAT_X_TBNXLATESHADOW`, which compares only where `TBNPAR` makes the subject plans). A
  shadow compares the kernel the default selects (for `LOSAT_X_DPSHADOW` the 256-bit, `DPW128`
  or base kernel that takes the call).
- The names are provisional (`x_` prefixes, `EXPERIMENT` comments). The 48 single variables are
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
export LOSAT_X_PAIRLUT=1 LOSAT_X_HSPINPLACE=1
export LOSAT_X_SEEDBUCKET=1
export LOSAT_X_TBNOIDCHUNK=1
export LOSAT_X_SEEDRELEASE=1
```

Value switches: `LOSAT_X_AHEAD` (2 to 65,536, built-in 16), `LOSAT_X_MBBATCH` (2 to 4,096,
built-in 64), `LOSAT_X_BXCHUNK` (1 or more, built-in 64).

`LOSAT_X_SEEDBUCKET` was in the block of the round 3-4 set and was left out when the set became
the default: it numbered the lookup hits of one chunk scan in a 32-bit counter, exact only up to
2^32 hits per scan, which a query and a subject of a few megabases pass. It now numbers the hits
per flush, which is exact for any number of hits, and joined the set on 2026-10-12. With the set
on it acts on the TBLASTX searches that `LOSAT_X_PAIRPAR` leaves to the scan loop and whose
diagonal table has at least 2^21 cells (a query of 524 kb or more): a query of 0.9-1.05,
1.8-2.1, 3.6-4.2 or 7.2-8.4 Mb against a subject longer than the table's spare cells
(`x_pair_par.rs` `guard`), and on threaded Wasm also a subject of more than one chunk (more than
15 Mb) at more than one thread. Measured there on 2026-10-12 (alternating pairs with and without
the switch, the same output in every run; Linux x86_64 with a 36 MiB L3, node 26 for Wasm), as
time with the switch over time without. On a quiet machine: a 1.0 Mb query against 1.0 Mb (the
smallest table it acts on, 8 MiB) at 8 threads 0.84-1.12 natively (median 0.88 of five pairs) and
0.77-0.98 on threaded Wasm, at one thread 0.85-1.09 natively (median 0.99) and 0.79-0.96 on
serial Wasm; a 1.0 Mb query against 5.0 Mb 0.99 and 1.08 (native, 8 threads); a 2.0 Mb query
against 5.0 Mb (16 MiB table) 0.62 and 0.67. With other work on the machine the scan order slows
down more than the bucketed one: most pairs at 1.0 Mb were then between 0.40 and 0.89, single
ones up to 1.44. A 0.6 Mb query against a 15.1 Mb subject on threaded Wasm with a 4 GiB memory:
0.76 and 0.85 (loaded machine). The Wasm memory reaches the same size with and without it.

`LOSAT_X_PAIRLUT` and `LOSAT_X_HSPINPLACE` were added after the default changed, as switches
that were off unless set, and joined the set on 2026-10-12. Both lower the memory of TBLASTX:
the unit path of `LOSAT_X_PAIRPAR` then builds only the tables it scans with, and the HSP list
of a subject is sorted and linked without a second whole copy. `LOSAT_X_PAIRLUT` acts only
while `LOSAT_X_PAIRPAR` is on.

`LOSAT_X_TBNOIDCHUNK` was added after the default changed, as a switch that was off unless set,
and joined the set on 2026-10-13. It changes only which subjects one entry into the tblastn
search pool takes (2 or more threads, several subjects) and needs no other switch. Since
2026-10-14 its chunks are bounded by residues, as NCBI's database reader cuts them for its
search threads.

`LOSAT_X_SEEDRELEASE` was added after the default changed, as a switch that was off unless set,
and joined the set on 2026-10-13. It lowers the peak memory of TBLASTX further: the tables that
only the seed stage reads are freed before the HSP list of a query batch's last subject is
merged and linked. It changes no value and no order of computation, acts on every TBLASTX path
(the unit path of `LOSAT_X_PAIRPAR` and the scan loop) and needs no other switch.

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
  blastn/megablast, `SEGSHARE`, `TBNBATCH`, `TBNOIDCHUNK`, `TBNOIDNEXT`, `TBNOIDUNITS`, `BXPAR`,
  `BXPOOL`, `BXBATCH`, `BXCHUNK`) leave the reference code running. The unit paths have no thread-count gate: the tblastn `TBNPAR` units
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
- `LOSAT_X_GREEDYFUSE=1` (not in the strict set; off unless set) the score-only megablast greedy
  alignment on an unpacked subject (the preliminary stage with `LOSAT_X_GREEDYPRELIM`): the whole
  loop over distances of `BLAST_GreedyAlign` in one function, one pass over the diagonals of a
  distance, eight bases per slide step, the same cells of the two rolling rows and of `max_score`
  written in the C order (`algorithm/blastn/alignment/x_greedy_fused.rs`, exactness argument in
  its module doc). It replaces the per-distance `LOSAT_X_GREEDYFAST` row there; traceback, packed
  subjects and the affine alignment keep their paths. x86_64 with AVX2/BMI runs a build with
  `tzcnt`/`lzcnt`, other CPUs and targets the plain build of the same code. Shadow:
  `LOSAT_X_GREEDYFUSESHADOW=1` (implies the switch) runs it on copies, runs the reference call as
  well and stops on the first difference in the return value, align lengths, seed, fence flag,
  rows or `max_score`; at exit `[X_SHADOW] GREEDYFUSE compared=N`.
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
  only with the reference diagonal order (no `LOSAT_X_BXSEEDBUCKET`, no `LOSAT_X_SEEDBUCKETSHADOW`). Shadow: `LOSAT_X_BXSCANSHADOW=1`
  (hit lists and diagonal table). Strict-set candidate (S-C, round 5).
- `LOSAT_X_SEEDBUCKET=1` (in the strict set since 2026-10-12) TBLASTX two-hit stage: hits buffered
  per diagonal range and processed
  range by range, HSPs restored to scan order (`tblastx/x_seed_bucket.rs`); only for diagonal tables
  of at least `LOSAT_X_SEEDBUCKET_MIN_CELLS` cells (default 2^21). Shadow:
  `LOSAT_X_SEEDBUCKETSHADOW=1` (HSP lists and the whole diagonal table; its summary line gives
  `max_scan_hits`, the most lookup hits one scan went through). Tuning:
  `LOSAT_X_SEEDBUCKET_CELLS` (cells per bucket), `LOSAT_X_SEEDBUCKET_BUDGET` (hits per flush;
  values above 2^31 are read as 2^31, and memory, 12 bytes per hit reserved up front, limits a
  real setting long before that). A shadow run compares something only where the scan loop runs:
  with the default set use `LOSAT_X_PAIRPAR=0` and check that `chunks_compared` is not 0. With
  `LOSAT_X_STATS=1` the line `[X_STATS] SEEDBUCKET bucketed_scans=... flushes=...` tells how many
  scans ran in the bucketed order (`LOSAT/tests/tblastx_seed_bucket.rs` uses it).
  BLASTX runs the bucketed order only with `LOSAT_X_BXSEEDBUCKET=1` as well, which is not in the
  set (off unless set to any value; no gain measured), so BLASTX keeps the reference order by default.

Rebuilding less:

- `LOSAT_X_SEGSHARE=1` blastp with 2+ threads: the subject SEG masks shared by the whole search,
  as in the one-thread path. Shadow: `LOSAT_X_SEGSHARESHADOW=1`.
- `LOSAT_X_TBNQSIDE=1` tblastn: query SEG and lookup table once per query set.
- `LOSAT_X_TBNSSIDE=1` tblastn: subject six-frame translations once per search of four or more
  query batches (20,000 query residues per batch), up to 256 MiB of translated residues. A
  search of one to three batches keeps none and translates a subject again for each batch, as
  NCBI does (`blast_engine.c:627-631`): one batch has no later reader, and two or three give
  the time the kept translations would save for their memory. Timed with two and three batches
  at 4 and 8 threads, with other work on the machine: on 28,036 16S subjects (40.6 Mb) the
  search without them has 115 MiB less (max RSS 218-226 -> 103-113 MiB) and takes a median of
  1.03-1.11 of the time, on one 606 kb genome the same memory and 1.05-1.08 of the time
  (0.01-0.06 s). Four batches at 8 threads gave 1.06 on the 16S set, so the saving does not
  start at four batches: the number is where the memory is preferred to it. `TBNBATCH` needs
  the kept translations, so two or three batches with at least one batch per thread are now
  searched one after another instead of side by side: 1.02-1.07 of the time on the 16S set and
  1.12-1.24 on the genome (0.04-0.09 s) with one batch per thread, 0.75-0.86 with three
  batches on two threads.
- `LOSAT_X_LUTSPLIT=1` tblastn: build only the half (contexts or table) a caller uses.
- `LOSAT_X_IDEALMEMO=1` BLOSUM62 ideal Karlin block once per process.
- `LOSAT_X_BXSEGMEMO=1` blastx redo stage: the subject SEG once per subject instead of once per window
  (per-thread memo of the last subject, byte-equal input; `blastx/x_seg_memo.rs`). Shadow:
  `LOSAT_X_BXSEGMEMOSHADOW=1` (the masker also runs on every call). Strict-set candidate (S-C, round 5).
- `LOSAT_X_HSPINPLACE=1` (TBLASTX, memory; in the strict set) the HSP list of a subject is sorted, cut
  into linking groups, put together and replayed into NCBI's output order without a second whole copy
  alive (`tblastx/x_hsp_inplace.rs`, exactness argument in its header). The index sorts of
  `sort_ungapped_hits_by_index_replay` (score order, link order) move the records inside the list, cycle by
  cycle, after the same `sort_by`; the linking call cuts its groups off the end of the sorted list as exactly
  sized lists, links them with the same per-group functions (serially or on the pool, as before), joins them
  in one allocation of the exact size, and applies the output replay as a permutation in place. A replay that
  leaves an HSP out or repeats one returns the reference's cloned list (`LOSAT_X_STATS=1` counts both:
  `[X_STATS] HSPINPLACE sorts= linkings= hits= groups= replay_inplace= replay_fallbacks= gate_reference=`).
  The per-group scratch of the literal port is not reserved unless that port can run (`LOSAT_LINKING_LEGACY`,
  a linking trace, `LOSAT_LINK_FAST_SHADOW`). A linking call with `LOSAT_DEBUG_CUTOFFS` or
  `LOSAT_DIAGNOSTICS` runs the reference function (`gate_reference`).
- `LOSAT_X_HSPINPLACESHADOW=1` shadow of `LOSAT_X_HSPINPLACE` (implies nothing else): every sort call and
  every linking call runs both ways (the reference on a clone of the input) and the results are compared HSP
  by HSP, every field, E-values by bits; the first difference stops the run and the search continues with the
  reference's result. One line at exit:
  `[X_SHADOW] HSPINPLACE sorts= linkings= hits_compared= replay_fallbacks=`.
- `LOSAT_X_BPLUT=1` (BLASTP, memory; not in the strict set: off unless set, to any value, and
  `LOSAT_X_DEFAULT` does not touch it) the lookup table of the default (non-compressed) BLASTP lookup is
  filled in place from per-cell counts by `LOSAT_X_BXLUT`'s builder (same append order as
  `s_AddNeighboringWords` and `BlastAaLookupFinalize`; `tblastx/lookup/x_lut_direct.rs`, exactness
  argument in its header), so the per-cell chains of the reference builder, about 1.5 times the
  finished table, are never held. In diagnostics mode the reference builder runs. Shadow:
  `LOSAT_X_BPLUTSHADOW=1` (both builders, every table compared field by field, the reference's table
  used); one line at exit: `[X_SHADOW] BPLUT tables_compared=`.
- `LOSAT_X_BPREDOFIT=1` (BLASTP, memory; not in the strict set: off unless set, to any value, and
  `LOSAT_X_DEFAULT` does not touch it) the composition-based redo builds each redone HSP list with room for
  exactly its alignments, and a list from which HSPs were removed afterwards gives back their room
  (`blastp/x_redo_fit.rs`): without it the list's vector grows by `push`, whose first allocation has room for
  four HSPs of 168 bytes, while most lists hold one. The lists hold the same HSPs in the same order, so there
  is no shadow. With `LOSAT_TIMING=1` each search of a query set prints
  `[TIMING] x_bpredofit: lists=<n> hsps=<m>`. Measured once (human proteome block 1, 566 proteins, 50,404
  letters, x human proteome, one thread, with `LOSAT_X_BPQBLOCK=50000` and the other memory switches; max RSS):
  121,304 -> 82,180 kB (NCBI -db 151,476 kB in the same run).
- `LOSAT_X_BPSEED=1` (BLASTP, speed; not in the strict set: off unless set, to any value, and
  `LOSAT_X_DEFAULT` does not touch it) the two-hit seed stage with an amino-acid lookup table
  (`blastp/x_seed_fused.rs`, exactness argument in its header): the subject scan and the diagonal test in one
  loop that reads each lookup cell's list in place, in the reference's order, and hands only the hits that need
  the context check and the ungapped extension to the reference's code; the context of a candidate comes from a
  two-level table (a coarse table of at most 16,384 entries and a walk over the frame bases) instead of the
  binary search; the diagonal entry of a later hit of the cell is prefetched `LOSAT_X_BPSEED_PREFETCH` hits
  ahead (default 16, 0 = none; tuning only). Measured 2026-10-13 on the workstation (OrthoFinder-style job,
  fly x human, every 10th protein, one thread, medians of 6 rotating rounds under load): 51.4 s -> 34.2 s,
  output unchanged. The switch also covers the compressed lookup table of blastp-fast
  (`blastp/x_seed_fused_compressed.rs`, exactness argument in its header): the scan (NCBI's
  `s_BlastCompressedAaScanSubject` statements) lists the subject offsets whose word passes the PV test, then
  reads those cells in scan order, the backbone cell 8 cells ahead prefetched, and hands each cell's hits in
  NCBI's order (the backbone's offsets, then the overflow list from the cell filled last) to the same test and
  closure; the offset-pair array's rounds are not repeated (a resumed round reads the index the rolling scan has,
  as NCBI's reciprocal division is exact below the table size). Measured 2026-10-14 on the WSL machine (blastp-fast,
  2M query letters x 1M subject letters, `LOSAT_X_BPQBLOCK=20000` on both sides, one thread, 8 alternating pairs,
  load 7.5-10.9): `blastp_scan_ungapped` summed over the blocks 2.05 -> 1.40 s (median ratio 0.67), wall
  8.94 -> 8.22 s (0.89), output unchanged.
- `LOSAT_X_BPSEEDSHADOW=1` shadow of `LOSAT_X_BPSEED` (turns the fused loop on), for both lookup tables: for
  every subject the fused loop and the reference loop run from the same diagonal table, and the initial HSPs
  (every field) and the whole table must be equal; the first difference stops the run (panic). For a subject
  shorter than a word less one letter (amino-acid table, one residue with word size 3) only the initial HSPs are
  compared: the reference's only word there ends one byte past the subject's buffer (NCBI reads the same byte),
  whose value can name a cell with hits; such a hit is never extended and leaves an entry that later subjects read
  as no previous hit (`x_shadow_compare`). The compressed scan reads letters past the buffer as the sentinel in
  both loops, so its tables are always compared.
  One line at exit: `[X_SHADOW] BPSEED subjects= initial_hsps= identical`, with
  `(tables not compared for N subjects shorter than a word)` when there were such subjects.
- `LOSAT_X_BPCLUT=1` (BLASTP `-task blastp-fast`, speed and memory; not in the strict set: off unless set, to
  any value, and `LOSAT_X_DEFAULT` does not touch it) the compressed lookup table (`CompressedAaLookupTable`) is
  built flat in three passes (`tblastx/lookup/x_clut_direct.rs`, exactness argument in its header): a count pass
  makes the reference builder's appends in its order (the words of `add_word_hits`, the neighbors' last letter
  taken as a prefix of the sorted row) and records the cell of each append; a layout pass places each cell's hits
  in one slice of one array (an array of starts, and the reference's PV array and longest chain from the counts);
  a fill pass writes the hits in the order of their appends and then moves the hits of each cell of more than 5
  to NCBI's scan order (`x_emission_position`: the first three, the overflow cell filled last, then the older
  cells newest first). The test for NCBI's overflow banks (the per-batch refusal) is the reference's, with the
  same message. Each block of `LOSAT_X_BPQBLOCK` builds its own flat table. The scans read a cell as its slice:
  the offset-pair scan (`scan_subject`'s statements) when `LOSAT_X_BPSEED` is off, the fused loop
  (`blastp/x_seed_fused_clut.rs`) when it is on. The table holds 4 bytes per hit plus 4 per backbone cell (the
  reference: 32 bytes per backbone cell plus 8 per overflowed hit), and its build holds another 4 bytes per hit.
  Measured 2026-10-14 on the WSL machine (blastp-fast, 2M query letters x 1M subject letters,
  `LOSAT_X_BPQBLOCK=20000` and `LOSAT_X_BPSEED=1` on both sides, one thread, 8 alternating pairs, load 5.9-8.2):
  `blastp_build_lookup` summed over the blocks 4.55 -> 2.27 s (median ratio 0.54), `blastp_scan_ungapped`
  1.41 -> 1.08 s (0.84), wall 8.21 -> 5.43 s (0.73), max RSS 77 -> 40 MB, output unchanged.
- `LOSAT_X_BPCLUTSHADOW=1` shadow of `LOSAT_X_BPCLUT` (for blastp-fast): both builders run, the refusal and the
  tables are compared (every cell's hits in the reference's scan order, the PV array, `longest_chain` and the
  scan's parameters) and the first difference stops the run (panic); the search uses the reference's table. One
  line at exit: `[X_SHADOW] BPCLUT tables_compared= refusals_compared= hits= identical`.
- `LOSAT_X_BPUNGAP=1` (BLASTP and BLASTX one-hit and two-hit BLOSUM62 ungapped extension, speed; not in the
  strict set: off unless set, to any value, and `LOSAT_X_DEFAULT` does not touch it) the word loop and the
  X-drop walks of `s_BlastAaExtendLeft` / `s_BlastAaExtendRight` read the scores from one flat 32 x 32 byte
  table (`q << 5 | s`; the sentinel's `BLAST_SCORE_MIN` stored as -128) with raw-pointer walks over the
  positions the reference indexes, NCBI's statements and exit tests unchanged (`blastp/x_ungapped_fast.rs`,
  exactness argument in its header). Acts in `extend_one_hit_blosum62` / `extend_two_hit_blosum62` (TBLASTN
  calls the matrix-generic functions and is not affected); a dropoff outside 1..=127, a word size outside
  1..=11 or a residue code >= 32 run the reference.
- `LOSAT_X_BPUNGAPSHADOW=1` shadow of `LOSAT_X_BPUNGAP` (turns the path on): every extension runs both ways and
  the results (score, offsets, length, `right_extend`, `s_last_off`, and the early return) must be equal; any
  difference panics with the inputs. One line at exit: `[X_SHADOW] BPUNGAP compared=`.
- `LOSAT_X_LINKSORT=1` (TBLASTX, speed; not in the strict set: off unless set, to any value, and
  `LOSAT_X_DEFAULT` does not touch it) the index sorts of the in-place HSP list path
  (`tblastx/x_link_sort.rs`, exactness argument in its header; it acts only where `LOSAT_X_HSPINPLACE` does).
  Two changes, both of how an index order is computed, not of the order: (1) the linking call sorts its
  output order once, stably, with `s_FwdCompareHSPsTransl` and ties by `s_RevCompareHSPsTransl`, where it
  sorts twice (a stable sort by B followed by a stable sort by A is one stable sort by "A, then B"); this
  also saves a sort at one thread. (2) An index sort of 32,768 indices or more runs as rayon's stable
  parallel sort (`par_sort_by`) when its call may use the search pool: the two sorts of a linking call
  under the call's `allow_parallel`, as its per-group linking; the other HSP sorts (the score sorts after
  each frame, after the first linking call and after the re-evaluation) when the calling thread is a worker
  of a rayon pool, which in a search started by the command line or the Web adapter it is exactly when
  `-num_threads` is 2 or more (the header names the one case this repository's code does not show: a
  program that calls a one-thread search from a thread of its own rayon pool). A stable sort with a
  comparator that is a total preorder has one result, so the index order is the same at any number of
  threads; the comparators are NCBI's four, each a lexicographic comparison of integer fields. With one thread, in a
  build without the `parallel` feature and on wasm32 without `wasm-threads` every sort is the serial
  statement. With the switch off the statements of `x_hsp_inplace.rs` run as before. `LOSAT_X_STATS=1`
  prints, only with the switch or its shadow on,
  `[X_STATS] LINKSORT sorts= indices= parallel_sorts= parallel_indices= composite_sorts=`
  (index sorts through the path and their lengths added; those that ran on the pool; final sorts made as
  one). With several subjects on a pool, a thread that waits for the parts of its sort can run another
  subject meanwhile, as it can while it waits for its linking groups, so the `sort`, `final_sorts`,
  `merge` and `reap` times of `LOSAT_TIMING` can then hold time of another subject. Measured once on
  2026-10-13 (medians of three alternating runs without and with the switch, one binary, a loaded
  machine; `[TIMING] x_link1` and `x_subject`, seconds, without -> with): AP027131 x AP027133 (genetic
  code 4, 338,859 HSPs in the first linking call) at 8 threads `sort` 0.085 -> 0.039, `final_sorts`
  0.067 -> 0.009, `reap` 0.073 -> 0.045, `merge` 0.167 -> 0.156, and at one thread `final_sorts`
  0.060 -> 0.035 (the other sorts are today's statements there); NZ_CP006932 against itself (516,395
  HSPs) at 8 threads `sort` 0.139 -> 0.075, `final_sorts` 0.100 -> 0.014, `reap` 0.098 -> 0.068, and at
  one thread `final_sorts` 0.105 -> 0.047. Max RSS did not change (164 -> 161 and 242 -> 243 MiB at 8
  threads, one run each). What is left of `sort` and `reap` includes the move of the records into the
  sorted order, which stays on one thread. In a search of 200 short subjects no list reached 32,768 HSPs
  and no sort went to the pool.
- `LOSAT_X_LINKSORTSHADOW=1` shadow of `LOSAT_X_LINKSORT` (turns the path on): every index sort runs both
  ways, the path's sort on a copy of the index vector and today's statement or statements on the vector
  itself, and the two vectors must be equal element for element; the first difference stops the run with
  the sort (first, final or dispatch), the length and the position, and the search continues with today's
  result. One line at exit: `[X_SHADOW] LINKSORT sorts= indices_compared=`. Together with
  `LOSAT_X_HSPINPLACESHADOW=1` the reference bodies are compared with the path's result end to end
  (`LOSAT/tests/tblastx_link_sort.rs` runs both).
- `LOSAT_X_BPQBLOCK=<letters>` (BLASTP, memory and speed; not in the strict set: off unless set to a positive
  integer, and `LOSAT_X_DEFAULT` does not touch it) the queries are searched in blocks of whole query batches
  (`blastp/x_query_blocks.rs`, exactness argument in its header): a block closes at the first batch end where
  its letters reach the value, and each block is a search of its own (its lookup table, its diagonal table,
  its hit lists) whose rows are written before the next block is searched, as NCBI searches and reports one
  batch at a time. Both lookup tables take it: the default one, and the compressed one of `-task blastp-fast`
  (word size 5), whose table each block builds for its own batches and tests for NCBI's overflow banks with
  them. Taken only where a block's report is its rows: one output in tabular format 6, no `-query_loc`, every
  batch read, more than one block, no subject whose tabular ID the reference path rejects, and none of the
  reference path's other refusals that come before its rows (a `-window_size` whose diagonal table with all
  queries passes 2^30, a query split that LOSAT rejects, or with the compressed table a query batch whose
  offsets add more than four times NCBI's overflow cells, about 859 million words, counted before the search);
  otherwise the search of all queries. The output file is created when the first block's rows are written; a
  run that fails in a later block leaves the rows of the blocks before.
  With `LOSAT_TIMING=1` each block prints its own `[TIMING]` lines.
- `LOSAT_X_BPPFREE=1` (BLASTP, memory; not in the strict set: off unless set, to any value, and
  `LOSAT_X_DEFAULT` does not touch it) the preliminary hit lists hold no more than their HSPs need
  (`blastp/x_prelim_free.rs`, exactness argument in its header): a subject's HSP list for a query gives back
  the room of the HSPs that the preliminary E-value reaped, and the serial composition-based redo (one thread)
  frees each query's preliminary list after its redo, as NCBI's redo frees each list it has redone. The
  parallel redo paths keep the lists, as before. No value changes, so there is no shadow. With
  `LOSAT_TIMING=1` each search of a query set prints `[TIMING] x_bppfree: freed_lists=<n>`. Measured once
  (C. elegans block of 150 proteins, 50,869 letters, x human proteome, one thread, max RSS): 153,744 ->
  117,968 kB (NCBI -db 109,080 in the same run).
- `LOSAT_X_BPPCOMPACT=1` (BLASTP, memory; not in the strict set: off unless set, and `LOSAT_X_DEFAULT` does
  not touch it) the subject scan keeps each preliminary HSP as a 64-byte element (the
  `BlastpPreliminaryHsp` and its E-value) instead of a 168-byte `BlastpHsp`, in the same hit-list helpers
  made generic over the element (`blastp/x_prelim_compact.rs`, exactness argument in its header); the redo
  reads lists rebuilt as the off path holds them: the one-thread redo rebuilds each query's list when it
  reaches it and drops it after its redo, a search with more than one thread rebuilds all lists before the
  redo, query chunks hand on rebuilt lists. `LOSAT_X_BPPCOMPACT=shadow` builds both lists, compares every
  rebuilt list field by field with the off path's (it panics on a difference) and prints
  `[X_SHADOW] BPPCOMPACT compared=<lists> hsps=<n>` at exit. With `LOSAT_TIMING=1` each search of a query
  set prints `[TIMING] x_bppcompact: expanded_lists=<n>`.
- `LOSAT_X_SEEDRELEASE=1` (TBLASTX, memory; in the strict set since 2026-10-13) the tables that only the
  scan reads are freed before the HSP list of the batch's last subject is merged and linked, which is where
  the heap of a large search peaks
  (`tblastx/x_seed_release.rs`, exactness argument in its header). A query batch keeps its tables in two
  holders, each reached only through one stage per subject: the unit path of `LOSAT_X_PAIRPAR` with its group
  tables, whose stage ends when the subject's units have run, before their lists are merged (the unit path
  runs as two calls for this, `run_subject_units` and `merge_subject` in `x_pair_par.rs`, the same statements
  in the same order); and the batch's lookup table, whose stage ends after the seed stage. The end of the
  batch's last stage of a holder frees it. The worker's diagonal table, which only the reference frame loop
  reads, is not created when the unit path takes the batch without its shadow, and is otherwise dropped by
  the worker that ended the last lookup stage. With the switch off the same holders keep the tables until the
  batch has been searched, as before it existed. Measured once on AP027131 x AP027133 (genetic code 4, one thread, bytes
  requested under massif), without -> with the switch: the end of the units 110.4 -> 102.4 MiB, the merge of
  their lists 155.8 -> 81.0 MiB, the first linking call 167.6 -> 107.3 MiB, which is the peak of the run in
  both; max RSS in one run each 167 -> 110 MiB at one thread and 221 -> 163 MiB at eight. Known limits: the
  gain is for the last subject of a batch (the only one of a pairwise search), since earlier subjects are
  merged and linked while later ones still need the tables; and with `LOSAT_X_PAIRLUT=0` the batch's full
  lookup table, which the frame loop names, stays alive through the merge, so that path keeps the higher
  seed-stage peak. `LOSAT_X_STATS=1` prints, only with the switch on,
  `[X_STATS] SEEDRELEASE batches= tables_released= unit_tables_released= batches_without_worker_table= worker_tables_released=`
  (batches searched; batches whose lookup table a stage end released; batches whose unit path a stage end
  released; batches whose workers held no diagonal table; worker tables dropped at the release). There is no
  shadow switch: the switch changes no value and no order of computation, only when three allocations are
  freed, so there is no second result to compare. A read after a release is excluded by ownership instead:
  a stage is the only way to its tables and ending it consumes it (the compiler rejects a later use), the
  merge after the units stage is a function that cannot name a group table, an open stage keeps its tables
  alive, a subject that arrives after a release stops the run, and an assertion that runs with and without
  the switch stops a frame loop that would run without the worker's table
  (`LOSAT/tests/tblastx_seed_release.rs` reads the counters).

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
- `LOSAT_X_TBNBATCH=1` (needs `TBNSSIDE` and its kept translations, so four or more query
  batches, and at least one batch per thread) tblastn query batches searched side by side,
  appended in input order.
- `LOSAT_X_TBNOIDCHUNK=1` (in the strict set since 2026-10-13) tblastn with 2+ threads and several
  subjects: one entry into the search pool takes four chunks of subjects per thread instead of one
  subject per thread, and each chunk is one task that the pool gives to whichever thread is free;
  inside it the subjects of the chunk are a parallel range, so a thread that finds no chunk left
  takes part of one still being searched. A chunk is cut as NCBI's database reader cuts
  them for its search threads (`seqdbimpl.cpp:684-708` `x_FillSeqBuffer`, `seqdbatlas.hpp:496-501`
  `GetSliceSize`): consecutive subjects while their bases stay within
  `min(2^30, F) / (4 * threads) + 1`, at least one subject. NCBI's F is the largest file of the
  database; for `-subject` LOSAT takes ceil(bases of all subjects / 4), the size of the `.nsq`
  that `makeblastdb` writes for the same subjects (LOSAT's tblastn has no `-db` yet). So a subject
  longer than that is a chunk of its own and the subjects around it go to the other threads.
  Until 2026-10-14 a chunk was a hundredth of the subjects (NCBI's one-thread iterator,
  `blast_engine.c:1409`); on the fly genome (1,870 records, the seven chromosome arms first) that
  put the arms in one chunk on one thread (20 zebrafish proteins at 8 threads: 23.96 s against
  6.10 s without the switch). The subjects are searched by the same call and their lists reach
  the collector in subject order, so only the scheduling changes. It was off unless set (to any
  value) until it joined the set. Without it 28,036 subjects at 8 threads are 3,505 entries, each
  shorter than the wake-up of the workers; with it 13 (a hundredth of the subjects per chunk),
  17 (chunks bounded by residues, one per thread: about 217 subjects of 317 kb at most each) or 5
  (four such chunks per thread). An entry of one chunk per thread, the first form of the residue
  rule, lasted as long as its slowest chunk with nothing left for a thread that finished early:
  on the workstation (NUMA node 1) the zebrafish genome took 102.6 s at 8 threads against 77.8 s
  with the hundredth, and the transcriptomes and 16S 5-22 % longer, with the same CPU time. Timed with
  chunks of a hundredth on those subjects (16S rRNA genes, 40.6 Mb together) at 8 threads with other work on the
  machine: one query batch of 55 proteins 1.43-1.64 s without it and 0.76-0.85 s with it (NCBI
  BLAST+ 2.17.0 `-db` 0.68-0.78 s), one protein 0.87-1.11 s and 0.22-0.27 s; the same output and
  the same peak memory. With the chunks bounded by residues, four rounds on a loaded machine (load
  4-25; the first form, one chunk per thread): medians 1.35 s (55 proteins) and 0.61 s (one
  protein) against 1.25 s and 0.60 s with the hundredth and 2.81 s and 1.50 s without the switch;
  on the fly genome 7.0 s against 26.2 s with the hundredth and 6.8 s without the switch (NCBI
  `-db` 6.6-7.4 s; on the workstation 16.4 s against 62.1 s, NCBI 16.1 s). Nothing changes at one thread, nor where `TBNBATCH` searches the batches
  side by side (each batch on one thread). With `LOSAT_WASI_THREADS_DEBUG=1` the
  `[losat-thread-activity]` line of the subject stage ends with `pool_entries=` (the entries of
  one query batch); `LOSAT/tests/tblastn_oid_chunk.rs` reads it.
- `LOSAT_X_TBNOIDNEXT=1` tblastn with 2+ threads and several subjects, where `TBNOIDCHUNK` makes
  the chunks: no wait between entries. The threads of the pool take the next chunk when they are
  free, as NCBI's search threads ask the database for the next chunk (`seqdbimpl.cpp:333-344`),
  instead of waiting at the end of each entry for its last chunk. The thread that owns the
  search takes the lists back chunk by chunk in OID order and enters them into the collector as
  before, and searches a chunk itself while the next one in order is not finished
  (`tblastn/x_oid_flow.rs`). At most `4 * threads` chunks are taken or finished beyond the oldest
  one not taken back; a thread waits before a new chunk when that window is full, so the memory of
  the waiting lists stays bounded. The chunks and the search of each subject are those of
  `TBNOIDCHUNK`, so the report is the same. Not in the set (off unless set to any value); with
  `LOSAT_X_TBNOIDCHUNK=0` (or `LOSAT_X_DEFAULT=0` alone) there are no chunks and it changes
  nothing. With `LOSAT_WASI_THREADS_DEBUG=1` the search is one entry (`pool_entries=1`) and the
  line adds `flow_chunks=` and `flow_peak_waiting=` (the most finished chunks that waited at
  once); `LOSAT/tests/tblastn_oid_chunk.rs` reads them.
- `LOSAT_X_TBNOIDUNITS=1` tblastn with 2+ threads and several subjects, where `TBNOIDCHUNK` makes
  the chunks (with or without `TBNOIDNEXT`): a subject longer than a chunk's budget (a chunk of its
  own, such as a chromosome) runs its (frame, chunk) units on the pool, as the single subject of a
  search does with `TBNPAR`, instead of one after another on the thread that took it. Threads that
  have no chunk to search (at the end of an entry, or when no chunk is left) take its units; its
  merge in NCBI order, linking, E-values and reaping run on its thread after its last unit, so the
  report is the same. NCBI's thread searches a subject whole, so at 16 and 32 threads the longest
  chromosome bounds both programs; this goes beyond NCBI. Other subjects keep their units on their
  thread. With `TBNOIDNEXT` the threads of the flow start together (each waits until all have
  entered the flow) so that a thread waiting for a unit can only take units, never a loop of the
  flow. Not in the set (off unless set to any value). With `LOSAT_WASI_THREADS_DEBUG=1` each subject
  whose units ran on the pool prints a `prelim_units` stage line with `parallel_selected=true`;
  `LOSAT/tests/tblastn_oid_chunk.rs` counts them.
- `LOSAT_X_TBNXLATE=1` tblastn: the six translated frames of a subject are read from tables, as
  NCBI makes them from the compressed subject (`blast_util.c:1071-1090`: one table of 64 residues
  per strand, `s_BlastGetTranslationTable`, and one table read per codon,
  `BLAST_TranslateCompressedSequence`). The subject's letters become ncbi2na in one pass; a gap or
  an ambiguous letter draws from the generator seeded with the subject's length, in subject
  order, as `s_Ncbi4naToNcbi2na` does (`blast_objmgr_tools.cpp:427-475`). The reference path
  reaches the same residues through letters (the subject as ncbi4na, as ncbi2na, as `ACGT`, its
  reverse complement, then `GeneticCode::get` per codon). The tables are filled from the
  reference path's own functions, so every residue is the reference's by construction. Not in the
  set (off unless set to any value). It acts where the subject plans of the unit path are made
  (`TBNPAR`, on by default); with `LOSAT_X_DEFAULT=0` the reference loop translates each subject
  itself and the variable changes nothing. Shadow: `LOSAT_X_TBNXLATESHADOW=1` makes both
  translations for every plan and stops (panic) if a frame, a unit or an error text differs; one
  `[X_SHADOW] TBNXLATE` line on stderr at the end gives the plans and residues compared, or says
  that nothing was compared. With `LOSAT_X_STATS=1` a `[X_STATS] TBNXLATE table_plans=` line
  counts the plans made from the tables; `LOSAT/tests/tblastn_subject_translate.rs` reads both.
- `LOSAT_X_TBNSTORE=1` tblastn `-subject` without `-subject_loc`: each subject record's letters are
  packed as soon as the reader has read it (`sequence/packed_subject.rs`, the shape of a BLAST
  database volume: ncbi2na 4 bases per byte with the remainder count in the last byte, the runs of
  letters other than `ACGT` as ncbi4na codes, the lower-case intervals), so the letters of the
  whole set are never held; the records keep their deflines. A subject's letters are made again
  into a buffer when it is searched (the preliminary stage), redone (composition mode 2),
  traced back (mode 0) or reported (outfmt 0), exactly as the reader gave them (the comparison
  tests in `packed_subject.rs`). The packed bases are those NCBI's subject set-up makes
  (`encode_subject_ncbi2na_packed`, `s_Ncbi4naToNcbi2na` in `blast_objmgr_tools.cpp:427-475`), the
  form a later step can scan without the letters. Memory: about a quarter of a byte per base held
  instead of one, plus per thread the letters of the subject in hand (a thread keeps a buffer of
  at most 8 MiB between subjects). A subject range keeps the letters (the records are cut before
  the search), and so does a set with a byte that has no ncbi4na code (only an API caller can
  make one); the library entry `run_local` reads the caller's records. Not in the set (off unless
  set to any value). With `LOSAT_X_STATS=1` a `[X_STATS] TBNSTORE subjects= stored_bytes=
  letters_made=` line gives the subjects held packed, the bytes of their stores and the times
  their letters were made; `LOSAT/tests/tblastn_subject_store.rs` reads it.
- `LOSAT_X_BPGAP=1` blastp preliminary stage: the restricted gapped alignment (NCBI
  `s_RestrictedGappedAlign`, `blast_gapalign.c:991-1300`, which `s_BlastProtGappedAlignment` runs
  for a query whose best ungapped score is below 0.68 of the cutoff) keeps the two score rows and
  the column gaps of the columns that are multiples of 10 in separate arrays, and on x86_64 with
  AVX2 computes eight cells of a row at a time (`blastp/x_restricted_gap.rs`, exactness argument in
  its header): in the rows where no gap opens in A the cells depend on each other only through the
  X-drop test, whose `best_score` is a prefix maximum; in the other rows the row gap is a decayed
  prefix maximum. 32-bit values as in the reference, the matrix rows read from 16-bit tables built
  from the reference's score function. Without AVX2 (other CPUs, aarch64, Wasm) the reference
  function runs (the scalar loops on the same arrays were slower than it). A letter
  code of 28 or more, a negative gap cost, a gap cost above 65,536 or an X-drop above 2^28 leaves
  the call to the reference function. Not in the set (off unless set to any value). Shadow:
  `LOSAT_X_BPGAPSHADOW=1` runs both for every restricted alignment and stops (panic) if
  `(a_offset, b_offset, best_score)` differ; one `[X_SHADOW] BPGAP compared= not_taken=` line on
  stderr at the end. With `LOSAT_X_STATS=1` a `[X_STATS] BPGAP fast= not_taken=` line.
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
- `LOSAT_X_PAIRLUT=1` (in the strict set; acts while `LOSAT_X_PAIRPAR` is on, not with its
  shadow) TBLASTX:
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
